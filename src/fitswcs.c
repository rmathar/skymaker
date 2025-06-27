/*
*				fitswcs.c
*
* Manage World Coordinate System data.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2016 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		08/03/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"fits/fitscat_defs.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"wcscelsys.h"
#include	"wcslib/wcs.h"
#include	"wcslib/wcshdr.h"

/******* create_wcs ***********************************************************
PROTO	wcsstruct *create_wcs(char **ctype, double *crval, double *crpix,
			double *cdelt, int *naxisn, int naxis)
PURPOSE	Generate a simple WCS (World Coordinate System) structure.
INPUT	Pointer to an array of char strings with WCS projection on each axis,
	pointer to an array of center coordinates (double),
	pointer to an array of device coordinates (double),
	pointer to an array of pixel scales (double),
	pointer to an array of image dimensions (int),
	number of dimensions.
OUTPUT	pointer to a WCS structure.
NOTES	If a pointer is set to null, the corresponding variables are set to
	default values.
AUTHOR	E. Bertin (IAP)
VERSION	09/08/2006
 ***/
wcsstruct	*create_wcs(char **ctype, double *crval, double *crpix,
			double *cdelt, int *naxisn, int naxis)

  {
   wcsstruct	*wcs;
   int		l;
   char *header =NULL, *inhdr = NULL;
   struct wcsprm * pihptr ;
   int pihrej, pihnwcs ;

  QCALLOC(wcs, wcsstruct, 1);
  wcs->naxis = naxis;
  QCALLOC(wcs->projp, double, naxis*100);
  wcs->nprojp = 0;

  wcs->longpole = wcs->latpole = 999.0;
  for (l=0; l<naxis; l++)
    {
    wcs->naxisn[l] = naxisn? naxisn[l] : 360.0;
/*-- The default WCS projection system is an all-sky Aitoff projection */
    if (ctype)
      strncpy(wcs->ctype[l], ctype[l], 8);
    else if (l==0)
      strncpy(wcs->ctype[l], "RA---AIT", 8);
    else if (l==1)
      strncpy(wcs->ctype[l], "DEC--AIT", 8);
    wcs->crval[l] = crval? crval[l]: 0.0;
    wcs->crpix[l] = crpix? crpix[l]: 0.0;
    wcs->cdelt[l] = 1.0;
    wcs->cd[l*(naxis+1)] = cdelt? cdelt[l] : 1.0;
    }

  wcs->epoch = wcs->equinox = 2000.0;
  QCALLOC(wcs->wcsprm, struct wcsprm, 1);

  /* represent the 5*naxis values above as string of header
  * cards to ensure that wcsprm is aligned with this
  */
  const int headlen = 5*naxis*sizeof(char[80]) ;
  header = malloc(headlen+1) ;
  header[headlen]='\0' ;
  for (l=0, inhdr=header; l<naxis; l++)
  {
	snprintf(inhdr,80, "CTYPE%d = %70s",l, wcs->ctype[l]) ;
	inhdr +=80 ;
	snprintf(inhdr,80, "CRVAL%d = %70f",l, wcs->crval[l]) ;
	inhdr +=80 ;
	snprintf(inhdr,80, "CRPIX%d = %70f",l, wcs->crpix[l]) ;
	inhdr +=80 ;
	snprintf(inhdr,80, "CDELT%d = %70f",l, wcs->cdelt[l]) ;
	inhdr +=80 ;
	snprintf(inhdr,80, "CD_%d_%d = %70f",l,l, wcs->cd[l*(naxis+1)]) ;
	inhdr +=80 ;
  }
  /* ensure all \0 are just blanks in the middle
  */
  for(int l=0 ; l < headlen ; l++)
     if ( header[l] == '\0')
        header[l] = 0x20 ;

  wcspih(header,5*naxis,WCSHDR_all,0,&pihrej,&pihnwcs,&pihptr) ;

   free(header) ;

/* Test if the WCS is recognized and a celestial pair is found */
//  wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm);
  wcsini(1,wcs->naxis,wcs->wcsprm) ;
  /* we keep only the first of these structures, although pihrej may be >1
  */
  wcssub(0,pihptr,0x0,0x0,wcs->wcsprm) ;
 
   wcsvfree(&pihnwcs,& pihptr) ;

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);  

  return wcs;
  } /* create_wcs */


/******* init_wcs ************************************************************
PROTO	void init_wcs(wcsstruct *wcs)
PURPOSE	Initialize astrometry and WCS (World Coordinate System) structures.
INPUT	WCS structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/05/2007
 ***/
void	init_wcs(wcsstruct *wcs)

  {

/* wcsprm structure */
  wcs->lng = wcs->wcsprm->lng;
  wcs->lat = wcs->wcsprm->lat;

/* Check-out chirality */
  wcs->chirality = wcs_chirality(wcs);

/* Initialize Equatorial <=> Celestial coordinate system transforms */
  init_wcscelsys(wcs);

  }


/******* init_wcscelsys *******************************************************
PROTO	void init_wcscelsys(wcsstruct *wcs)
PURPOSE	Initialize Equatorial <=> Celestial coordinate system transforms.
INPUT	WCS structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2006
 ***/
void	init_wcscelsys(wcsstruct *wcs)

  {
  double	*mat,
		a0,d0,ap,dp,ap2,y;
  int		s,lng,lat;

  lng = wcs->wcsprm->lng;
  lat = wcs->wcsprm->lat;
/* Is it a celestial system? If not, exit! */
  if (lng==lat)
    {
    wcs->celsysconvflag = 0;
    return;
    }
/* Find the celestial system */
  for (s=0; *celsysname[s][0] && strncmp(wcs->ctype[lng], celsysname[s][0], 4);
	s++);
/* Is it a known, non-equatorial system? If not, exit! */
  if (!s || !*celsysname[s][0])
    {
    wcs->celsysconvflag = 0;
    return;
    }
  wcs->celsys = (celsysenum)s;
/* Some shortcuts */
  a0 = celsysorig[s][0]*DEG;
  d0 = celsysorig[s][1]*DEG;
  ap = celsyspole[s][0]*DEG;
  dp = celsyspole[s][1]*DEG;
/* First compute in the output referential the longitude of the south pole */
  y = sin(ap - a0);
/*
  x = cos(d0)*(cos(d0)*sin(dp)*cos(ap-a0)-sin(d0)*cos(dp));
  ap2 = atan2(y,x);
*/
  ap2 = asin(cos(d0)*y) ;
/* Equatorial <=> Celestial System transformation parameters */
  mat = wcs->celsysmat;
  mat[0] = ap;
  mat[1] = ap2;
  mat[2] = cos(dp);
  mat[3] = sin(dp);

  wcs->celsysconvflag = 1;
  return;
  }


/******* read_wcs *************************************************************
PROTO	wcsstruct *read_wcs(tabstruct *tab)
PURPOSE	Read WCS (World Coordinate System) info in the FITS header.
INPUT	tab structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/11/2013
 ***/
wcsstruct	*read_wcs(tabstruct *tab)

  {
#define	FITSREADF(buf, k, val, def) \
		{if (fitsread(buf,k, &val, H_FLOAT,T_DOUBLE) != RETURN_OK) \
		   val = def; \
		}

#define	FITSREADI(buf, k, val, def) \
		{if (fitsread(buf,k, &val, H_INT,T_LONG) != RETURN_OK) \
		   val = def; \
		}

#define	FITSREADS(buf, k, str, def) \
		{if (fitsread(buf,k,str, H_STRING,T_STRING) != RETURN_OK) \
		   strcpy(str, (def)); \
		}
   char		str[MAXCHARS];
   //char		wstr1[TNX_MAXCHARS], wstr2[TNX_MAXCHARS];

   wcsstruct	*wcs;
   double	drota;
   int		j, l, naxis;
   char		name[16],
		*buf, *filename, *ptr;

  buf = tab->headbuf;
  filename = (tab->cat? tab->cat->filename : strcpy(name, "internal header"));

  FITSREADS(buf, "OBJECT  ", str, "Unnamed");

  QCALLOC(wcs, wcsstruct, 1);
  if (tab->naxis > NAXIS)
    {
    warning("Maximum number of dimensions supported by this version of the ",
	"software exceeded\n");
    tab->naxis = 2;
    }

  wcs->naxis = naxis = tab->naxis;
  QCALLOC(wcs->projp, double, naxis*100);

  for (l=0; l<naxis; l++)
    {
    wcs->naxisn[l] = tab->naxisn[l];
    sprintf(str, "CTYPE%-3d", l+1);
    FITSREADS(buf, str, str, "");
    strncpy(wcs->ctype[l], str, 8);
    sprintf(str, "CUNIT%-3d", l+1);
    FITSREADS(buf, str, str, "deg");
    strncpy(wcs->cunit[l], str, 32);
    sprintf(str, "CRVAL%-3d", l+1);
    FITSREADF(buf, str, wcs->crval[l], 0.0);
    sprintf(str, "CRPIX%-3d", l+1);
    FITSREADF(buf, str, wcs->crpix[l], 1.0);
    sprintf(str, "CDELT%-3d", l+1);
    FITSREADF(buf, str, wcs->cdelt[l], 1.0);
    sprintf(str, "CRDER%-3d", l+1);
    FITSREADF(buf, str, wcs->crder[l], 0.0);
    sprintf(str, "CSYER%-3d", l+1);
    FITSREADF(buf, str, wcs->csyer[l], 0.0);
    if (fabs(wcs->cdelt[l]) < 1e-30)
      error(EXIT_FAILURE, "*Error*: CDELT parameters out of range in ",
	filename);
    }

  if (fitsfind(buf, "CD?_????")!=RETURN_ERROR)
/*-- If CD keywords exist, use them for the linear mapping terms... */
    for (l=0; l<naxis; l++)
      for (j=0; j<naxis; j++)
        {
        sprintf(str, "CD%d_%d", l+1, j+1);
        FITSREADF(buf, str, wcs->cd[l*naxis+j], l==j?1.0:0.0)
        }
  else if (fitsfind(buf, "PC?_????")!=RETURN_ERROR)
/*-- ...If PC keywords exist, use them for the linear mapping terms... */
    for (l=0; l<naxis; l++)
      for (j=0; j<naxis; j++)
        {
        sprintf(str, "PC%d_%d", l+1, j+1);
        FITSREADF(buf, str, wcs->cd[l*naxis+j], l==j?1.0:0.0)
        wcs->cd[l*naxis+j] *= wcs->cdelt[l];
        }
  else if (fitsfind(buf, "PC0??0??")!=RETURN_ERROR)
/*-- ...If PC keywords exist, use them for the linear mapping terms... */
    for (l=0; l<naxis; l++)
      for (j=0; j<naxis; j++)
        {
        sprintf(str, "PC%03d%03d", l+1, j+1);
        FITSREADF(buf, str, wcs->cd[l*naxis+j], l==j?1.0:0.0)
        wcs->cd[l*naxis+j] *= wcs->cdelt[l];
        }
  else
    {
/*-- ...otherwise take the obsolete CROTA2 parameter */
    FITSREADF(buf, "CROTA2  ", drota, 0.0)
    wcs->cd[3] = wcs->cd[0] = cos(drota*DEG);
    wcs->cd[1] = -(wcs->cd[2] = sin(drota*DEG));
    wcs->cd[0] *= wcs->cdelt[0];
    wcs->cd[2] *= wcs->cdelt[0];
    wcs->cd[1] *= wcs->cdelt[1];
    wcs->cd[3] *= wcs->cdelt[1];
    }
  QCALLOC(wcs->wcsprm, struct wcsprm, 1);

/* Test if the WCS is recognized and a celestial pair is found */
  if (!wcsini(1,wcs->naxis,wcs->wcsprm))
    {
     char	*pstr;
     double	date;
     int	biss, dpar[3];

/*-- Coordinate reference frame */
/*-- Search for an observation date expressed in Julian days */
    FITSREADF(buf, "MJD-OBS ", date, -1.0);
    if (date<0.0)
      FITSREADF(buf, "MJDSTART", date, -1.0);
/*-- Precession date (defined from Ephemerides du Bureau des Longitudes) */
/*-- in Julian years from 2000.0 */
    if (date>0.0)
      wcs->obsdate = 2000.0 - (MJD2000 - date)/365.25;
    else
      {
/*---- Search for an observation date expressed in "civilian" format */
      FITSREADS(buf, "DATE-OBS", str, "");
      if (*str)
        {
/*------ Decode DATE-OBS format: DD/MM/YY or YYYY-MM-DD */
        for (l=0; l<3 && (pstr = strtok_r(l?NULL:str,"/- ", &ptr)); l++)
          dpar[l] = atoi(pstr);
        if (l<3 || !dpar[0] || !dpar[1] || !dpar[2])
          {
/*-------- If DATE-OBS value corrupted or incomplete, assume 2000-1-1 */
          warning("Invalid DATE-OBS value in header: ", str);
          dpar[0] = 2000; dpar[1] = 1; dpar[2] = 1;
          }
        else if (strchr(str, '/') && dpar[0]<32 && dpar[2]<100)
          {
          j = dpar[0];
          dpar[0] = dpar[2]+1900;
          dpar[2] = j;
          }

        biss = (dpar[0]%4)?0:1;
/*------ Convert date to MJD */
        date = -678956 + (365*dpar[0]+dpar[0]/4) - biss
			+ ((dpar[1]>2?((int)((dpar[1]+1)*30.6)-63+biss)
		:((dpar[1]-1)*(63+biss))/2) + dpar[2]);
        wcs->obsdate = 2000.0 - (MJD2000 - date)/365.25;
        }
      else
/*------ Well if really no date is found */
        wcs->obsdate = 0.0;
      }

    FITSREADF(buf, "EPOCH", wcs->epoch, 2000.0);
    FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->epoch);
    if (fitsread(buf, "RADESYS", str, H_STRING,T_STRING) != RETURN_OK)
      FITSREADS(buf, "RADECSYS", str,
	wcs->equinox >= 2000.0? "ICRS" : (wcs->equinox<1984.0? "FK4" : "FK5"));
    if (!strcmp(str, "ICRS"))
      wcs->radecsys = RDSYS_ICRS;
    else if (!strcmp(str, "FK5"))
      wcs->radecsys = RDSYS_FK5;
    else if (!strcmp(str, "FK4"))
      {
      if (wcs->equinox == 2000.0)
        {
        FITSREADF(buf, "EPOCH  ", wcs->equinox, 1950.0);
        FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->equinox);
        }
      wcs->radecsys = RDSYS_FK4;
      warning("FK4 precession formulae not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else if (!strcmp(str, "FK4-NO-E"))
      {
      if (wcs->equinox == 2000.0)
        {
        FITSREADF(buf, "EPOCH", wcs->equinox, 1950.0);
        FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->equinox);
        }
      wcs->radecsys = RDSYS_FK4_NO_E;
      warning("FK4 precession formulae not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else if (!strcmp(str, "GAPPT"))
      {
      wcs->radecsys = RDSYS_GAPPT;
      warning("GAPPT reference frame not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else
      {
      warning("Using ICRS instead of unknown astrometric reference frame: ",
		str);
      wcs->radecsys = RDSYS_ICRS;
      }

#if 0
/*-- Projection parameters */
    if (!strcmp(wcs->wcsprm->pcode, "TNX"))
      {
/*---- IRAF's TNX projection: decode these #$!?@#!! WAT parameters */
      if (fitsfind(buf, "WAT?????") != RETURN_ERROR)
        {
/*------ First we need to concatenate strings */
        pstr = wstr1;
        sprintf(str, "WAT1_001");
        for (j=2; fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK; j++)
	  {
          sprintf(str, "WAT1_%03d", j);
          pstr += strlen(pstr);
	  }
        pstr = wstr2;
        sprintf(str, "WAT2_001");
        for (j=2; fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK; j++)
	  {
          sprintf(str, "WAT2_%03d", j);
          pstr += strlen(pstr);
	  }
/*------ LONGPOLE defaulted to 180 deg if not found */
        if ((pstr = strstr(wstr1, "longpole"))
		|| (pstr = strstr(wstr2, "longpole")))
          pstr = strpbrk(pstr, "1234567890-+.");
        wcs->longpole = pstr? atof(pstr) : 999.0;
        wcs->latpole = 999.0;
/*------ RO defaulted to 180/PI if not found */
        if ((pstr = strstr(wstr1, "ro"))
		|| (pstr = strstr(wstr2, "ro")))
          pstr = strpbrk(pstr, "1234567890-+.");
        wcs->r0 = pstr? atof(pstr) : 0.0;
/*------ Read the remaining TNX parameters */
        if ((pstr = strstr(wstr1, "lngcor"))
		|| (pstr = strstr(wstr2, "lngcor")))
          wcs->tnx_lngcor = read_tnxaxis(pstr);
        if (!wcs->tnx_lngcor)
          error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ",
			filename);
        if ((pstr = strstr(wstr1, "latcor"))
		|| (pstr = strstr(wstr2, "latcor")))
          wcs->tnx_latcor = read_tnxaxis(pstr);
        if (!wcs->tnx_latcor)
          error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ",
			filename);
        }
      }
    else
      {
      if (fitsread(buf, "LONPOLE",&wcs->longpole,H_FLOAT,T_DOUBLE) != RETURN_OK)
        FITSREADF(buf, "LONGPOLE", wcs->longpole, 999.0);
      FITSREADF(buf, "LATPOLE ", wcs->latpole, 999.0);
/*---- Old convention */
      if (fitsfind(buf, "PROJP???") != RETURN_ERROR)
        for (j=0; j<10; j++)
          {
          sprintf(str, "PROJP%-3d", j);
          FITSREADF(buf, str, wcs->projp[j], 0.0);
          }
/*---- New convention */
      if (fitsfind(buf, "PV?_????") != RETURN_ERROR)
        for (l=0; l<naxis; l++)
          for (j=0; j<100; j++)
            {
            sprintf(str, "PV%d_%d ", l+1, j);
            FITSREADF(buf, str, wcs->projp[j+l*100], 0.0);
            }
      }
#endif /* 0 */
    }

/* Initialize other WCS structures */
  init_wcs(wcs);

/* Find the range of coordinates */
  range_wcs(wcs);

#undef FITSREADF
#undef FITSREADI
#undef FITSREADS

  return wcs;
  } /* read_wcs */


/******* end_wcs **************************************************************
PROTO	void end_wcs(wcsstruct *wcs)
PURPOSE	Free WCS (World Coordinate System) infos.
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	2025-06-26
 ***/
void	end_wcs(wcsstruct *wcs)

  {
  if (wcs)
    {
    wcsfree(wcs->wcsprm);
    free(wcs);
    }

  } /* end_wcs */


/******* range_wcs ***********************************************************
PROTO	void range_wcs(wcsstruct *wcs)
PURPOSE	Find roughly the range of WCS coordinates on all axes,
	and typical pixel scales.
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	24/08/2010
 ***/
void	range_wcs(wcsstruct *wcs)

  {
   double		step[NAXIS], raw[NAXIS], rawmin[NAXIS],
			world[NAXIS], world2[NAXIS];
   double		*worldmin, *worldmax, *scale, *worldc,
			rad, radmax, lc;
   int			linecount[NAXIS];
   int			i,j, naxis, npoints, lng,lat;

  naxis = wcs->naxis;

/* World range */
  npoints = 1;
  worldmin = wcs->wcsmin;
  worldmax = wcs->wcsmax;
/* First, find the center and use it as a reference point for lng */
  lng = wcs->lng;
  lat = wcs->lat;
  for (i=0; i<naxis; i++)
    raw[i] = (wcs->naxisn[i]+1.0)/2.0;
  if (raw_to_wcs(wcs, raw, world))
    {
/*-- Oops no mapping there! So explore the image in an increasingly large */
/*-- domain to find  a better "center" (now we know there must be angular */
/*-- coordinates) */
    for (j=0; j<100; j++)
      {
      for (i=0; i<naxis; i++)
        raw[i] += wcs->naxisn[i]/100.0*(0.5-(double)rand()/RAND_MAX);      
      if (!raw_to_wcs(wcs, raw, world))
        break;
      }
    }

  if (lng!=lat)
    lc = world[lng];
  else
    {
    lc = 0.0;   /* to avoid gcc -Wall warnings */
    lng = -1;
    }

/* Pixel scales at image center */
  scale = wcs->wcsscale;
  for (i=0; i<naxis; i++)
    {
    if ((i==lng || i==lat) && lng!=lat)
      wcs->pixscale = scale[i] = sqrt(wcs_scale(wcs, raw));
    else
      {
      raw[i] += 1.0;
      raw_to_wcs(wcs, raw, world2);
      scale[i] = fabs(world2[i] - world[i]);
      raw[i] -= 1.0;
      if (lng==lat)
        wcs->pixscale = scale[i];
      }
    wcs->wcsscalepos[i] = world[i];
    }


/* Find "World limits" */
  for (i=0; i<naxis; i++)
    {
    raw[i] = rawmin[i] = 0.5;
    step[i] = wcs->naxisn[i]/(WCS_NRANGEPOINTS-1.0);
    npoints *= WCS_NRANGEPOINTS;
    worldmax[i] = -(worldmin[i] = 1e31);
    linecount[i] = 0;
    }

  radmax = 0.0;
  worldc = wcs->wcsscalepos;

  for (j=npoints; j--;)
    {
    raw_to_wcs(wcs, raw, world);
/*-- Compute maximum distance to center */
    if ((rad=wcs_dist(wcs, world, worldc)) > radmax)
      radmax = rad;
    for (i=0; i<naxis; i++)
      {
/*---- Handle longitudes around 0 */
      if (i==lng)
        {
        world[i] -= lc;
        if (world[i]>180.0)
          world[i] -= 360.0;
        else if (world[i] <= -180.0)
          world[i] += 360.0;
        }
      if (world[i]<worldmin[i])
        worldmin[i] = world[i];
      if (world[i]>worldmax[i])
        worldmax[i] = world[i];
      }


    for (i=0; i<naxis; i++)
      {
      raw[i] += step[i];
      if (++linecount[i]<WCS_NRANGEPOINTS)
        break;
      else
        {
        linecount[i] = 0;       /* No need to initialize it to 0! */
        raw[i] = rawmin[i];
        }
      }
    }

  wcs->wcsmaxradius = radmax;

  if (lng!=lat)
    {
    worldmin[lng] = fmod_0_p360(worldmin[lng]+lc);
    worldmax[lng] = fmod_0_p360(worldmax[lng]+lc);
    if (worldmax[lat]<-90.0)
      worldmax[lat] = -90.0;
    if (worldmax[lat]>90.0)
      worldmax[lat] = 90.0;
    }

  } /* range_wcs */



/******* celsys_to_eq *********************************************************
PROTO	int celsys_to_eq(wcsstruct *wcs, double *wcspos)
PURPOSE	Convert arbitrary celestial coordinates to equatorial.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	celsys_to_eq(wcsstruct *wcs, double *wcspos)

  {
   double	*mat,
		a2,d2,sd2,cd2cp,sd,x,y;
   int		lng, lat;

  mat = wcs->celsysmat;
  a2 = wcspos[lng = wcs->wcsprm->lng]*DEG - mat[1];
  d2 = wcspos[lat = wcs->wcsprm->lat]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd2 = sin(d2);
  cd2cp = cos(d2)*mat[2];
  sd = sd2*mat[3]-cd2cp*cos(a2);
/* ...and the longitude */
  y = cd2cp*sin(a2);
  x = sd2 - sd*mat[3];
  wcspos[lng] = fmod((atan2(y,x) + mat[0])/DEG+360.0, 360.0);
  wcspos[lat] = asin(sd)/DEG;

  return RETURN_OK;
  }


/******* eq_to_celsys *********************************************************
PROTO	int eq_to_celsys(wcsstruct *wcs, double *wcspos)
PURPOSE	Convert equatorial to arbitrary celestial coordinates.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	eq_to_celsys(wcsstruct *wcs, double *wcspos)

  {
   double	*mat,
		a,d,sd2,cdcp,sd,x,y;
   int		lng, lat;

  mat = wcs->celsysmat;
  a = wcspos[lng = wcs->wcsprm->lng]*DEG - mat[0];
  d = wcspos[lat = wcs->wcsprm->lat]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd = sin(d);
  cdcp = cos(d)*mat[2];
  sd2 = sd*mat[3]+cdcp*cos(a);
/* ...and the longitude */
  y = cdcp*sin(a);
  x = sd2*mat[3]-sd;
  wcspos[lng] = fmod((atan2(y,x) + mat[1])/DEG+360.0, 360.0);
  wcspos[lat] = asin(sd2)/DEG;

  return RETURN_OK;
  }


/******* raw_to_wcs ***********************************************************
PROTO	int raw_to_wcs(wcsstruct *, double *, double *)
PURPOSE	Convert raw (pixel) coordinates to WCS (World Coordinate System).
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	2025-06-26
 ***/
int	raw_to_wcs(wcsstruct *wcs, double *pixpos, double *wcspos)

  {
   double	imgcrd[NAXIS],
		phi,theta;
   int stat[NAXIS] ;

  if (wcsp2s(wcs->wcsprm,wcs->naxis,1, pixpos, imgcrd, & phi, & theta, wcspos, stat) )
    {
    for (int i=0; i<wcs->naxis; i++)
      wcspos[i] = WCS_NOCOORD;
    return RETURN_ERROR;
    }

/* If needed, convert from a different coordinate system to equatorial */
  if (wcs->celsysconvflag)
    celsys_to_eq(wcs, wcspos);

  return RETURN_OK;
  }


/******* wcs_to_raw ***********************************************************
PROTO	int wcs_to_raw(wcsstruct *, double *, double *)
PURPOSE	Convert WCS (World Coordinate System) coords to raw (pixel) coords.
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	2025-06-26
 ***/
int	wcs_to_raw(wcsstruct *wcs, double *wcspos, double *pixpos)

  {
   double	imgcrd[NAXIS],
		phi,theta;
   int stat[NAXIS] ;

/* If needed, convert to a coordinate system different from equatorial */
  if (wcs->celsysconvflag)
    eq_to_celsys(wcs, wcspos);

  if (wcss2p(wcs->wcsprm,wcs->naxis,1,wcspos,&phi,&theta,imgcrd,pixpos,stat) )
    {
    for (int i=0; i<wcs->naxis; i++)
      pixpos[i] = WCS_NOCOORD;
    return RETURN_ERROR;
    }

  return RETURN_OK;
  }


/******* wcs_dist ***********************************************************
PROTO	double wcs_dist(wcsstruct *wcs, double *wcspos1, double *wcspos2)
PURPOSE	Compute the angular distance between 2 points on the sky.
INPUT	WCS structure,
	Pointer to the first array of world coordinates,
	Pointer to the second array of world coordinates.
OUTPUT	Angular distance (in degrees) between points.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/07/2002
 ***/
double	wcs_dist(wcsstruct *wcs, double *wcspos1, double *wcspos2)

  {
  double	d, dp;
  int		i, lng, lat;

  lng = wcs->lng;
  lat = wcs->lat;
  if (lat!=lng)
    {
/*-- We are operating in angular coordinates */
    d = sin(wcspos1[lat]*DEG)*sin(wcspos2[lat]*DEG)
	+ cos(wcspos1[lat]*DEG)*cos(wcspos2[lat]*DEG)
		*cos((wcspos1[lng]-wcspos2[lng])*DEG);
    return d>-1.0? (d<1.0 ? acos(d)/DEG : 0.0) : 180.0;
    }
  else
    {
    d = 0.0;
    for (i=0; i<wcs->naxis; i++)
      {
      dp = wcspos1[i] - wcspos2[i];
      d += dp*dp;
      }
    return sqrt(d);
    }
  }


/******* wcs_scale ***********************************************************
PROTO	double wcs_scale(wcsstruct *wcs, double *pixpos)
PURPOSE	Compute the sky area equivalent to a local pixel.
INPUT	WCS structure,
	Pointer to the array of local raw coordinates,
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2008
 ***/
double	wcs_scale(wcsstruct *wcs, double *pixpos)

  {
   double	wcspos[NAXIS], wcspos1[NAXIS], wcspos2[NAXIS], pixpos2[NAXIS];
   double	dpos1,dpos2;
   int		lng, lat;

  if (raw_to_wcs(wcs, pixpos, wcspos))
    return 0.0;

  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }

/* Compute pixel scale */
  pixpos2[lng] = pixpos[lng] + 1.0;
  pixpos2[lat] = pixpos[lat];
  if (raw_to_wcs(wcs, pixpos2, wcspos1))
    return 0.0;
  pixpos2[lng] -= 1.0;
  pixpos2[lat] += 1.0;
  if (raw_to_wcs(wcs, pixpos2, wcspos2))
    return 0.0;
  dpos1 = wcspos1[lng]-wcspos[lng];
  dpos2 = wcspos2[lng]-wcspos[lng];
  if (wcs->lng!=wcs->lat)
    {
    if (dpos1>180.0)
      dpos1 -= 360.0;
    else if (dpos1<-180.0)
      dpos1 += 360.0;
    if (dpos2>180.0)
      dpos2 -= 360.0;
    else if (dpos2<-180.0)
      dpos2 += 360.0;
    return fabs((dpos1*(wcspos2[lat]-wcspos[lat])
		-(wcspos1[lat]-wcspos[lat])*dpos2)*cos(wcspos[lat]*DEG));
    }
  else
    return fabs((dpos1*(wcspos2[lat]-wcspos[lat])
		-(wcspos1[lat]-wcspos[lat])*dpos2));
  }


/****** wcs jacobian *********************************************************
PROTO	double wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
PURPOSE	Compute the local Jacobian matrix of the astrometric deprojection.
INPUT	WCS structure,
	Pointer to the array of local raw coordinates,
	Pointer to the jacobian array (output).
OUTPUT	Determinant over spatial coordinates (=pixel area), or -1.0 if mapping
	was unsuccesful.
NOTES   Memory must have been allocated (naxis*naxis*sizeof(double)) for the
        Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION	11/10/2007
 ***/
double	wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
  {
   double	pixpos0[NAXIS], wcspos0[NAXIS], wcspos[NAXIS],
		dpos;
   int		i,j, lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    pixpos0[i] = pixpos[i];
  if (raw_to_wcs(wcs, pixpos0, wcspos0) == RETURN_ERROR)
    return -1.0;
  for (i=0; i<naxis; i++)
    {
    pixpos0[i] += 1.0;
    if (raw_to_wcs(wcs, pixpos0, wcspos) == RETURN_ERROR)
      return -1.0;
    pixpos0[i] -= 1.0;
    for (j=0; j<naxis; j++)
      {
      dpos = wcspos[j]-wcspos0[j];
      if (lng!=lat && j==lng)
        {
        if (dpos>180.0)
          dpos -= 360.0;
        else if (dpos<-180.0)
          dpos += 360.0;
        dpos *= cos(wcspos0[lat]*DEG);
        }
      jacob[j*naxis+i] = dpos;
      }
    }

  if (lng==lat)
    {
    lng = 0;
    lat = 1;
    }

  return fabs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat]
		- jacob[lat+naxis*lng]*jacob[lng+naxis*lat]);
  }


/******* wcs_chirality *******************************************************
PROTO	int wcs_chirality(wcsstruct *wcs)
PURPOSE	Compute the chirality of a WCS projection.
INPUT	WCS structure.
OUTPUT	+1 if determinant of matrix is positive, -1 if negative, 0 if null.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/09/2006
 ***/
int	wcs_chirality(wcsstruct *wcs)

  {
   double	a;
   int		lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  if (lng==lat && naxis>=2)
    {
    lng = 0;
    lat = 1;
    }

  a = wcs->cd[lng*naxis+lng]*wcs->cd[lat*naxis+lat]
	- wcs->cd[lng*naxis+lat]*wcs->cd[lat*naxis+lng];
  return a>TINY? 1 : (a<-TINY? -1 : 0);
  }


/******************************** fmod_0_p360 *******************************/
/*
Fold input angle in the [0,+360[ domain.
*/
double  fmod_0_p360(double angle)
  {
  return angle>0.0? fmod(angle,360.0) : fmod(angle,360.0)+360.0;
  }
