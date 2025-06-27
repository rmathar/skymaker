#pragma once
/*
*				fitswcs.h
*
* Include file for fitswcs.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2016 IAP/CNRS/UPMC
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

/*-------------------------------- macros -----------------------------------*/

/*----------------------------- Internal constants --------------------------*/

#define		NAXIS	2		/* Max number of FITS axes */

#define		DEG	(M_PI/180.0)	/* 1 deg in radians */
#define		ARCMIN	(DEG/60.0)	/* 1 arcsec in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MAS	(ARCSEC/1000.0)	/* 1 mas in radians */
#define		YEAR	(365.25*DAY)	/* 1 year in seconds */
#define		DAY	(24.0*HOUR)	/* 1 day in seconds */
#define		HOUR	(60.0*MINUTE)	/* 1 hour in seconds */
#define		MINUTE	60.0		/* 1 minute in seconds */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		JD2000	(2400000.5+MJD2000)	/* Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JD1950	(2400000.5+MJD1950)	/* Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		WCS_NOCOORD	1e31	/* Code for non-existing coordinates */

#define		WCS_NGRIDPOINTS	12	/* Number of WCS grid points / axis */
#define		WCS_NGRIDPOINTS2	(WCS_NGRIDPOINTS*WCS_NGRIDPOINTS)
#define		WCS_INVMAXDEG	9	/* Maximum inversion polynom degree */
#define		WCS_INVACCURACY	0.001	/* Maximum inversion error (pixels) */
#define		WCS_NRANGEPOINTS 32	/* Number of WCS range points / axis */

/*-------------------------------- typedefs ---------------------------------*/

typedef  enum {CELSYS_NATIVE, CELSYS_PIXEL, CELSYS_EQUATORIAL, CELSYS_GALACTIC,
	CELSYS_ECLIPTIC, CELSYS_SUPERGALACTIC}	celsysenum;

/*------------------------------- structures --------------------------------*/

typedef struct wcs
  {
  int		naxis;			/* Number of image axes */
  int		naxisn[NAXIS];		/* FITS NAXISx parameters */
  char		ctype[NAXIS][9];	/* FITS CTYPE strings */
  char		cunit[NAXIS][32];	/* FITS CUNIT strings */
  double	crval[NAXIS];		/* FITS CRVAL parameters */
  double	cdelt[NAXIS];		/* FITS CDELT parameters */
  double	crpix[NAXIS];		/* FITS CRPIX parameters */
  double	crder[NAXIS];		/* FITS CRDER parameters */
  double	csyer[NAXIS];		/* FITS CSYER parameters */
  double	cd[NAXIS*NAXIS];	/* FITS CD matrix */
  double	*projp;			/* FITS PV/PROJP mapping parameters */
  int		nprojp;			/* number of useful projp parameters */
  double	longpole,latpole;	/* FITS LONGPOLE and LATPOLE */
  double	wcsmin[NAXIS];		/* minimum values of WCS coords */
  double	wcsmax[NAXIS];		/* maximum values of WCS coords */
  double	wcsscale[NAXIS];	/* typical pixel scale at center */
  double	wcsscalepos[NAXIS];	/* WCS coordinates of scaling point */
  double	wcsmaxradius;		/* Maximum distance to wcsscalepos */
  int		outmin[NAXIS];		/* minimum output pixel coordinate */
  int		outmax[NAXIS];		/* maximum output pixel coordinate */
  int		lat,lng;		/* longitude and latitude axes # */
  double	r0;			/* projection "radius" */
  double	lindet;			/* Determinant of the local matrix */
  int		chirality;		/* Chirality of the CD matrix */
  double	pixscale;		/* (Local) pixel scale */
  double	ap2000,dp2000;		/* J2000 coordinates of pole */
  double	ap1950,dp1950;		/* B1950 coordinates of pole */
  double	obsdate;		/* Date of observations */
  double	equinox;		/* Equinox of observations */
  double	epoch;			/* Epoch of observations (deprec.) */
  enum {RDSYS_ICRS, RDSYS_FK5, RDSYS_FK4, RDSYS_FK4_NO_E, RDSYS_GAPPT}
		radecsys;		/* FITS RADECSYS reference frame */
  celsysenum	celsys;			/* Celestial coordinate system */
  double	celsysmat[4];		/* Equ. <=> Cel. system parameters */
  int		celsysconvflag;		/* Equ. <=> Cel. conversion needed? */
  struct wcsprm	*wcsprm;		/* WCSLIB's wcsprm structure */
  }	wcsstruct;

/*------------------------------- functions ---------------------------------*/

extern wcsstruct	*create_wcs(char **ctype, double *crval, double *crpix,
				double *cdelt, int *naxisn, int naxis),
			*read_wcs(tabstruct *tab);

extern double		fmod_0_p360(double angle),
			wcs_dist(wcsstruct *wcs,
				double *wcspos1, double *wcspos2),
			wcs_jacobian(wcsstruct *wcs, double *pixpos,
				double *jacob),
			wcs_scale(wcsstruct *wcs, double *pixpos);

extern int		celsys_to_eq(wcsstruct *wcs, double *wcspos),
			eq_to_celsys(wcsstruct *wcs, double *wcspos),
			raw_to_wcs(wcsstruct *wcs,
				double *pixpos, double *wcspos),
			wcs_chirality(wcsstruct *wcs),
			wcs_supproj(char *name),
			wcs_to_raw(wcsstruct *wcs,
				double *wcspos, double *pixpos);


extern void		end_wcs(wcsstruct *wcs),
			init_wcs(wcsstruct *wcs),
			init_wcscelsys(wcsstruct *wcs),
			range_wcs(wcsstruct *wcs);
