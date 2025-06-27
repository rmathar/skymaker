/*
*				fitswrite.c
*
* Low-level functions for writing LDAC FITS catalogs.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1995-2024 CEA/AIM/UParisSaclay
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
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef	HAVE_CONFIG_H
#include "config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"


/****** save_cat **************************************************************
PROTO	void save_cat(catstruct *cat, char *filename)
PURPOSE	Save a FITS catalog with name filename.
INPUT	catalog structure,
	filename.
OUTPUT	-.
NOTES	Any preexisting file with name filename is overwritten.
AUTHOR	E. Bertin (IAP)
VERSION	09/09/2003
 ***/
void	save_cat(catstruct *cat, char *filename, struct wcsprm *wcs)

  {
   tabstruct	*tab;
   int		i;

  strcpy(cat->filename, filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);

  tab = cat->tab;
/*Go through each segment in the right order to save data*/
  for (i=0; i<cat->ntab; i++)
    {
/*-- Make sure that the tab header is primary or extension! */
    if (i)
      ext_head(tab);
    else
      prim_head(tab);
    save_tab(cat, tab);
    while (!((tab=tab->nexttab)->nseg))
       ;
    }

  if (close_cat(cat) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Problem while closing", cat->filename);

  } /* save_cat */


/****** save_tab **************************************************************
PROTO	void save_tab(catstruct *cat, tabstruct *tab)
PURPOSE	Save a FITS table.
INPUT	pointer to the catalog structure,
	pointer to the table structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/02/2020
 ***/
void	save_tab(catstruct *cat, tabstruct *tab)

  {
   catstruct		*tabcat;
   keystruct		*key;
   tabstruct		*keytab;
   KINGSIZE_T		tabsize;
   long			size;
   char			*buf, *inbuf, *outbuf, *fptr,*ptr;
   unsigned short	ashort = 1;
   int			b,j,k,o, nbytes,nkey,nobj,spoonful,
			tabflag, larrayin,larrayout, esize, bswapflag;

  bswapflag = *((char *)&ashort);	// Byte-swapping flag
/*  The header itself*/
  tabflag = save_head(cat, tab)==RETURN_OK?1:0;
/*  Allocate memory for the output buffer */
  tabsize = 0;
  tabcat = NULL;	/* to satisfy gcc -Wall */
  inbuf = NULL;		/* to satisfy gcc -Wall */
  if (tabflag)
    {
/*-- If segment is a binary table, save it row by row */
    QMALLOC(outbuf, char, (larrayout = tab->naxisn[0]));
    nkey = tab->nkey;
    tabsize = larrayin = 0;
    for (j=tab->nseg; j--;)
      {
      update_tab(tab);
/*---- Scan keys to find the reference tab and other things*/
      keytab = NULL;
      key = tab->key;
      for (k=nkey; k--; key = key->nextkey)
        if (!key->ptr)
          {
          keytab = key->tab;
          tabcat = keytab->cat;
          }
/*---- If table contains some keys with no ptrs, we have to access a file */
      if (keytab)
        {
        QMALLOC(inbuf, char, (larrayin = keytab->naxisn[0]));
        if (open_cat(tabcat, READ_ONLY) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Cannot access ", tabcat->filename);
        QFSEEK(tabcat->file, keytab->bodypos, SEEK_SET, tabcat->filename);
        }
      nobj = tab->naxisn[1];
      for (o=0; o<nobj; o++)
        {
        if (keytab)
          QFREAD(inbuf, larrayin, tabcat->file, tabcat->filename);
        fptr = outbuf;
        for (k=nkey; k--; key = key->nextkey)
          {
          nbytes = key->nbytes;
          ptr = key->ptr? (char *)key->ptr+nbytes*o:inbuf+key->pos;
          for (b=nbytes; b--;)
            *(fptr++) = *(ptr++);
          if (bswapflag)
            if (key->ptr)
              {
              esize = t_size[key->ttype];
              swapbytes(fptr-nbytes, esize, nbytes/esize);
              }
          }
        QFWRITE(outbuf, larrayout, cat->file, cat->filename);
        }
      if (keytab)
        {
        free(inbuf);
        if (close_cat(tabcat) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Problem while closing",
		tabcat->filename);
        }
      tabsize += tab->tabsize;
      tab = tab->nexttab;
      }
    free(outbuf);
    }
  else
    {
/*-- If segment is not a binary table, save it ``as it is'' */
/*-- We use a limited-size buffer ``in case of'' */
    size = tabsize = tab->tabsize;
    if (tabsize)
      {
      if (tab->bodybuf)
        {
/*------ A body is present in memory and needs to be written */
        if (bswapflag)
          swapbytes(tab->bodybuf, tab->bytepix, tabsize/tab->bytepix);
        QFWRITE(tab->bodybuf, (size_t)tabsize, cat->file, cat->filename);
        if (bswapflag)
          swapbytes(tab->bodybuf, tab->bytepix, tabsize/tab->bytepix);
        }
      else
/*------ The body should be copied from the source tab */
        {
        tabcat = tab->cat;
        spoonful = size<DATA_BUFSIZE?size:DATA_BUFSIZE;
        QMALLOC(buf, char, spoonful);
        if (open_cat(tabcat, READ_ONLY) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Cannot access ", tabcat->filename);
        QFSEEK(tabcat->file, tab->bodypos, SEEK_SET, tabcat->filename);
        for (;size>0; size -= spoonful)
          {
          if (spoonful>size)
          spoonful = size;
          QFREAD(buf, spoonful, tabcat->file, tabcat->filename);
          QFWRITE(buf, spoonful, cat->file, cat->filename);
          }
        free(buf);
        if (close_cat(tabcat) != RETURN_OK)
          error(EXIT_FAILURE, "*Error*: Problem while closing",
		tabcat->filename);
        }
      }
    }

/* FITS padding*/
  pad_tab(cat, tabsize);

  } /* save_tab */


/****** save_head *************************************************************
PROTO	int save_head(catstruct *cat, tabstruct *tab)
PURPOSE	Save a FITS table header.
INPUT	catalog structure,
	table structure.
OUTPUT	RETURN_OK if tab is a binary table, or RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/12/2004
 ***/
int	save_head(catstruct *cat, tabstruct *tab)

  {
   int		tabflag;


/*  Make the table parameters reflect its content*/
  update_tab(tab);
/*  The header itself*/
  tabflag = update_head(tab);
  QFTELL(cat->file, tab->headpos, cat->filename);
  QFWRITE(tab->headbuf, tab->headnblock*FBSIZE, cat->file, cat->filename);

  return tabflag;
  }


/******* pad_tab *************************************************************
PROTO	int pad_tab(catstruct *cat, KINGSIZE_T size)
PURPOSE	Pad the FITS body of a tab with 0's to FBSIZE.
INPUT	A pointer to the cat structure,
	the number of elements that have been written.
OUTPUT	RETURN_OK if padding necessary, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	23/01/2003
 ***/
int pad_tab(catstruct *cat, KINGSIZE_T size)
  {
   static char  padbuf[FBSIZE];
   int		padsize;

  padsize = PADEXTRA(size);
  if (padsize)
    {
    QFWRITE(padbuf, padsize, cat->file, cat->filename);
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }


