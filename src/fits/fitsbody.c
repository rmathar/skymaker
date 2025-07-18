/*
*				fitsbody.c
*
* Handle memory allocation for FITS bodies.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1994,1997 ESO
*	          		(C) 1995,1996 Leiden Observatory 
*	          		(C) 1998-2021 IAP/CNRS/SorbonneU
*	          		(C) 2021-2023 CFHT/CNRS
*	          		(C) 2023-2025 CEA/AIM/UParisSaclay
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
*	Last modified:		25/03/2025
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<unistd.h>
#include	<sys/types.h>

#ifdef	HAVE_SYS_MMAN_H
#include	<sys/mman.h>
#endif
#include	"fitscat_defs.h"
#include	"fitscat.h"

#ifdef	HAVE_CFITSIO
#include    CFITSIO_H
#endif

size_t	body_maxram = BODY_DEFRAM,
	body_maxvram = BODY_DEFVRAM,
	body_ramleft, body_vramleft, body_ramflag;

int	body_vmnumber;

char	body_swapdirname[MAXCHARS] = BODY_DEFSWAPDIR;

/******* alloc_body ***********************************************************
PROTO	PIXTYPE *alloc_body(tabstruct *tab,
		void (*func)(PIXTYPE *ptr, int npix))
PURPOSE	Allocate memory for and read a FITS data body (read-only). If not
	enough RAM is available, a swap file is created.
INPUT	Table (tab) structure.
OUTPUT	Pointer to the mapped data if OK, or NULL otherwise.
NOTES	The file pointer must be positioned at the beginning of the data.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	25/03/2025
 ***/
PIXTYPE	*alloc_body(tabstruct *tab, void (*func)(PIXTYPE *ptr, int npix))
  {
   FILE		*file;
   PIXTYPE	*buffer;
   int  n;
   size_t	npix, size, sizeleft, spoonful;

  if (!body_ramflag)
    {
    body_ramleft = body_maxram;
    body_vramleft = body_maxvram;
    body_ramflag = 1;
    }

/* Return a NULL pointer if size is zero */
  if (!tab->tabsize)
    return (PIXTYPE *)NULL;

/* Check that there is a cat parent structure and that the file is open */
   if (tab->cat && !tab->cat->file)
     error(EXIT_FAILURE, "*Internal Error*: Cannot access table: ",
			tab->extname);

/* Decide if the data will go in physical memory or on swap-space */
#ifdef	HAVE_CFITSIO
  if (tab->isTileCompressed) {
    npix = (size_t)tab->naxisn[0];
    for (n=1; n<tab->naxis; n++)
      npix *= (size_t)tab->naxisn[n];
  } else
  npix = tab->tabsize/tab->bytepix;
#else
  npix = tab->tabsize/tab->bytepix;
#endif
  size = npix*sizeof(PIXTYPE);
  if (size < body_ramleft)
    {
/*-- There should be enough RAM left: try to do a malloc() */
    if ((tab->bodybuf = malloc(size)))
      {
      QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
#ifdef	HAVE_CFITSIO
      tab->cfitsio_currentElement = 1;
#endif
      read_body(tab, (PIXTYPE *)tab->bodybuf, npix);
/*---- Apply pixel processing */
      if (func)
        (*func)((PIXTYPE *)tab->bodybuf, npix);
      body_ramleft -= size;

      return (PIXTYPE *)tab->bodybuf;
      }
    else
      tab->bodybuf = NULL;
    }

  if (size < body_vramleft)
    {
/*-- Convert and copy the data to a swap file, and mmap() it */
    if (!(buffer = malloc(DATA_BUFSIZE)))
      return NULL;
    sprintf(tab->swapname, "%s/vm%05ld_%05x.tmp",
		body_swapdirname, (long)getpid(),
		(unsigned int)++body_vmnumber) ;
    if (!(file=fopen(tab->swapname, "wb+")))
      error(EXIT_FAILURE, "*Error*: cannot create swap-file ", tab->swapname);
    add_cleanupfilename(tab->swapname);
    spoonful = (size%DATA_BUFSIZE);
    if (!spoonful)
      spoonful = DATA_BUFSIZE;
    QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
#ifdef	HAVE_CFITSIO
    tab->cfitsio_currentElement = 1;
#endif
    read_body(tab, buffer, spoonful/sizeof(PIXTYPE));
/*-- Apply pixel processing */
    if (func)
      (*func)(buffer, spoonful/sizeof(PIXTYPE));
    QFWRITE(buffer, spoonful, file, tab->swapname);
    for (sizeleft = size; sizeleft -= spoonful;)
      {
      read_body(tab, buffer, (spoonful=DATA_BUFSIZE)/sizeof(PIXTYPE));
/*--- Apply pixel processing */
      if (func)
        (*func)(buffer, spoonful/sizeof(PIXTYPE));
      QFWRITE(buffer, spoonful, file, tab->swapname);
      }
    free(buffer);
    tab->bodybuf = mmap(NULL,size,PROT_READ,MAP_SHARED,fileno(file),(off_t)0);
    fclose(file);
    tab->swapflag = 1;
    body_vramleft -= size;

/*-- Memory mapping problem */
    if (tab->bodybuf == (void *)-1)
      return NULL;
    return (PIXTYPE *)tab->bodybuf;
    }

/* If no memory left at all: forget it! */
  return NULL;
  }


/******* alloc_ibody ***********************************************************
PROTO	FLAGTYPE *alloc_ibody(tabstruct *tab,
			void (*func)(FLAGTYPE *ptr, int npix))
PURPOSE	Allocate memory for and read a FITS integer data body (read-only).
	If not enough RAM is available, a swap file is created.
INPUT	Table (tab) structure.
OUTPUT	Pointer to the mapped data if OK, or NULL otherwise.
NOTES	The file pointer must be positioned at the beginning of the data.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
FLAGTYPE	*alloc_ibody(tabstruct *tab,
			void (*func)(FLAGTYPE *ptr, int npix))
  {
   FILE		*file;
   FLAGTYPE	*buffer;
   size_t	npix, size, sizeleft, spoonful;

  if (!body_ramflag)
    {
    body_ramleft = body_maxram;
    body_vramleft = body_maxvram;
    body_ramflag = 1;
    }

/* Return a NULL pointer if size is zero */
  if (!tab->tabsize)
    return (FLAGTYPE *)NULL;

/* Check that there is a cat parent structure and that the file is open */
   if (tab->cat && !tab->cat->file)
     error(EXIT_FAILURE, "*Internal Error*: Cannot access table: ",
			tab->extname);

/* Decide if the data will go in physical memory or on swap-space */
  npix = tab->tabsize/tab->bytepix;
  size = npix*sizeof(FLAGTYPE);
  if (size < body_ramleft)
    {
/*-- There should be enough RAM left: try to do a malloc() */
    if ((tab->bodybuf = malloc(size)))
      {
      QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
#ifdef	HAVE_CFITSIO
      tab->cfitsio_currentElement = 1;
#endif
      read_ibody(tab, (FLAGTYPE *)tab->bodybuf, npix);
/*---- Apply pixel processing */
      if (func)
        (*func)((FLAGTYPE *)tab->bodybuf, npix);
      body_ramleft -= size;

      return (FLAGTYPE *)tab->bodybuf;
      }
    else
      tab->bodybuf = NULL;
    }

  if (size < body_vramleft)
    {
/*-- Convert and copy the data to a swap file, and mmap() it */
    if (!(buffer = malloc(DATA_BUFSIZE)))
      return NULL;
    sprintf(tab->swapname, "%s/vm%05ld_%05x.tmp",
		body_swapdirname, (long)getpid(),
		(unsigned int)++body_vmnumber) ;
    if (!(file=fopen(tab->swapname, "wb+")))
      error(EXIT_FAILURE, "*Error*: cannot create swap-file ", tab->swapname);
    add_cleanupfilename(tab->swapname);
    spoonful = (size%DATA_BUFSIZE);
    if (!spoonful)
      spoonful = DATA_BUFSIZE;
    QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
#ifdef	HAVE_CFITSIO
    tab->cfitsio_currentElement = 1;
#endif
    read_ibody(tab, buffer, spoonful/sizeof(FLAGTYPE));
/*-- Apply pixel processing */
    if (func)
      (*func)(buffer, spoonful/sizeof(FLAGTYPE));
    QFWRITE(buffer, spoonful, file, tab->swapname);
    for (sizeleft = size; sizeleft -= spoonful;)
      {
      read_ibody(tab, buffer, (spoonful=DATA_BUFSIZE)/sizeof(FLAGTYPE));
/*--- Apply pixel processing */
      if (func)
        (*func)(buffer, spoonful/sizeof(FLAGTYPE));
      QFWRITE(buffer, spoonful, file, tab->swapname);
      }
    free(buffer);
    tab->bodybuf = mmap(NULL,size,PROT_READ,MAP_SHARED,fileno(file),(off_t)0);
    fclose(file);
    tab->swapflag = 1;
    body_vramleft -= size;

/*-- Memory mapping problem */
    if (tab->bodybuf == (void *)-1)
      return NULL;
    return (FLAGTYPE *)tab->bodybuf;
    }

/* If no memory left at all: forget it! */
  return NULL;
  }


/******* free_body ************************************************************
PROTO	void free_body(tabstruct *tab)
PURPOSE	Free FITS body data.
INPUT	Tab structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	04/03/2000
 ***/
void	free_body(tabstruct *tab)

  {
   size_t	size;

/* Free the body! (if allocated) */
  if (tab->bodybuf)
    {
    size = (tab->tabsize/tab->bytepix)*sizeof(PIXTYPE);
    if (tab->swapflag)
      {
      if (munmap(tab->bodybuf, size))
        warning("Can't unmap ", tab->cat->filename);
      tab->swapflag = 0;
      tab->bodybuf = NULL;
      body_vramleft += size;
      if (unlink(tab->swapname))
        warning("Can't delete ", tab->swapname);
      remove_cleanupfilename(tab->swapname);
      *tab->swapname = '\0';
      }
    else
      {
      QFREE(tab->bodybuf);
      body_ramleft += size;
      }
    }

/* Free the decompression buffer if allocated */
  if (tab->compress_buf)
    QFREE(tab->compress_buf);

  return;
  }

#ifdef	HAVE_CFITSIO
/******* readTileCompressed ***************************************************
 *
 * Function to read a chunk of a tile-compressed FITS image
 *
 ***/
void readTileCompressed(tabstruct *tab,  size_t	spoonful, void *bufdata0) {

   catstruct    *cat;
   int  status, hdutype;

  status = 0;
  // Exit if no parent catalog or CFITSIO information
  if (!(cat = tab->cat) || !cat->cfitsio_infptr)
    return;

  // Move to correct HDU
  fits_movabs_hdu(cat->cfitsio_infptr, tab->cfitsio_hdunum, &hdutype, &status);
  if (status != 0) {
    fprintf(stderr, "Error moving to HDU %d\n", tab->cfitsio_hdunum);
    fits_report_error(stderr, status);
  }

  // pixels count from 1
  if (!tab->cfitsio_currentElement)
    tab->cfitsio_currentElement = 1;

  // now read section of image
   int datatype;
  switch(tab->bitpix){
    case BYTE_IMG:
      datatype = TBYTE;
      break;
    case SHORT_IMG:
      datatype = TSHORT;
      break;
    case LONG_IMG:
      datatype = TLONG;
      break;
    case FLOAT_IMG:
      datatype = TFLOAT;
      break;
    case DOUBLE_IMG:
      datatype = TDOUBLE;
      break;
    default:
      datatype = TFLOAT;
      break;
  }

   int anynul;
   double bscale = 1.0, bzero = 0.0, nulval = 0.;

  // turn off any scaling so that we copy raw pixel values
  status = 0;
  fits_set_bscale(cat->cfitsio_infptr,  bscale, bzero, &status);

  // now read the image
  status = 0;
  fits_read_img(cat->cfitsio_infptr, datatype, tab->cfitsio_currentElement,
    spoonful, &nulval, bufdata0, &anynul, &status);

  // report reading error
  if (status) {
    fprintf(stderr, "CFITSIO ERROR reading start=%d end=%d absolute end=%d\n",
	tab->cfitsio_currentElement,
	(tab->cfitsio_currentElement + spoonful),
	(tab->naxisn[0]*tab->naxisn[1]));
    fits_report_error(stderr, status);
  }

  // update file 'pointer'
  tab->cfitsio_currentElement += spoonful;
}

#endif // HAVE_CFITSIO

/******* read_body ************************************************************
PROTO	read_body(tabstruct *tab, PIXTYPE *ptr, long size)
PURPOSE	Read floating point values from the body of a FITS table.
INPUT	A pointer to the tab structure,
	a pointer to the array in memory,
	the number of elements to be read.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
void	read_body(tabstruct *tab, PIXTYPE *ptr, size_t size)
  {
  catstruct		*cat;
  static double		bufdata0[DATA_BUFSIZE/sizeof(double)];
  unsigned char		cuval, cublank;
  char			*bufdata,
			cval, cblank;
  unsigned short	suval, sublank, ashort=1;
  short			val16, sval, sblank;
#ifdef HAVE_LONG_LONG_INT
  ULONGLONG		lluval, llublank;
  SLONGLONG		llval, llblank;
#endif
  unsigned int		iuval, iublank;
  int			curval, dval, blankflag, bswapflag, ival, iblank;
  
  size_t	i, bowl, spoonful, npix;
  double	bs,bz;

/* a NULL cat structure indicates that no data can be read */
  if (!(cat = tab->cat))
    return;

  bs = tab->bscale;
  bz = tab->bzero;

  blankflag = tab->blankflag;
  bswapflag = *((char *)&ashort);	// Byte-swapping flag

  switch(tab->compress_type)
    {
/*-- Uncompressed image */
    case COMPRESS_NONE:
      bowl = DATA_BUFSIZE/tab->bytepix;
      spoonful = size<bowl?size:bowl;
      for(; size>0; size -= spoonful)
        {
        if (spoonful>size)
          spoonful = size;
        bufdata = (char *)bufdata0;

#ifdef	HAVE_CFITSIO
        if (tab->isTileCompressed && cat->cfitsio_infptr)
       	  readTileCompressed(tab, spoonful, (void *)bufdata0);
        else
          QFREAD(bufdata, spoonful*tab->bytepix, cat->file, cat->filename);
#else
        QFREAD(bufdata, spoonful*tab->bytepix, cat->file, cat->filename);
#endif // HAVE_CFITSIO
        switch(tab->bitpix)
          {
          case BP_BYTE:
            if (blankflag)
	      {
              if (tab->bitsgn)
	        {
                cblank = (char)tab->blank;
#pragma ivdep
                for (i=spoonful; i--;)
                  *(ptr++) = ((cval = *(bufdata++)) == cblank)?
			-BIG : cval*bs + bz;
		}
              else
	        {
                cublank = (unsigned char)tab->blank;
#pragma ivdep
                for (i=spoonful; i--;)
                  *(ptr++) = ((cuval=*((unsigned char *)bufdata++))==cublank)?
			-BIG : cuval*bs + bz;
		}
	      }
            else
	      {
              if (tab->bitsgn)
#pragma ivdep
                for (i=spoonful; i--;)
                  *(ptr++) = *(bufdata++)*bs + bz;
              else
#pragma ivdep
                for (i=spoonful; i--;)
                  *(ptr++) = *((unsigned char *)bufdata++)*bs + bz;
	      }
            break;

          case BP_SHORT:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 2, spoonful);
            if (blankflag)
	      {
              if (tab->bitsgn)
                {
                sblank = (short)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(short))
                  *(ptr++) = ((sval = *((short *)bufdata)) == sblank)?
			-BIG : sval*bs + bz;
                }
              else
                {
                sublank = (unsigned short)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(unsigned short))
                  *(ptr++) = ((suval=*((unsigned short *)bufdata)) == sublank)?
			-BIG : suval*bs + bz;
                }
              }
            else
	      {
              if (tab->bitsgn)
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(short))
                  *(ptr++) = *((short *)bufdata)*bs + bz;
              else
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(unsigned short))
                  *(ptr++) = *((unsigned short *)bufdata)*bs + bz;
	      }
            break;

          case BP_LONG:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 4, spoonful);
            if (blankflag)
	      {
              if (tab->bitsgn)
                {
                iblank = (int)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(int))
                  *(ptr++) = ((ival = *((int *)bufdata)) == iblank)?
			-BIG : ival*bs + bz;
                }
              else
                {
                iublank = (unsigned int)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(unsigned int))
                  *(ptr++) = ((iuval = *((unsigned int *)bufdata)) == iublank)?
			-BIG : iuval*bs + bz;
                }
	      }
            else
	      {
              if (tab->bitsgn)
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(int))
                  *(ptr++) = *((int *)bufdata)*bs + bz;
              else
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(unsigned int))
                  *(ptr++) = *((unsigned int *)bufdata)*bs + bz;
	      }
            break;

#ifdef HAVE_LONG_LONG_INT
          case BP_LONGLONG:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 8, spoonful);
            if (blankflag)
	      {
              if (tab->bitsgn)
                {
                llblank = (SLONGLONG)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(SLONGLONG))
                  *(ptr++) = ((llval = *((SLONGLONG *)bufdata)) == llblank)?
			-BIG : llval*bs + bz;
                }
              else
                {
                llublank = (ULONGLONG)tab->blank;
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(ULONGLONG))
                  *(ptr++) = ((lluval = *((ULONGLONG *)bufdata)) == llublank)?
			-BIG : lluval*bs + bz;
                }
	      }
            else
	      {
              if (tab->bitsgn)
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(SLONGLONG))
                  *(ptr++) = *((SLONGLONG *)bufdata)*bs + bz;
              else
#pragma ivdep
                for (i=spoonful; i--; bufdata += sizeof(ULONGLONG))
                  *(ptr++) = *((ULONGLONG *)bufdata)*bs + bz;
	      }
            break;
#endif
          case BP_FLOAT:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 4, spoonful);
#pragma ivdep
            for (i=spoonful; i--; bufdata += sizeof(float))
              *(ptr++) = ((0x7f800000&*(unsigned int *)bufdata) == 0x7f800000)?
			-BIG : *((float *)bufdata)*bs + bz;
            break;
          case BP_DOUBLE:
            if (bswapflag)
	      {
#ifdef	HAVE_CFITSIO
              if (!tab->isTileCompressed)
                swapbytes(bufdata, 8, spoonful);
#else
              swapbytes(bufdata, 8, spoonful);
#endif
#pragma ivdep
              for (i=spoonful; i--; bufdata += sizeof(double))
                *(ptr++) = ((0x7ff00000 & *(unsigned int *)(bufdata+4))
			== 0x7ff00000)?
			-BIG : *((double *)bufdata)*bs + bz;
              }
            else
              {
#pragma ivdep
              for (i=spoonful; i--; bufdata += sizeof(double))
                *(ptr++) = ((0x7ff00000 & *(unsigned int *)bufdata)
			== 0x7ff00000)?
			-BIG : *((double *)bufdata)*bs + bz;
	      }
            break;

          default:
            error(EXIT_FAILURE,"*FATAL ERROR*: unknown BITPIX type in ",
                                "read_body()");
            break;
          }
        }
      break;

/*-- Compressed image */
    case COMPRESS_BASEBYTE:
      if (!tab->compress_buf)
        QMALLOC(tab->compress_buf, char, FBSIZE);
      bufdata = tab->compress_bufptr;
      curval = tab->compress_curval;
      npix = tab->compress_npix;
      while (size--)
        {
        if (!(npix--))
          {
          if (curval != tab->compress_checkval)
            error(EXIT_FAILURE, "*Error*: invalid BASEBYTE checksum in ",
                cat->filename);
          bufdata = tab->compress_buf;
          QFREAD(bufdata, FBSIZE, cat->file, cat->filename);
          curval = 0;
          if (bswapflag)
            swapbytes(bufdata, 4, 1);
          tab->compress_checkval = *((int *)bufdata);
         bufdata += 4;
          if (bswapflag)
            swapbytes(bufdata, 2, 1);
          npix = (int)(*((short *)bufdata))-1;
          bufdata+=2;
          }
        if ((dval=(int)*(bufdata++))==-128)
          {
          if (bswapflag)
            swapbytes(bufdata, 2, 1);
          memcpy(&val16, bufdata, 2);
          dval = (int)val16;
          if (dval==-32768)
            {
            bufdata += 2;
            if (bswapflag)
              swapbytes(bufdata, 4, 1);
            memcpy(&dval,bufdata,4);
            bufdata += 4;
            }
          else
            bufdata += 2;
          }
        *(ptr++) = dval*bs + bz;
        curval += dval;
        }
      tab->compress_curval = curval;
      tab->compress_bufptr = bufdata;
      tab->compress_npix = npix;
      break;

    case COMPRESS_PREVPIX:
      if (!tab->compress_buf)
        QMALLOC(tab->compress_buf, char, FBSIZE);
      bufdata = tab->compress_bufptr;
      curval = tab->compress_curval;
      npix = tab->compress_npix;
      while (size--)
        {
        if (!(npix--))
          {
          if (curval != tab->compress_checkval)
            error(EXIT_FAILURE, "*Error*: invalid PREV_PIX checksum in ",
                tab->cat->filename);
          bufdata = tab->compress_buf;
          QFREAD(bufdata, FBSIZE, cat->file, cat->filename);
          if (bswapflag)
            swapbytes(bufdata, 2, 3);
          curval = (int)*(short *)bufdata;
          npix = (int)*(short *)(bufdata+=2)-1;
          tab->compress_checkval = (int)(*(short *)(bufdata+=2));
          bufdata+=4;
          }
        if ((dval=(int)*(bufdata++))==-128)
          {
          if (bswapflag)
            swapbytes(bufdata, 2, 1);
          memcpy(&val16, bufdata, 2);
          curval = (int)val16;
          bufdata += 2;
          }
        else
          curval += dval;
        *(ptr++) = curval*bs + bz;
        }
      tab->compress_curval = curval;
      tab->compress_bufptr = bufdata;
      tab->compress_npix = npix;
      break;

    default:
      error(EXIT_FAILURE,"*Internal Error*: unknown compression mode in ",
                                "read_body()");
    }

  return;
  }


/******* read_ibody ***********************************************************
PROTO	read_ibody(tabstruct *tab, FLAGTYPE *ptr, long size)
PURPOSE	Read integer values from the body of a FITS table.
INPUT	A pointer to the tab structure,
	a pointer to the array in memory,
	the number of elements to be read.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	26/08/2020
 ***/
void	read_ibody(tabstruct *tab, FLAGTYPE *ptr, size_t size)
  {
   catstruct		*cat;
   static int		bufdata0[DATA_BUFSIZE/sizeof(int)];
   char			*bufdata;
   short		val16;
   unsigned short	ashort = 1;
   int			i, bowl, spoonful, npix, curval, dval, bswapflag;

/* a NULL cat structure indicates that no data can be read */
  if (!(cat = tab->cat))
    return;

  bswapflag = *((char *)&ashort);	// Byte-swapping flag

  switch(tab->compress_type)
    {
/*-- Uncompressed image */
    case COMPRESS_NONE:
      bowl = DATA_BUFSIZE/tab->bytepix;
      spoonful = size<bowl?size:bowl;
      for(; size>0; size -= spoonful)
        {
        if (spoonful>size)
          spoonful = size;
        bufdata = (char *)bufdata0;

#ifdef	HAVE_CFITSIO
        if (tab->isTileCompressed)
          readTileCompressed(tab, spoonful, (void *)bufdata0);
        else
          QFREAD(bufdata, spoonful*tab->bytepix, cat->file, cat->filename);
#else
        QFREAD(bufdata, spoonful*tab->bytepix, cat->file, cat->filename);
#endif
        switch(tab->bitpix)
          {
          case BP_BYTE:
#pragma ivdep
            for (i=spoonful; i--;)
              *(ptr++) = (FLAGTYPE)*((unsigned char *)bufdata++);
            break;

          case BP_SHORT:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 2, spoonful);
#pragma ivdep
            for (i=spoonful; i--; bufdata += sizeof(unsigned short))
              *(ptr++) = (FLAGTYPE)*((unsigned short *)bufdata);
            break;

          case BP_LONG:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 4, spoonful);
#pragma ivdep
            for (i=spoonful; i--; bufdata += sizeof(unsigned int))
              *(ptr++) = (FLAGTYPE)*((unsigned int *)bufdata);
            break;

#ifdef HAVE_LONG_LONG_INT
          case BP_LONGLONG:
#ifdef	HAVE_CFITSIO
            if (!tab->isTileCompressed && bswapflag)
#else
            if (bswapflag)
#endif
              swapbytes(bufdata, 8, spoonful);
#pragma ivdep
            for (i=spoonful; i--; bufdata += sizeof(ULONGLONG))
              *(ptr++) = (FLAGTYPE)*((ULONGLONG *)bufdata);
            break;
#endif
          case BP_FLOAT:
          case BP_DOUBLE:
            error(EXIT_FAILURE,"*Error*: expected integers, not floats, in ",
				cat->filename);
            break;
          default:
            error(EXIT_FAILURE,"*FATAL ERROR*: unknown BITPIX type in ",
				"readdata()");
            break;
          }
        }
      break;

/*-- Compressed image */
    case COMPRESS_BASEBYTE:
      if (!tab->compress_buf)
        QMALLOC(tab->compress_buf, char, FBSIZE);
      bufdata = tab->compress_bufptr;
      curval = tab->compress_curval;
      npix = tab->compress_npix;
      while (size--)
        {
        if (!(npix--))
          {
          if (curval != tab->compress_checkval)
            error(EXIT_FAILURE, "*Error*: invalid BASEBYTE checksum in ",
		cat->filename);
          bufdata = tab->compress_buf;
          QFREAD(bufdata, FBSIZE, cat->file, cat->filename);
          curval = 0;
          if (bswapflag)
            swapbytes(bufdata, 4, 1);
          tab->compress_checkval = *((int *)bufdata);
         bufdata += 4;
         if (bswapflag)
           swapbytes(bufdata, 2, 1);
          npix = (int)(*((short *)bufdata))-1;
          bufdata+=2;
          }
        if ((dval=(int)*(bufdata++))==-128)
          {
          if (bswapflag)
            swapbytes(bufdata, 2, 1);
          memcpy(&val16, bufdata, 2);
          dval = (int)val16;
          if (dval==-32768)
            {
            bufdata += 2;
            if (bswapflag)
              swapbytes(bufdata, 4, 1);
            memcpy(&dval,bufdata,4);
            bufdata += 4;
            }
          else
            bufdata += 2;
          }
        *(ptr++) = (FLAGTYPE)dval;
        curval += dval;
        }
      tab->compress_curval = curval;
      tab->compress_bufptr = bufdata;
      tab->compress_npix = npix;
      break;

    case COMPRESS_PREVPIX:
      if (!tab->compress_buf)
        QMALLOC(tab->compress_buf, char, FBSIZE);
      bufdata = tab->compress_bufptr;
      curval = tab->compress_curval;
      npix = tab->compress_npix;
      while (size--)
        {
        if (!(npix--))
          {
          if (curval != tab->compress_checkval)
            error(EXIT_FAILURE, "*Error*: invalid PREV_PIX checksum in ",
		cat->filename);
          bufdata = tab->compress_buf;
          QFREAD(bufdata, FBSIZE, cat->file, cat->filename);
          if (bswapflag)
            swapbytes(bufdata, 2, 3);
          curval = (int)*(short *)bufdata;
          npix = (int)*(short *)(bufdata+=2)-1;
          tab->compress_checkval = (int)(*(short *)(bufdata+=2));
          bufdata+=4;
          }
        if ((dval=(int)*(bufdata++))==-128)
          {
          if (bswapflag)
            swapbytes(bufdata, 2, 1);
          memcpy(&val16, bufdata, 2);
          curval = (int)val16;
          bufdata += 2;
          }
        else
          curval += dval;
        *(ptr++) = (FLAGTYPE)curval;
        }
      tab->compress_curval = curval;
      tab->compress_bufptr = bufdata;
      tab->compress_npix = npix;
      break;

    default:
      error(EXIT_FAILURE,"*Internal Error*: unknown compression mode in ",
				"readdata()");
    }

  return;
  }



/******* write_ibody ***********************************************************
PROTO	write_ibody(tabstruct *tab, FLAGTYPE *ptr, long size)
PURPOSE	Write integer values to a FITS body.
INPUT	A pointer to the tab structure,
	a pointer to the array in memory,
	the number of elements to be written.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	11/02/2020
 ***/
void	write_ibody(tabstruct *tab, FLAGTYPE *ptr, size_t size)
  {
   static FLAGTYPE	bufdata0[DATA_BUFSIZE/sizeof(FLAGTYPE)];
   catstruct		*cat;
   char			*cbufdata0;
   size_t		i, bowl, spoonful;
   unsigned short	ashort = 1;
   double		bs,bz;
   int			bswapflag;

  bswapflag = *((char *)&ashort);	// Byte-swapping flag

  bs = tab->bscale;
  bz = tab->bzero;
  cat = tab->cat;
  if (!cat)
    error(EXIT_FAILURE, "*Internal Error*: no parent cat structure for table ",
		tab->extname);

  cbufdata0 = (char *)bufdata0;	/* A trick to remove gcc aliasing warnings */
 
  switch(tab->compress_type)
    {
/*-- Uncompressed image */
    case COMPRESS_NONE:
      bowl = DATA_BUFSIZE/tab->bytepix;
      spoonful = size<bowl?size:bowl;
      for(; size>0; size -= spoonful)
        {
        if (spoonful>size)
          spoonful = size;
        switch(tab->bitpix)
          {
          case BP_BYTE:
            if (tab->bitsgn)
              {
               char	*bufdata = (char *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (char)*(ptr++);
              }
            else
              {
               unsigned char	*bufdata = (unsigned char *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (unsigned char)*(ptr++);
              }
            break;

          case BP_SHORT:
            if (tab->bitsgn)
              {
               short	*bufdata = (short *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (short)*(ptr++);
              }
            else
              {
               unsigned short	*bufdata = (unsigned short *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (unsigned short)*(ptr++);
              }
            if (bswapflag)
              swapbytes(cbufdata0, 2, spoonful);
            break;

          case BP_LONG:
           if (tab->bitsgn)
              {
               int	*bufdata = (int *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (int)*(ptr++);
              }
            else
              {
               unsigned int	*bufdata = (unsigned int *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (unsigned int)*(ptr++);
              }
            if (bswapflag)
              swapbytes(cbufdata0, 4, spoonful);
            break;

#ifdef HAVE_LONG_LONG_INT
          case BP_LONGLONG:
           if (tab->bitsgn)
              {
               SLONGLONG	*bufdata = (SLONGLONG *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (SLONGLONG)*(ptr++);
              }
            else
              {
               ULONGLONG	*bufdata = (ULONGLONG *)cbufdata0;
#pragma ivdep
              for (i=spoonful; i--;)
                *(bufdata++) = (ULONGLONG)*(ptr++);
              }
            if (bswapflag)
              swapbytes(cbufdata0, 8, spoonful);
            break;
#endif
          case BP_FLOAT:
            {
             float	*bufdata = (float *)cbufdata0;
#pragma ivdep
            for (i=spoonful; i--;)
              *(bufdata++) = (float)((double)*(ptr++)-bz)/bs;
            if (bswapflag)
              swapbytes(cbufdata0, 4, spoonful);
            }
            break;

          case BP_DOUBLE:
            {
             double	*bufdata = (double *)cbufdata0;
#pragma ivdep
            for (i=spoonful; i--;)
              *(bufdata++) = ((double)*(ptr++)-bz)/bs;
            if (bswapflag)
              swapbytes(cbufdata0, 8, spoonful);
            }
            break;

          default:
            error(EXIT_FAILURE,"*FATAL ERROR*: unknown BITPIX type in ",
                                "read_body()");
            break;
          }
        QFWRITE(cbufdata0, spoonful*tab->bytepix, cat->file, cat->filename);
        }
      break;

/*-- Compressed image */
    case COMPRESS_BASEBYTE:
      break;

    case COMPRESS_PREVPIX:
      break;

    default:
      error(EXIT_FAILURE,"*Internal Error*: unknown compression mode in ",
                                "read_body()");
    }

  }
