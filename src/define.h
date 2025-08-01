#pragma once

#include <math.h>
/*
*				define.h
*
* Global definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2016 IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SkyMaker is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SkyMaker is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SkyMaker. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		15/11/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"3.x"
#define		DATE		"2005-09-19"
#endif
/*------------------------ what, who, when and where ------------------------*/

#define         BANNER          "SkyMaker"
#define         EXECUTABLE      "sky"
#define		MYVERSION	VERSION
#define         COPYRIGHT       "2010-2012 IAP/CNRS/UPMC"
#define		DISCLAIMER	BANNER " comes with ABSOLUTELY NO WARRANTY\n" \
		"You may redistribute copies of " BANNER "\n" \
		"under the terms of the GNU General Public License."
#define		AUTHORS		"Emmanuel BERTIN <bertin@iap.fr>"
#define		WEBSITE		"http://astromatic.net/software/skymaker"
#define		INSTITUTE	"IAP  http://www.iap.fr"

/*----------------------------- Internal constants --------------------------*/
#define		OUTPUT		stdout		/* where all msgs are sent */
#define		BIG		1e+30		/* a huge number */
#define		SMALL		(1/BIG)		/* A very small number */
#define		MAXCHAR		1024		/* max. number of characters */
#ifndef M_PI
#define M_PI	3.1415926535898			/* never met before? */
#endif
#ifndef M_PI_2f
#define M_PI_2f	1.57079632679489661923f
#endif
#define C               2.9979250e8             /* speed of light in MKS */
#define	PSF_NORDER	15			/* Max size = 2^15 */

/*----------------------------- Unit conversions ----------------------------*/
#define	DEG		(M_PI/180.0)		/* one degree in rad */
#define	ARCSEC		(DEG/3600.0)		/* one arsec in rad */
#define	MICRON		1e-6			/* one micron in MKS */
#define	MM		1e-3			/* one mm in MKS */
#define	CM		1e-2			/* one cm in MKS */
#define	KM		1000.0			/* one km in MKS */
#define	PC		3.085678e16		/* one parsec in MKS */
#define	MPC		(1.0e6*PC)		/* one Mpc */

/*------------ Set defines according to machine's specificities -------------*/
#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*--------------------- in case of missing constants ------------------------*/

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*------------------------------- Other Macros -----------------------------*/

#ifdef _GNU_SOURCE
#define	DEXP(x)		exp10(x)		/* 10^x */
#define	DEXPF(x)	exp10f(x)	/* 10^x */
#else
#ifdef __USE_GNU
#define	DEXP(x)		exp(M_LN10*(x))		/* 10^x */
#define	DEXPF(x)	expf(M_LN10f*(x))	/* 10^x */
#else
#define	DEXP(x)		exp(2.30258509299*(x))		/* 10^x */
#define	DEXPF(x)	expf(2.30258509299f*(x))	/* 10^x */
#endif
#endif

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (fseek(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(pos, afile, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QFREE(x)	{free(x); x = NULL;}

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QCALLOC16(ptr, typ, nel) \
		{if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                   memset(ptr, 0, (size_t)(nel)*sizeof(typ)); \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC16(ptr, typ, nel) \
		{if (posix_memalign((void **)&ptr, 16, (size_t)(nel)*sizeof(typ))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define QPRINTF		if (prefs.verbose_type != QUIET)	fprintf
