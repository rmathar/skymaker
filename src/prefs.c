/*
*				prefs.c
*
* Run-time configuration functions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		13/11/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(USE_THREADS) \
&& (defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD))	/* BSD, Apple */
 #include <sys/types.h>
 #include <sys/sysctl.h>
#elif defined(USE_THREADS) && defined(HAVE_MPCTL)		/* HP/UX */
 #include <sys/mpctl.h>
#endif

#define USE_CPUREDISPATCH 1

#if defined(__INTEL_COMPILER) && defined (USE_CPUREDISPATCH)
 #include <cpuid.h>
#endif

#ifdef HAVE_MKL
 #include MKL_H
#endif

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "preflist.h"

prefstruct	prefs;

/********************************* dumpprefs ********************************/
/*
Print the default preference parameters.
*/
void	dumpprefs(int state)
  {
   char **dp;

  dp = default_prefs;
  while (**dp)
    if (**dp != '*')
      printf("%s\n",*(dp++));
    else if (state)
      printf("%s\n",*(dp++)+1);
    else
      dp++;

  return;
  }


/********************************* readprefs ********************************/
/*
Read a configuration file in ``standard'' format (see the SExtractor
documentation)
*/
void    readprefs(char *filename, char **argkey, char **argval, int narg)

  {
   FILE		*infile;
   char		*cp, str[MAXCHAR], *keyword, *value, **dp;
   int		i, ival, nkey, warn, argi, flagc, flagd, flage, flagz;
   double	dval;
#ifndef	NO_ENVVAR
   static char	value2[MAXCHAR],envname[MAXCHAR];
   char		*dolpos;
#endif


  if ((infile = fopen(filename,"r")) == NULL)
    {
    flage = 1;
    warning(filename, " not found, using internal defaults");
    }
  else
    flage = 0;

/*Build the keyword-list from pkeystruct-array */

  for (i=0; key[i].name[0]; i++)
    strcpy(keylist[i], key[i].name);
  keylist[i][0] = '\0';

/*Scan the configuration file*/

  argi=0;
  flagc = 0;
  flagd = 1;
  dp = default_prefs;
  for (warn=0;;)
    {
    if (flagd)
      {
      if (**dp)
        {
        if (**dp=='*')
          strcpy(str, *(dp++)+1);
        else
          strcpy(str, *(dp++));
        }
      else
        flagd = 0;
      }
    if (!flagc && !flagd)
      if (flage || !fgets(str, MAXCHAR, infile))
        flagc=1;

    if (flagc)
      {
      if (argi<narg)
        {
        sprintf(str, "%s %s", argkey[argi], argval[argi]);
        argi++;
        }
      else
        break;
      }
    keyword = strtok(str, notokstr);
    if (keyword && keyword[0]!=0 && keyword[0]!=(char)'#')
      {
     if (warn>=10)
        error(EXIT_FAILURE, "*Error*: No valid keyword found in ", filename);
      nkey = findkeys(keyword, keylist, FIND_STRICT);
      if (nkey!=RETURN_ERROR)
        {
        value = strtok((char *)NULL, notokstr);
#ifndef	NO_ENVVAR
/*------ Expansion of environment variables (preceded by '$') */
        if (value && (dolpos=strchr(value, '$')))
          {
           int	nc;
           char	*valuet,*value2t, *envval;

          value2t = value2;
          valuet = value;
          while (dolpos)
            {
            while (valuet<dolpos)
              *(value2t++) = *(valuet++);	/* verbatim copy before '$' */
            if (*(++valuet) == (char)'{')
              valuet++;
            strncpy(envname, valuet, nc=strcspn(valuet,"}/:\"\'\\"));
            *(envname+nc) = (char)'\0';
            if (*(valuet+=nc) == (char)'}')
              valuet++;
            if (!(envval=getenv(envname)))
              error(EXIT_FAILURE, "Environment variable not found: ",
				envname);
            while(*envval)			/* Copy the ENV content */
              *(value2t++) = *(envval++);
            while(*valuet && *valuet!=(char)'$')/* Continue verbatim copy */
              *(value2t++) = *(valuet++);
            if (*valuet)
              dolpos = valuet;
            else
              {
              dolpos = NULL;
              *value2t = (char)'\0';
              }
	    }

          value = value2;
          }
#endif
        switch(key[nkey].type)
          {
          case P_FLOAT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            dval = atof(value);
            if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
              *(double *)(key[nkey].ptr) = dval;
            else
              error(EXIT_FAILURE, keyword," keyword out of range");
            break;

          case P_INT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            ival = atoi(value);
            if (ival>=key[nkey].imin && ival<=key[nkey].imax)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " keyword out of range");
            break;

          case P_STRING:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," string is empty!");
            strcpy((char *)key[nkey].ptr, value);
            break;

          case P_BOOL:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if ((cp = strchr("yYnN", (int)value[0])))
              *(int *)(key[nkey].ptr) = (tolower((int)*cp)=='y')?1:0;
            else
              error(EXIT_FAILURE, keyword, " value must be Y or N");
            break;

          case P_KEY:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if ((ival = findkeys(value, key[nkey].keylist,FIND_STRICT))
			!= RETURN_ERROR)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " set to an unknown keyword");
            break;

          case P_BOOLLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              if ((cp = strchr("yYnN", (int)value[0])))
                ((int *)(key[nkey].ptr))[i] = (tolower((int)*cp)=='y')?1:0;
              else
                error(EXIT_FAILURE, keyword, " value must be Y or N");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_INTLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              ival = strtol(value, NULL, 0);
              if (ival>=key[nkey].imin && ival<=key[nkey].imax)
                ((int *)key[nkey].ptr)[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_FLOATLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              dval = atof(value);
              if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
                ((double *)key[nkey].ptr)[i] = dval;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_KEYLIST:
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
	      if ((ival = findkeys(value, key[nkey].keylist, FIND_STRICT))
			!= RETURN_ERROR)
                ((int *)(key[nkey].ptr))[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " set to an unknown keyword");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_STRINGLIST:
            if (!value || value[0]==(char)'#')
              {
              value = "";
              flagz = 1;
              }
            else
              flagz = 0;
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              free(((char **)key[nkey].ptr)[i]);
              QMALLOC(((char **)key[nkey].ptr)[i], char, MAXCHAR);
              strcpy(((char **)key[nkey].ptr)[i], value);
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = flagz?0:i;
            break;

          default:
            error(EXIT_FAILURE, "*Internal ERROR*: Type Unknown",
				" in readprefs()");
            break;
          }
        key[nkey].flag = 1;
        }
      else
        {
        warning(keyword, " keyword unknown");
        warn++;
        }
      }
    }

  for (i=0; key[i].name[0]; i++)
    if (!key[i].flag)
      error(EXIT_FAILURE, key[i].name, " configuration keyword missing");
  if (!flage)
    fclose(infile);

  return;
  }


/********************************* findkey **********************************/
/*
find an item within a list of keywords.
*/
int	findkeys(char *str, char keyw[][32], int mode)

  {
  int i;

  for (i=0; keyw[i][0]; i++)
    if (!cistrcmp(str, keyw[i], mode))
      return i;

  return RETURN_ERROR;
  }


/******************************* cistrcmp ***********************************/
/*
case-insensitive strcmp.
*/
int	cistrcmp(char *cs, char *ct, int mode)

  {
   int  i, diff;

  if (mode)
    {
    for (i=0; cs[i]&&ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }
  else
    {
    for (i=0; cs[i]||ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }

  return 0;
  }



/********************************* useprefs *********************************/
/*
Update various structures according to the prefs.
*/
void	useprefs()

  {
   char		str[MAXCHAR], *end;
#if defined(__INTEL_COMPILER) && defined (USE_CPUREDISPATCH)
   extern int	__intel_cpu_indicator;
   unsigned int	eax, ebx, ecx, edx;
#endif
#ifdef USE_THREADS
   int                  nproc;
#endif

/* Multithreading */
#ifdef USE_THREADS
  if (prefs.nthreads<=0)
    {
/*-- Get the number of processors for parallel builds */
/*-- See, e.g. http://ndevilla.free.fr/threads */
    nproc = -1;
#if defined(_SC_NPROCESSORS_ONLN)               /* AIX, Solaris, Linux */
    nproc = (int)sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    nproc = (int)sysconf(_SC_NPROCESSORS_CONF);
#elif defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD) /* BSD, Apple */
    {
     int        mib[2];
     size_t     len;

     mib[0] = CTL_HW;
     mib[1] = HW_NCPU;
     len = sizeof(nproc);
     sysctl(mib, 2, &nproc, &len, NULL, 0);
     }
#elif defined (_SC_NPROC_ONLN)                  /* SGI IRIX */
    nproc = sysconf(_SC_NPROC_ONLN);
#elif defined(HAVE_MPCTL)                       /* HP/UX */
    nproc =  mpctl(MPC_GETNUMSPUS_SYS, 0, 0);
#endif

    if (nproc>0)
       prefs.nthreads = ((prefs.nthreads) && nproc>(-prefs.nthreads))?
			-prefs.nthreads
			: nproc;
    else
      {
      prefs.nthreads = prefs.nthreads? -prefs.nthreads : 2;
      sprintf(str, "NTHREADS defaulted to %d", prefs.nthreads);
      warning("Cannot find the number of CPUs on this system:", str);
      }
    }
#if defined(HAVE_MKL)
  mkl_set_num_threads(prefs.nthreads);
#elif !defined(HAVE_FFTWF_MP)
   if (prefs.nthreads>1)
     warning("This executable has been compiled using a version of the FFTW "
	"library without support for multithreading. ",
	"Performance will be degraded.");
#endif

#else
  if (prefs.nthreads != 1)
    {
    prefs.nthreads = 1;
    warning("NTHREADS != 1 ignored: ",
	"this build of " BANNER " is single-threaded");
    }
#endif

/* Override INTEL CPU detection routine to help performance on 3rd-party CPUs */
#if defined(__INTEL_COMPILER) && defined(USE_CPUREDISPATCH)
  __get_cpuid(1, &eax, &ebx, &ecx, &edx);
  if (ebx&bit_AVX2)
    __intel_cpu_indicator = 0x400000;
  else if (ecx&bit_AVX)
    __intel_cpu_indicator = 0x20000;
  else if (ecx&bit_SSE4_2)
    __intel_cpu_indicator = 0x8000;
  else if (ecx&bit_SSE4_1)
    __intel_cpu_indicator = 0x2000;
  else if (ecx&bit_SSSE3)
    __intel_cpu_indicator = 0x1000;
  else if (ecx&bit_SSE3)
    __intel_cpu_indicator = 0x0800;
  else if (edx&bit_SSE2)
    __intel_cpu_indicator = 0x0200;
  else if (edx&bit_SSE)
    __intel_cpu_indicator = 0x0080;
  else if (edx&bit_MMX)
    __intel_cpu_indicator = 0x0008;
  else
    __intel_cpu_indicator = 0x0001;
#endif

/* Get rid of the original extension in the image name, if present */
  strcpy(prefs.outlistname, prefs.filename);
  if ((end=strrchr(prefs.outlistname, '.')))
    *end = 0;
  strcat(prefs.outlistname, ".list");

/* Input and output lists must have different filenames */
  if (!strcmp(prefs.outlistname, prefs.inlistname))
    error(EXIT_FAILURE, "*Error*: input and output lists would have the "
	"same filename: ", prefs.outlistname);

  if (prefs.nimasize < 2)
    prefs.imasize[1] = prefs.imasize[0];
  if (prefs.npixscale < 2)
    prefs.pixscale[1] = prefs.pixscale[0];
  if (prefs.nmscan < 2)
    prefs.mscan[1] = prefs.mscan[0];
  if (prefs.npsfsize < 2)
    prefs.psfsize[1] = prefs.psfsize[0];

/* Ensure the image dimensions are multiples of the mscan step number */
  if (prefs.imasize[0]%prefs.mscan[0])
    {
    sprintf(str, "%d",prefs.mscan[0]);
    error(EXIT_FAILURE, "*Error*: IMAGE_SIZE[0] must be a multiple of ", str);
    }
  if (prefs.imasize[1]%prefs.mscan[1])
    {
    sprintf(str, "%d",prefs.mscan[1]);
    error(EXIT_FAILURE, "*Error*: IMAGE_SIZE[1] must be a multiple of ", str);
    }

  return;
  }
