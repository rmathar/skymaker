/*
*				random.c
*
* Random number generators.
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

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "define.h"
#include "globals.h"
#include "random.h"

#ifndef	THREADS_NMAX
#define	THREAD_NMAX 16
#endif

#define	LOGGAMMA(x)	lgamma_r(x,&signp)	/* Generally slower */

#ifdef USE_THREADS
static unsigned int	glob_seed[THREADS_NMAX];
#endif

/****** random_gauss ********************************************************
PROTO   double random_gauss(double sigma)
PURPOSE Generate a random number with a (centered) Gaussian pdf.
INPUT   Standard deviation.
OUTPUT  Gauss-distributed random number.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
double	random_gauss(double sigma, int p)
  {
   double	x,y,z, r;

  while((z=pow(x=random_double(p)-0.5,2.0) + pow(y=random_double(p)-0.5,2.0))
	> 0.25);
  while ((r=random_double(p)) <= 0.0);

  return sigma*sqrt(-2.0*log(r)/z)*x;
  }


/****** random_int **********************************************************
PROTO   int random_int(void)
PURPOSE Generate a random integer over the range [0,RAND_MAX].
INPUT   Process index (ignored for single threads).
OUTPUT  Random integer number with uniform distribution.
NOTES   The actual upper bound of the range is implementation-dependent.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
int	random_int(int p)
  {
#ifdef USE_THREADS
  return (int)rand_r(&glob_seed[p]);
#else
  return (int)rand();
#endif
  }


/****** random_double ********************************************************
PROTO   double random_double(void)
PURPOSE Generate a random number with uniform distribution over [0.0,1.0[
INPUT   Process index (ignored for single threads).
OUTPUT  Random double with uniform distribution.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
double	random_double(int p)
  {
#ifdef USE_THREADS
  return (double)rand_r(&glob_seed[p]) / RAND_MAX;
#else
  return (double)rand() / RAND_MAX;
#endif
  }


/****** init_random **********************************************************
PROTO	void init_random(int seed)
PURPOSE	Initialize the random number generator.
INPUT	Seed.
OUTPUT	-.
NOTES	The seed is used to initialize the random sequence at a particular
	position, which is implementation-dependent. If seed = 0, then the
	actual seed is taken from the time() function (which varies each
	second). For multithreaded code, a series of secondary seeds is
	generated.
AUTHOR	E. Bertin (IAP)
VERSION	07/03/2018
*/
void	init_random(int seed)
  {
   struct timeval	t; 
#ifdef USE_THREADS
   int			p;
#endif

  if (seed)
    srand((unsigned int)seed);
  else {
    gettimeofday(&t, NULL);
    srand((unsigned int)(t.tv_usec));
  }
#ifdef USE_THREADS
  for (p=0; p<THREADS_NMAX; p++)
    glob_seed[p] = (unsigned int)rand();
#endif

  return;
  }


/****** random_poisson *******************************************************
PROTO   double random_poisson(double xm, int p)
PURPOSE Returns a random number with Poisson deviate (based on algorithm
	described in Numerical Recipes in C, chap. 7) centered on xm.
INPUT   Mean of the Poisson distribution,
	process index (ignored for single threads).
OUTPUT  A double containing the integer (!) variable with Poisson deviate.
NOTES   I am still searching for a faster algorithm!!
AUTHOR  E. Bertin (IAP)
VERSION 05/04/2013
*/
double	random_poisson(double xm, int p)
  {
   double		sq,alxm,g,oldm,em,t,y;
   int			signp;

  sq = alxm = g = 0.0;
  oldm = -1.0;
  if (xm < 12.0)
    {
    if (xm != oldm)
      {
      oldm=xm;
      g=exp(-xm);
      }
    em = -1.0;
    t=1.0;
    do
      {
      em += 1.0;
      t *= random_double(p);
      } while (t > g);
    }
  else
    {
    if (xm != oldm)
      {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-LOGGAMMA(xm+1.0);
      }
    do
      {
      do
        {
        y=tan(M_PI*random_double(p));
        em=sq*y+xm;
        } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-LOGGAMMA(em+1.0)-g);
      } while (random_double(p) > t);
    }

  return em;
  }

