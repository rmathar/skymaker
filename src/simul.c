/*
*				simul.c
*
* Manage simulations.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		01/12/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "prefs.h"
#include "fitswcs.h"
#include "corr.h"
#include "imaout.h"
#include "random.h"
#include "simul.h"

/****** sim_init *************************************************************
PROTO   simstruct *sim_init(void)
PURPOSE Initialize a simulation from the prefs.
INPUT   -.
OUTPUT  Pointer to an allocated and filled sim structure.
NOTES   Global prefs variables are used.
AUTHOR  E. Bertin (IAP)
VERSION 01/12/2020
*/
simstruct	*sim_init(void)
  {
   simstruct	*sim;
   double	motfact[2], pos[2];
   int		i,o, nx,ny,xoffset,yoffset;

  QCALLOC(sim, simstruct, 1);
  strcpy(sim->filename, prefs.filename);
  strcpy(sim->headname, prefs.headname);
  strcpy(sim->inlistname, prefs.inlistname);
  strcpy(sim->outlistname, prefs.outlistname);

  sim->imatype = prefs.imatype;
  sim->imasize[0] = prefs.imasize[0];
  sim->imasize[1] = prefs.imasize[1];
  sim->imasize[2] = sim->imasize[3] = 1;
  sim->mscan[0] = prefs.mscan[0];
  sim->mscan[1] = prefs.mscan[1];
  sim->nmscan = sim->mscan[0]*sim->mscan[1];

  if (sim->imatype==GRID || sim->imatype==GRID_NONOISE)
    {
    sim->gridstep = prefs.grid_size;
    nx = sim->imasize[0]/prefs.grid_size;
    if (nx<1)
      {
      nx = 1;
      xoffset = sim->imasize[0]/2;
      }
    else
      xoffset = prefs.grid_size/2;

    ny = sim->imasize[1]/prefs.grid_size;
    if (ny<1)
      {
      ny = 1;
      yoffset = sim->imasize[1]/2;
      }
    else
      yoffset = prefs.grid_size/2;

    sim->ngrid[0] = nx;
    sim->ngrid[1] = ny;
    sim->gridoffset[0] = xoffset;
    sim->gridoffset[1] = yoffset;
    sim->gridindex = 0;
    }
  else
    /* set the flag if the configure value in LIST_COORDTYPE of 
    * the sky.conf file is not the LISTCOORD_PIXEL of src/list.h
    * Keeping the wcsflag=0 means the wcslib is never actually used!
    */
    sim->wcsflag = (prefs.listcoord_type != LISTCOORD_PIXEL);

  sim->cat = imaout_inithead(sim);

  if (sim->wcsflag && sim->wcs) {
    pos[0] = (sim->imasize[0] + 1.0) / 2.0;
    pos[1] = (sim->imasize[1] + 1.0 )/ 2.0;
    sim->pixscale[0] = sim->pixscale[1] = sqrt(wcs_scale(sim->wcs, pos)) * DEG
						/ ARCSEC;
  } else {
    sim->pixscale[0] = prefs.pixscale[0];
    sim->pixscale[1] = prefs.pixscale[1];
  }

  sim->lambdaeq = prefs.lambdaeq;
  sim->wellcap = prefs.wellcap;
  sim->satlev = prefs.satlev;
  sim->gain = prefs.gain;
  sim->ron = prefs.ron;
  sim->magzero = prefs.magzero;
  sim->expotime = prefs.expotime;
  sim->magback = prefs.magback;
  sim->psftype = prefs.psftype;
  strcpy(sim->psfname, prefs.psfname);
  sim->seeing = prefs.seeing;
  sim->psfseeingtype = prefs.psfseeingtype;
  sim->psfhalosb = prefs.psfhalosb;
  sim->psfoversamp = prefs.psfoversamp;
  sim->psfsize[0] = prefs.psfsize[0];
  sim->psfsize[1] = prefs.psfsize[1];
  sim->psfnarms = prefs.psfnarms;
  sim->psfdm1 = prefs.psfdm1;
  sim->psfdm2 = prefs.psfdm2;
  sim->psfarmw = prefs.psfarmw;
  sim->psfarmang = prefs.psfarmang;
  for (o=0; o<PSF_NVARORDER; o++)
    {
    sim->psfd80defoc[o] = prefs.psfd80defoc[o];
    sim->psfd80spher[o] = prefs.psfd80spher[o];
    sim->psfd80comax[o] = prefs.psfd80comax[o];
    sim->psfd80comay[o] = prefs.psfd80comay[o];
    sim->psfd80ast00[o] = prefs.psfd80ast00[o];
    sim->psfd80ast45[o] = prefs.psfd80ast45[o];
    sim->psfd80tri00[o] = prefs.psfd80tri00[o];
    sim->psfd80tri30[o] = prefs.psfd80tri30[o];
    sim->psfd80qua00[o] = prefs.psfd80qua22[o];
    }
  for (i=0; i<2; i++)
    {
    sim->psfdefocc[i] = prefs.psfdefocc[i];
    sim->psfspherc[i] = prefs.psfspherc[i];
    sim->psfcomac[i] = prefs.psfcomac[i];
    sim->psfastc[i] = prefs.psfastc[i];
    sim->psftric[i] = prefs.psftric[i];
    sim->psftric[i] = prefs.psftric[i];
    sim->psfquac[i] = prefs.psfquac[i];
    }
  sim->psftracktype = prefs.psftracktype;
  sim->psftrackmaj = prefs.psftrackmaj;
  sim->psftrackmin = prefs.psftrackmin;
  sim->psftrackang = prefs.psftrackang;
  sim->aurange = prefs.aurange;
  sim->scountdens = prefs.scountdens;
  sim->scountslope = prefs.scountslope;
  sim->maglim[0] = prefs.maglim[0];
  sim->maglim[1] = prefs.maglim[1];
  sim->psfmotionseed = prefs.psfmotionseed;
  sim->starposseed = prefs.starposseed;

/* Force PSF dimensions to be powers of 2 */
/*
  for (i = sim->psfsize[0]-1, j=1; i!=0; i/=2, j*=2);
  sim->psfsize[0] = j;
  for (i = sim->psfsize[1]-1, j=1; i!=0; i/=2, j*=2);
  sim->psfsize[1] = j;
*/
/* Pad the pixel map in order to bypass aliasing with the aureole FFTs */
  sim->margin[0] = sim->mscan[0]*(sim->aurange/sim->mscan[0]);
  sim->margin[1] = sim->mscan[1]*(sim->aurange/sim->mscan[1]);

/* Prepare noise correlation function */
  if (prefs.corr_type != CORR_NONE) {
    corr_generate(sim, prefs.corr_func_type, prefs.corr_scale);
    if (sim->margin[0] < sim->corrsize[0] / 2)
      sim->margin[0] = sim->corrsize[0] / 2;
    if (sim->margin[1] < sim->corrsize[1] / 2)
      sim->margin[1] = sim->corrsize[1] / 2;
  }

  sim->fimasize[0] = sim->imasize[0]+2*sim->margin[0];
  sim->fimasize[1] = sim->imasize[1]+2*sim->margin[1];

/* Obtain Fried parameter ro from seeing FWHM (Dierickx P. 1992, JMO 39, 3, 569) */
  if (sim->psfseeingtype != NO_SEEING)
    sim->ro = 0.975863*sim->lambdaeq*MICRON/(sim->seeing*ARCSEC);
  else
    sim->ro = sim->seeing = 0.0;

/*-- Compute PSF motion (Martin H.M. 1987, PASP 99, 1360)*/
  if (sim->psfseeingtype == SHORT_EXPOSURE)
    {
    sim->psfmotion = 0.41207*pow(sim->lambdaeq*MICRON/sim->psfdm1, 0.166667)
		*pow(sim->lambdaeq*MICRON/sim->ro, 0.833333)/ARCSEC;
    init_random(sim->psfmotionseed);
    QMALLOC(sim->psfmot[0], double, sim->nmscan);
    QMALLOC(sim->psfmot[1], double, sim->nmscan);
    motfact[0] = sim->psfmotion/sim->pixscale[0]/sim->mscan[0];
    motfact[1] = sim->psfmotion/sim->pixscale[1]/sim->mscan[1];
    for (i=0; i<sim->nmscan; i++)
      {
      sim->psfmot[0][i] = random_gauss(motfact[0], 0);
      sim->psfmot[1][i] = random_gauss(motfact[1], 0);
      }
    }
  else 
    {
    sim->psfmotion = 0.0;
    sim->psfmot[0] = sim->psfmot[1] = NULL;
    }

/* Compute the zero-point in e-/s/sub-image */

  if (sim->expotime && sim->gain)
    sim->magzero2 = sim->magzero
	+2.5*log10(sim->expotime*sim->gain/(sim->mscan[0]*sim->mscan[1]));
  else
    sim->magzero2 = -100.0;

/* Compute a rough estimate of the noise level in e- per pixel */
  sim->minquant = (sim->imatype==SKY_NONOISE || sim->imatype==GRID_NONOISE)?
			  QUANT_ACCURACY
			: QUANT_ACCURACY*sqrt(sim->ron*sim->ron
			  + (sim->pixscale[0]*sim->pixscale[1]
			  * DEXP(0.4*(sim->magzero-sim->magback))
				*sim->expotime*sim->gain
				/(sim->mscan[0]*sim->mscan[1])));
  if (sim->minquant<QUANT_ACCURACY)
    sim->minquant = QUANT_ACCURACY;

  return sim;
  }

/****** sim_end **************************************************************
PROTO   void sim_end(void)
PURPOSE Terminate a sim structure.
INPUT   Pointer to the sim structure.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 18/04/2017
*/
void    sim_end(simstruct *sim)

  {
  end_wcs(sim->wcs);

  free_cat(&sim->cat, 1);

  free(sim->image);
  free(sim->psfmot[0]);
  free(sim->psfmot[1]);
  free(sim->weight);
  free(sim->noise);
  free(sim->corr);
  free(sim);

  return;
  }

