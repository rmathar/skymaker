#pragma once
/*
*				galaxies.h
*
* Include file for galaxies.c.
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

#include "simul.h"

#include "list.h"

/*---------------------------- Internal constants ---------------------------*/
#define	SERSIC_SMOOTHR	4.0	/* Profile smoothing radius (pixels) */
#define	VDKCUTRAD	5.0	/* van der Kruit disk truncation radius in r_h*/

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern int	make_galaxy(simstruct *sim, objstruct *obj);
extern double	raster_sersic(simstruct *sim, objstruct *obj, PIXTYPE *pix,
			int width, int height,
			double reff, double aspect, double posang, double n);

extern PIXTYPE	trunc_prof(PIXTYPE *pix, double xcenter, double ycenter,
			int width, int height);

