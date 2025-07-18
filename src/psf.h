#pragma once
/*
*				psf.h
*
* Include file for psf.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		12/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "simul.h"

/*---------------------------- Internal constants ---------------------------*/
/*-------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	center_psf(simstruct *sim),
		freepsf(simstruct *sim),
		makeaureole(simstruct *sim),
		makepsf(simstruct *sim),
		readpsf(simstruct *sim);

extern PIXTYPE	*interp_psf(simstruct *sim, double *pos, double *dpos),
		*interp_dft(simstruct *sim, int order, double *pos,
			double *dpos);

extern int	pos_to_indices(simstruct *sim, double *pos,
			int *index, PIXTYPE *weight);

