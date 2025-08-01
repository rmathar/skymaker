#pragma once
/*
*				imaout.h
*
* Include file for imaout.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		24/05/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "fits/fitscat.h"

#include "prefs.h"

#include "simul.h"

extern char	prefsname[MAXCHAR];

/*--------------------------- structure definitions -------------------------*/

typedef struct structheadkey
  {
  char          name[80];               /* name */
  char          comment[80];            /* a comment */
  void          *ptr;                   /* pointer to the data */
  h_type        htype;                  /* standard ``h_type'' (display) */
  t_type        ttype;                  /* standard ``t_type'' (storage) */
  char		prefsname[80];		/* true name in the preferences */
  }             headkeystruct;

/*---------------------------------- protos --------------------------------*/
catstruct       *imaout_inithead(simstruct *sim);

void		imaout_write(simstruct *sim);

