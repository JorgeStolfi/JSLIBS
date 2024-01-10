/*
  Copyright (C) <2003> <Alexandre Xavier Falcão>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  please see full copyright in COPYING file.
  -------------------------------------------------------------------------
  written by A.X. Falcão <afalcao@ic.unicamp.br>, March 27th 2003

  This program is a collection of functions to measure the performance
  of other functions.

*/

#include "comptime.h"

timer *Tic(){ /* It marks the initial time */
  timer *tic=NULL;
  tic = (timer *)malloc(sizeof(timer));
  gettimeofday(tic,NULL); 
  return(tic);
}

timer *Toc(){ /* It marks the final time */
  timer *toc=NULL;
  toc = (timer *)malloc(sizeof(timer));
  gettimeofday(toc,NULL);
  return(toc);
}

float CTime(timer *tic, timer *toc) /* It computes the time difference */
{ 
  float t=0.0;
  if ((tic!=NULL)&&(toc!=NULL)){
    t = (toc->tv_sec-tic->tv_sec)*1000.0 + 
      (toc->tv_usec-tic->tv_usec)*0.001;
    free(tic);free(toc);
    tic=NULL; toc=NULL;
  }
  return(t);
}



