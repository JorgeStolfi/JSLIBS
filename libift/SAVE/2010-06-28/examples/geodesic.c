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

  This program computes the geodesic path from a source set to a
  target set of pixels. One can optimize it by stopping the IFT when
  the first pixel of the target set is removed from queue.

*/

#include "gift.h"

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *img,*lambda,*aux;
  int p,n,dst=0,mincost;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */

  
  img      = ReadImage("../data/geomask.pgm");  
  aux      = Threshold(img,1,255);
  A        = Circular(1.5);  
  lambda   = LambdaInit(img->ncols,img->nrows);
  n        = img->ncols*img->nrows;
  for(p=0; p < n; p++) {
    if (img->val[p]==184) /* take as seed */
      lambda->val[p] = 1;
  }

  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(aux,lambda,NULL,A,Fgeo);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Geodesic in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);
  
  mincost = INT_MAX;
  for(p=0; p < n; p++) 
    if (img->val[p]==68) /* destination */
      if (aimg->cost->val[p] < mincost) {
	mincost = aimg->cost->val[p];
	dst     = p;
      }
  DrawPath(img,aimg->pred,dst,0); 
  WriteImage(img,"geodesic.pgm");

  DestroyImage(&lambda);
  DestroyAnnImg(&aimg);
  DestroyImage(&img);
  DestroyImage(&aux);
  DestroyAdjRel(&A);

  /* The following block must the remarked when using non-linux machines 

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  */

  fclose(ftime);

  return(0);
}
