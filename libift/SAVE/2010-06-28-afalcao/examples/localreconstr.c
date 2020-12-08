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

  This program illustrates image segmentations using local morphological
  reconstruction.

*/

#include "gift.h"

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *img,*lambda,*handicap;
  Image *aux;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */

  /* Local superior reconstruction */

  img      = ReadImage("../data/knee.pgm");  
  aux      = CopyImage(img);
  DrawPoint(aux,65,71,3.0,255);
  WriteImage(aux,"knee_a.pgm");
  DestroyImage(&aux);
  A        = Circular(1.5);  
  lambda   = LambdaInit(img->ncols,img->nrows);
  handicap = CreateImage(img->ncols,img->nrows);
  handicap->val[65+handicap->tbrow[71]]=141;
  SetLambda(lambda,65,71,1);

  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(img,lambda,handicap,A,Flsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"LocalSupRec in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  aux = DrawBorder(aimg->cost,aimg->label,0);
  WriteImage(aux,"knee_b.pgm");

  DestroyImage(&aux);
  DestroyImage(&img);
  DestroyAdjRel(&A);
  DestroyImage(&lambda);
  DestroyImage(&handicap);
  DestroyAnnImg(&aimg);

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
