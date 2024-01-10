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

  This program illustrates image segmentations using morphological
  reconstructions.
*/

#include "gift.h"

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *img,*lambda,*handicap;
  Image *aux1,*aux2;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */
  
  img      = ReadImage("../data/wrist.pgm");  
  A        = Circular(1.5);  
  aux1      = CopyImage(img);
  DrawPoint(aux1,102,72,3.0,255);
  WriteImage(aux1,"wrist_a.pgm");
  DestroyImage(&aux1);

  /* Inferior reconstruction */

  aux1      = Complement(img);
  lambda    = LambdaInit(img->ncols,img->nrows);
  SetLambda(lambda,102,72,1);
  handicap = CreateImage(img->ncols,img->nrows);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(aux1,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"InfRec in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);
  
  DestroyImage(&aux1);
  aux1    = Complement(aimg->cost);
  DestroyImage(&lambda);
  DestroyAnnImg(&aimg);
  lambda = LambdaFrame(aux1->ncols,aux1->nrows,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(aux1,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"SupRec in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  WriteImage(aimg->cost,"wrist_b.pgm");
  DestroyImage(&aux1);
  aux1 = Threshold(aimg->cost,60,225);
  aux2 = DrawBorder(img,aux1,255);
  WriteImage(aux2,"wrist_c.pgm");

  DestroyImage(&handicap);
  DestroyImage(&aux1);
  DestroyImage(&aux2);
  DestroyImage(&lambda);
  DestroyAnnImg(&aimg);
  DestroyImage(&img);
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
