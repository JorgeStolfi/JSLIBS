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

  This program illustrates image segmentations using boundary
  tracking. One can further optimize it by stoping the IFT as soon as
  the target pixel is removed from queue. Another possible
  optimization is to precompute the gradient function instead of
  computing it on-the-fly.

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
  Pixel p1,pn;

  /* The following block must the remarked when using non-linux machines */

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /* Image preparation for contour tracking */

  aux1     = ReadImage("../data/caudate.pgm");  
  aux2     = CopyImage(aux1);
  A        = Circular(1.5);  
  img      = GaussStretch(aux1,140.0,30.0);
  lambda   = LambdaInit(img->ncols,img->nrows);
  handicap = CreateImage(img->ncols,img->nrows);

  /* First segment: (67,68) to (16,20) */

  p1.x = 67; p1.y = 68;
  pn.x = 16; pn.y = 20;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  DestroyImage(&handicap);
  handicap = CopyImage(aimg->cost);
  SetLambda(lambda,p1.x,p1.y,NIL);
  DestroyAnnImg(&aimg);

  /* Second segment: (16,20) to (67,68) */

  p1.x = 16; p1.y = 20;
  pn.x = 67; pn.y = 68;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  WriteImage(aux2,"ctrack_caudate.pgm");

  DestroyImage(&aux1);
  DestroyImage(&aux2);
  DestroyImage(&lambda);
  DestroyImage(&handicap);
  DestroyAnnImg(&aimg);
  DestroyImage(&img);
  DestroyAdjRel(&A);

  /* Image preparation for contour tracking */

  aux1     = ReadImage("../data/cerebellum.pgm");  
  aux2     = CopyImage(aux1);
  A        = Circular(1.5);  
  img      = GaussStretch(aux1,60.0,30.0);
  lambda   = LambdaInit(img->ncols,img->nrows);
  handicap = CreateImage(img->ncols,img->nrows);

  /* First segment: (271,51) to (22,83) */

  p1.x = 271; p1.y = 51;
  pn.x = 22;  pn.y = 83;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  DestroyImage(&handicap);
  handicap = CopyImage(aimg->cost);
  SetLambda(lambda,p1.x,p1.y,NIL);
  DestroyAnnImg(&aimg);

  /* Second segment: (22,83) to (64,140) */

  p1.x = 22; p1.y = 83;
  pn.x = 64; pn.y = 140;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  DestroyImage(&handicap);
  handicap = CopyImage(aimg->cost);
  SetLambda(lambda,p1.x,p1.y,NIL);
  DestroyAnnImg(&aimg);

  /* Third segment: (64,140) to (163,109) */

  p1.x = 64;  p1.y = 140;
  pn.x = 163; pn.y = 109;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  DestroyImage(&handicap);
  handicap = CopyImage(aimg->cost);
  SetLambda(lambda,p1.x,p1.y,NIL);
  DestroyAnnImg(&aimg);

  /* Fourth segment: (163,109) to (286,120) */

  p1.x = 163; p1.y = 109;
  pn.x = 286; pn.y = 120;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255);
  DestroyImage(&handicap);
  handicap = CopyImage(aimg->cost);
  SetLambda(lambda,p1.x,p1.y,NIL);
  DestroyAnnImg(&aimg);

  /* Last segment: (286,120) to (271,51) */

  p1.x = 286;  p1.y = 120;
  pn.x = 271;  pn.y = 51;

  SetLambda(lambda,p1.x,p1.y,1);
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,handicap,A,Fctrack);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Contourtrack in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  DrawPoint(aux2,pn.x,pn.y,2.7,255);
  DrawPath(aux2,aimg->pred,pn.x+img->tbrow[pn.y],255); 
  WriteImage(aux2,"ctrack_cerebellum.pgm");

  DestroyImage(&aux1);
  DestroyImage(&aux2);
  DestroyImage(&lambda);
  DestroyImage(&handicap);
  DestroyAnnImg(&aimg);
  DestroyImage(&img);
  DestroyAdjRel(&A);

  /* The following block must the remarked when using non-linux machines */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  fclose(ftime);

  return(0);
}
