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

  This program illustrates image segmentations using the watershed transform.
*/

#include "gift.h"

Image *Gradient(Image *img)
{
  Image *grad;
  Pixel s,u,v,g;
  int s1,u1,v1;

  grad = CreateImage(img->ncols,img->nrows);
  for (u.y=0; u.y < grad->nrows; u.y++) 
    for (u.x=0; u.x < grad->ncols; u.x++) {
      v.x = ((u.x + 1)==grad->ncols)?u.x:u.x+1; 
      v.y = u.y;
      s.x = ((u.x-1)==-1)?u.x:u.x-1;
      s.y = u.y;
      u1  = u.x + grad->tbrow[u.y];
      v1  = v.x + grad->tbrow[v.y];
      s1  = s.x + grad->tbrow[s.y];
      g.x = ((img->val[u1]+img->val[v1])/2 - 
	     (img->val[u1]+img->val[s1])/2);  
      v.x = u.x;
      v.y = ((u.y + 1)==grad->nrows)?u.y:u.y+1; 
      s.x = u.x;
      s.y = ((u.y-1)==-1)?u.y:u.y-1;
      u1  = u.x + grad->tbrow[u.y];
      v1  = v.x + grad->tbrow[v.y];
      s1  = s.x + grad->tbrow[s.y];
      g.y = ((img->val[u1]+img->val[v1])/2 - 
	     (img->val[u1]+img->val[s1])/2);  
      grad->val[u1] = (int)sqrt(g.x*g.x + g.y*g.y);
    }
  return(grad);
}

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *img,*lambda,*handicap;
  Image *aux1,*aux2,*edt,*grad;
  int p,n,Imax;
  Pixel u;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */
  /* Example 1: 
     Segmenting the left caudate nucleus in an MR image of the brain. */

  img  = ReadImage("../data/caudate.pgm");
  WriteImage(img,"caudate_a.pgm");

  /* Preprocessing by Gaussian stretching followed by gradient. */

  aux1 = GaussStretch(img,140.0,30.0);
  A    = Circular(1.5);
  grad = Gradient(aux1);

  /* Set labeling functions lambda and handicaps */

  lambda   = LambdaInit(img->ncols,img->nrows);
  handicap = CreateImage(img->ncols,img->nrows);
  lambda->val[45+img->tbrow[37]] =1;
  lambda->val[82+img->tbrow[44]] =2;

  aux2 = CopyImage(grad);
  DrawPoint(aux2,45,37,2.7,255);
  DrawPoint(aux2,82,44,2.7,255);
  WriteImage(aux2,"caudate_b.pgm");
  DestroyImage(&aux2);
  
  /* compute the IFT-watershed with FIFO tiebreaking */
  
  gettimeofday(&tic,NULL);  
  aimg     = IFTFIFO(grad,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Watershed in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  aux2     = DrawBorder(img,aimg->label,255);
  WriteImage(aux2,"caudate_c.pgm");
  DestroyImage(&img);
  DestroyImage(&aux1);
  DestroyImage(&grad);
  DestroyImage(&aux2);
  DestroyAdjRel(&A);
  DestroyImage(&lambda);
  DestroyImage(&handicap);
  DestroyAnnImg(&aimg);

  /* Example 2: 
     Separating toching cells by a geometric criterion */

  img      = ReadImage("../data/cells.pgm");  
  WriteImage(img,"cells_a.pgm");

  /* Preprocessing by thresholding and morphological closing */

  A        = Circular(1.5);  
  aux1     = Threshold(img,0,221);
  aux2     = Close(aux1,A);

  /* Set labeling functions lambda and handicaps are not needed here. */

  lambda   = LambdaContPixel(aux2,0);
  DestroyImage(&aux1);

  /* Compute IFT-Euclidean distance transform */

  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(aux2,lambda,NULL,A,Fedt);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"EDT in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  aux1 = SQRT(aimg->cost);
  edt  = Complement(aux1);

  DestroyAnnImg(&aimg);
  DestroyImage(&lambda);
  DestroyImage(&aux2);
  DestroyImage(&aux1);

  /* Find the regional minima of the EDT complement to use as seeds
     for the watershed transform. */

  lambda = LambdaTrivial(img->ncols,img->nrows);
  gettimeofday(&tic,NULL);  
  aimg   = IFTLIFO(edt,lambda,NULL,A,Fini);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"RegMin in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);
   
  n = img->ncols*img->nrows;
  DestroyImage(&lambda);
  lambda=LambdaInit(img->ncols,img->nrows);
  aux1 = CopyImage(edt);
  Imax = MaximumValue(aux1);
  for (p=0; p < n; p++) 
    if (aimg->pred->val[p]==NIL){ /* p belongs to a regional minimum */
      u.x = p%img->ncols;
      u.y = p/img->ncols;
      DrawPoint(aux1,u.x,u.y,3.0,Imax);
      if ((u.x==0)||(u.y==0)||
	  (u.x==img->ncols-1)||
	  (u.y==img->nrows-1)){ /* the minimum belongs to the image
                                  frame */
	SetLambda(lambda,u.x,u.y,1);
      }else /* the minimum belongs to a cell */
	SetLambda(lambda,u.x,u.y,p);
    }
  WriteImage(aux1,"cells_b.pgm");
  DestroyImage(&aux1);
  DestroyAnnImg(&aimg);

  /* watershed transform (superior reconstruction) */

  handicap = CopyImage(edt);
  gettimeofday(&tic,NULL);  
  aimg   = IFTFIFO(edt,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"SupRec(Watershed) in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  aux1 = DrawBorder(img,aimg->label,0);
  WriteImage(aux1,"cells_c.pgm");
  DestroyImage(&aux1);  
  DestroyImage(&handicap);
  DestroyAnnImg(&aimg);
  DestroyImage(&edt);
  DestroyAnnImg(&aimg);
  DestroyImage(&img);
  DestroyAdjRel(&A);
  DestroyImage(&lambda);

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
