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

  This program illustrates multiscale skeletons based on the quasi-DVD
  (discrete Voronoi diagram), which consists of 8-connected influence
  zones for a set of contour pixels. One can optimize it by
  accumulating the absolute displacements dx and dy along the path,
  and by using a look-up table for all possible squared values inside
  the image. In other words, one can use path-cost function:

  (sum |x_{i+1}-x_i|)^2 + (sum |y_{i+1}-y_i|)^2, 

  for a path <(x_1,y_1),(x_2,y_2),...,(x_n,y_n)>, where (x_i,y_i) is
  the coordinate of pixel i.
*/

#include "gift.h"

Image *Perimeter(Image *bin)
{
  int p,n;
  Image *cont,*perim;
  int *hist;

  cont  = LambdaContour(bin);
  n     = cont->ncols*cont->nrows;
  for (p=0; p < n; p++) { 
    if(cont->val[p]==NIL) 
      cont->val[p]=0;
  }

  hist  = AllocIntArray(MaximumValue(cont)+1);
  perim = CreateImage(cont->ncols,cont->nrows);
  for (p=0; p < n; p++)
    hist[cont->val[p]]++;
  for (p=0; p < n; p++)
    if (cont->val[p] > 0)
      perim->val[p] = hist[cont->val[p]];

  free(hist);
  DestroyImage(&cont);

  return(perim);
}

Image *MSSkel(AnnImg *aimg)
{
  int i,p,q,n,maxd1,maxd2,d1,d2,MaxD;
  Pixel u,v;
  Image *msskel,*cont=NULL,*perim=NULL,*label=NULL;
  AdjRel *A;

  /* Compute MS Skeletons */

  cont   = LambdaContour(aimg->img);
  label  = LambdaContPixel(aimg->img,1);
  n      = aimg->label->ncols*aimg->label->nrows;
  for (p=0; p < n; p++) { 
    if(cont->val[p]==NIL) 
      cont->val[p]=0;
    if(label->val[p]==NIL) 
      label->val[p]=0;
  }
  perim  = Perimeter(aimg->img);
  A      = Circular(1.0);
  msskel = CreateImage(aimg->label->ncols,aimg->label->nrows);

  MaxD = INT_MIN;
  for (p=0; p < n; p++) {
    if (aimg->pred->val[p] != p) { 
      u.x = p%aimg->label->ncols;
      u.y = p/aimg->label->ncols;
      maxd1 = maxd2 = INT_MIN;
      for (i=1; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(aimg->label,v.x,v.y)){
	  q = v.x + aimg->label->tbrow[v.y];
	  if (cont->val[aimg->label->val[p]] == cont->val[aimg->label->val[q]]){ 
	    d2   = label->val[aimg->label->val[q]]-
	      label->val[aimg->label->val[p]];
	    if (d2 > (perim->val[aimg->label->val[p]]-1-d2)){
	      d2 = (perim->val[aimg->label->val[p]]-1-d2);
	    } 
	    if (d2 > maxd2){
	      maxd2 = d2;
	    }
	  } else {
	    d1 = cont->val[aimg->label->val[q]] - cont->val[aimg->label->val[p]];
	    if (d1 > maxd1) 
	      maxd1 = d1;
	  }
	}
      }
      if (maxd1 > 0) {
	msskel->val[p] = INT_MAX;
      } else {
	msskel->val[p] = MAX(maxd2,0);
	if (msskel->val[p] > MaxD)
	  MaxD = msskel->val[p];    
      }
    }
  }

  for (p=0; p < n; p++) { /* Set up SKIZ */
    if (msskel->val[p] == INT_MAX)
      msskel->val[p] = MaxD + 1;
  }

  DestroyImage(&cont);
  DestroyImage(&perim);
  DestroyImage(&label);
  DestroyAdjRel(&A);

  return(msskel);
}

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *img,*lambda,*aux,*mssk,*label;
  int maxval, p, n;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */
  /*----------------------------------------------------------------------*/
  
  aux = ReadImage("../data/msskel.pgm");
  img = Threshold(aux,200,255); /* binary must be 0/1 */
  WriteImage(img,"msskel_a.pgm");
  A   = Circular(1.5);
  DestroyImage(&aux);
  lambda = LambdaContPixel(img,0);
  
  gettimeofday(&tic,NULL);  
  aimg=IFTFIFO(img,lambda,NULL,A,Fedt);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"EDT in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  label = LambdaContPixel(aimg->img,1);
  n     = label->ncols*label->nrows;
  for (p=0; p < n; p++){
    if (aimg->img->val[p]==1)
      label->val[p] = label->val[aimg->label->val[p]];
    else
      label->val[p]=0;
  }
  WriteImage(label,"msskel_b.pgm");
  DestroyImage(&label);

  aux  = MSSkel(aimg);
  mssk = Mult(img,aux);
  DestroyImage(&aux);

  WriteImage(mssk,"msskel_c.pgm");


  maxval = MaximumValue(mssk);
  aux = Threshold(mssk,(int)(0.05*maxval),maxval);
  WriteImage(aux,"msskel_d.pgm");
  DestroyImage(&aux);

  aux = Threshold(mssk,(int)(0.35*maxval),maxval);
  WriteImage(aux,"msskel_e.pgm");
  DestroyImage(&aux);

  aux = Threshold(mssk,(int)(0.75*maxval),maxval);
  WriteImage(aux,"msskel_f.pgm");
  DestroyImage(&aux);

  DestroyImage(&aux);
  DestroyImage(&img);
  DestroyImage(&mssk);
  DestroyAnnImg(&aimg);
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
