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

  This program is a collection of labeling functions. Seed pixels are
  considered those whose label is different from 0.
*/

#include "seeds.h"

Image *LambdaInit(int ncols, int nrows)
{
  Image *lambda=CreateImage(ncols,nrows);
  int p,n=ncols*nrows;

  for (p=0; p < n; p++) 
    lambda->val[p]=NIL;

  return(lambda);
}

Image *LambdaTrivial(int ncols, int nrows)
{
  Image *lambda=CreateImage(ncols,nrows);
  int p,n;

  n = ncols*nrows;
  for (p=0; p < n; p++) 
    lambda->val[p] = p;

  return(lambda);
}

Image  *LambdaRegMin(Image *img, AdjRel *A)
{
  AnnImg *aimg;
  AdjRel *A4;
  Image *lambda=LambdaTrivial(img->ncols,img->nrows);
  int p,n;
    
  aimg=IFTLIFO(img,lambda,NULL,A,Fini);
  n = img->ncols*img->nrows;
  for(p=0; p < n; p++) 
    if (aimg->pred->val[p]==NIL)
      aimg->label->val[p]=1;
    else
      aimg->label->val[p]=0;
  DestroyImage(&lambda);
  A4     = Circular(1.0);
  lambda = LambdaComp(aimg->label,A4);
  DestroyAdjRel(&A4);
  DestroyAnnImg(&aimg);

  return(lambda);
}

Image *LambdaFrame(int ncols, int nrows, int value) 
{
  int x,y; 
  Image *lambda=LambdaInit(ncols,nrows);

  for(x=0; x < ncols; x++) 
    lambda->val[x]=lambda->val[x+lambda->tbrow[nrows-1]]=value;
  for(y=0; y < nrows; y++) 
    lambda->val[lambda->tbrow[y]]=lambda->val[ncols-1+lambda->tbrow[y]]=value;
  return(lambda);
}

Image *LambdaContour(Image *bin)
{
  Image *bndr=NULL;
  Image *color=NULL,*pred=NULL,*lambda=NULL;
  int p=0,q,r,i,j,left=0,right=0,n,*LIFO,last,l=1;
  AdjRel *A,*L,*R;
  Pixel u,v,w;
  
  A     = Circular(1.0);
  n     = bin->ncols*bin->nrows;
  bndr  = CreateImage(bin->ncols,bin->nrows);
  for (p=0; p < n; p++){
    if (bin->val[p]==1){
      u.x = p%bin->ncols;
      u.y = p/bin->ncols;
      for (i=1; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(bin,v.x,v.y)){
	  q = v.x + bin->tbrow[v.y];
	  if (bin->val[q]==0){
	    bndr->val[p]=1;
	    break;
	  }
	} else {
	    bndr->val[p]=1;
	    break;
	}
      }
    }
  }
  DestroyAdjRel(&A);

  A      = Circular(1.5);
  L      = LeftSide(A);
  R      = RightSide(A);
  lambda = LambdaInit(bndr->ncols,bndr->nrows);
  color  = CreateImage(bndr->ncols,bndr->nrows);
  pred   = CreateImage(bndr->ncols,bndr->nrows);
  LIFO   = AllocIntArray(n);
  last   = NIL;
  for (j=0; j < n; j++){
    if ((bndr->val[j]==1)&&
	(color->val[j]!=BLACK)&&
	ValidContPoint(bin,L,R,j)){      
      last++; LIFO[last]    = j;
      color->val[j] = GRAY;
      pred->val[j] = j;
      while(last != NIL){
	p = LIFO[last];	last--;	
	color->val[p]=BLACK;
	u.x = p%bndr->ncols;
	u.y = p/bndr->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if (ValidPixel(bndr,v.x,v.y)){
	    q = v.x + bndr->tbrow[v.y];
	    if ((q==j)&&(pred->val[p]!=j)){
	      last = NIL;
	      break;
	    }
	    w.x = u.x + L->dx[i]; 
	    w.y = u.y + L->dy[i];
	    if (ValidPixel(bndr,w.x,w.y))
	      left = w.x + bndr->tbrow[w.y];
	    else
	      left = -1;
	    w.x = u.x + R->dx[i]; 
	    w.y = u.y + R->dy[i];
	    if (ValidPixel(bndr,w.x,w.y))
	      right = w.x + bndr->tbrow[w.y];
	    else
	      right = -1;
	    
	    if ((bndr->val[q]==1)&&
		(color->val[q] != BLACK)&&
		(((left!=-1)&&(right!=-1)&&(bin->val[left]!=bin->val[right]))||
		 ((left==-1)&&(right!=-1)&&(bin->val[right]==1)) ||
		 ((right==-1)&&(left!=-1)&&(bin->val[left]==1))) ) {
	      pred->val[q] = p;
	      if (color->val[q] == WHITE){
		last++; LIFO[last] = q;
		color->val[q]=GRAY;
	      }
	    } 
	  }
	}	
      }
      r = p;
      while(pred->val[p]!=p){
	lambda->val[p] = l;
	p = pred->val[p];
      } 
      if (r != p){
	lambda->val[p] = l;
	l++;
      }
    }
  }

  DestroyAdjRel(&A);
  DestroyAdjRel(&L);
  DestroyAdjRel(&R);
  DestroyImage(&bndr);
  DestroyImage(&color);
  DestroyImage(&pred);
  free(LIFO);
  return(lambda);
}

Image *LambdaContPixel(Image *bin, char label)
{
  Image *bndr=NULL;
  Image *color=NULL,*pred=NULL,*lambda=NULL;
  int p=0,q,r,i,j,n,left=0,right=0,*LIFO,last,l;
  AdjRel *A,*L,*R;
  Pixel u,v,w;
  
  A     = Circular(1.0);
  n     = bin->ncols*bin->nrows;
  bndr  = CreateImage(bin->ncols,bin->nrows);
  for (p=0; p < n; p++){
    if (bin->val[p]==1){
      u.x = p%bin->ncols;
      u.y = p/bin->ncols;
      for (i=1; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(bin,v.x,v.y)){
	  q = v.x + bin->tbrow[v.y];
	  if (bin->val[q]==0){
	    bndr->val[p]=1;
	    break;
	  }
	} else {
	    bndr->val[p]=1;
	    break;
	}
      }
    }
  }  
  DestroyAdjRel(&A);

  A      = Circular(1.5);
  L      = LeftSide(A);
  R      = RightSide(A);
  lambda = LambdaInit(bndr->ncols,bndr->nrows);
  color  = CreateImage(bndr->ncols,bndr->nrows);
  pred   = CreateImage(bndr->ncols,bndr->nrows);
  n      = bndr->ncols*bndr->nrows;
  LIFO   = AllocIntArray(n);
  last   = NIL;
  for (j=0; j < n; j++){
    if ((bndr->val[j]==1)
	&&(color->val[j]!=BLACK)
	&&ValidContPoint(bin,L,R,j)){
      last++;
      LIFO[last]    = j;
      color->val[j] = GRAY;
      pred->val[j] = j;
      while(last != NIL){
	p = LIFO[last]; last--;	
	color->val[p]=BLACK;
	u.x = p%bndr->ncols;
	u.y = p/bndr->ncols;
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if (ValidPixel(bndr,v.x,v.y)){
	    q = v.x + bndr->tbrow[v.y];
	    if ((q==j)&&(pred->val[p]!=j)){
	      last = NIL;
	      break;
	    }
	    
	    w.x = u.x + L->dx[i]; 
	    w.y = u.y + L->dy[i];
	    if (ValidPixel(bndr,w.x,w.y))
	      left = w.x + bndr->tbrow[w.y];
	    else
	      left = -1;
	    w.x = u.x + R->dx[i]; 
	    w.y = u.y + R->dy[i];
	    if (ValidPixel(bndr,w.x,w.y))
	      right = w.x + bndr->tbrow[w.y];
	    else
	      right = -1;
	    
	    if ((bndr->val[q]==1)&&
		(color->val[q] != BLACK)&&
		(((left!=-1)&&(right!=-1)&&(bin->val[left]!=bin->val[right]))||
		 ((left==-1)&&(right!=-1)&&(bin->val[right]==1)) ||
		 ((right==-1)&&(left!=-1)&&(bin->val[left]==1)))){ 
	      pred->val[q] = p;
	      if (color->val[q] == WHITE){
		last++;
		LIFO[last] = q;
		color->val[q]=GRAY;
	      }
	    } 
	  }
	}	
      }
      r = p;
      if (label==1){
	l = 1;
	while(pred->val[p]!=p){
	  lambda->val[p] = l;
	  p = pred->val[p];
	  l++;
	}
	if (r != p) {
	  lambda->val[p] = l;
	}
      }else{ /* use trivial labeling function */
	while(pred->val[p]!=p){
	  lambda->val[p] = p;
	  p = pred->val[p];
	}
	if (r != p) {
	  lambda->val[p] = p;
	}
      }
    }
  }

  DestroyAdjRel(&A);
  DestroyAdjRel(&L);
  DestroyAdjRel(&R);
  DestroyImage(&bndr);
  DestroyImage(&color);
  DestroyImage(&pred);
  free(LIFO);
  return(lambda);
} 

Image *LambdaComp(Image *bin, AdjRel *A)
{
  Image *lambda=NULL;
  int i,j,n,p,q,l=1;
  int *FIFO=NULL;
  int first=0,last=0;
  Pixel u,v;
  
  lambda = LambdaInit(bin->ncols,bin->nrows);
  n     = bin->ncols*bin->nrows;
  FIFO  = AllocIntArray(n);
  for (j=0; j < n; j++){
    if ((bin->val[j]==1)&&(lambda->val[j]==NIL)){
      lambda->val[j]=l;
      FIFO[last]=j;      
      last++;      
      while(first != last){
	p = FIFO[first];
	first++;
	u.x = p%bin->ncols;
	u.y = p/bin->ncols;
	for (i=0; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if (ValidPixel(bin,v.x,v.y)){
	    q = v.x + bin->tbrow[v.y];
	    if ((bin->val[q]==1)&&(lambda->val[q] == NIL)){
	      lambda->val[q] = lambda->val[p];
	      FIFO[last] = q;
	      last++;
	    }
	  }
	}
      }
      l++;
      first=last=0;
    }
  }
  
  free(FIFO);

  return(lambda);
}

void SetLambda(Image *lambda, int x, int y, int value)
{
  if (ValidPixel(lambda,x,y))
    lambda->val[x+lambda->tbrow[y]]=value;
}


bool ValidContPoint(Image *bin, AdjRel *L, AdjRel *R, int p)
{
  int i,q,n,left,right;
  Pixel u,v,l,r;
  bool found=false;

  u.x = p%bin->ncols;
  u.y = p/bin->ncols;
  n   = L->n;

  for (i=0; i < n; i++) {
    v.x = u.x + L->dx[i];
    v.y = u.y + L->dy[i];
    if (ValidPixel(bin,v.x,v.y)){
      q = v.x + bin->tbrow[v.y];
      if ((bin->val[q]==1)&&(p!=q)){
	l.x = u.x + L->dx[i]; 
	l.y = u.y + L->dy[i];
	r.x = u.x + R->dx[i]; 
	r.y = u.y + R->dy[i];	
	if (ValidPixel(bin,l.x,l.y))
	  left = l.x + bin->tbrow[l.y];
	else
	  left = -1;
	if (ValidPixel(bin,r.x,r.y))
	  right = r.x + bin->tbrow[r.y];
	else
	  right = -1;
	if (((left!=-1)&&(right!=-1)&&(bin->val[left]!=bin->val[right]))||
	    ((left==-1)&&(right!=-1)&&(bin->val[right]==1)) ||
	    ((right==-1)&&(left!=-1)&&(bin->val[left]==1))){
	  found = true;
	  break;
	}
      }
    }
  }
  
  return(found);
}
