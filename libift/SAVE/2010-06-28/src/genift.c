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

  This program contains path-cost functions and two algorithms for
  image foresting transform with different tiebreaking policies: FIFO
  and LIFO.
*/

#include "genift.h"

int Fsuprec(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q)
{
  if (p==q) 
    if (lambda->val[p] >= 0) /* seed pixel */      
      return(handicap->val[p]);
    else
      return(INT_MAX);
  else
    return(MAX(aimg->cost->val[p],aimg->img->val[q]));
}

int Flsuprec(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q)
{
  if (p==q) 
    if (lambda->val[p] >= 0) /* seed pixel */      
      return(handicap->val[p]);
    else
      return(INT_MAX);
  else{
    if (aimg->cost->val[p] > aimg->img->val[q])
      return(aimg->cost->val[p]);
    else
      return(INT_MAX);
  }
}

int Fini(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q) 
{
  if (p==q)
      return(aimg->img->val[p]);
  else{
    if (aimg->img->val[p] <= aimg->img->val[q])
      return(aimg->cost->val[p]);
    else
      return(INT_MAX);
  }
}

int Fedt(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q)
{
  int dx,dy,r;
  Pixel u,v;
 
  if (p==q) {
    if (lambda->val[p] >= 0) /* seed pixel */      
      return(0);
    else
      return(INT_MAX);
  }else{
    r = aimg->label->val[p];
    u.x = (r%aimg->img->ncols);
    u.y = (r/aimg->img->ncols);
    v.x = q%aimg->img->ncols;
    v.y = q/aimg->img->ncols;
    dx = v.x-u.x;
    dy = v.y-u.y;
    return(dx*dx + dy*dy); 
  }
}

int Fctrack(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q)
{  
  int d,left,right,dx,dy;
  float dxl,dyl,dxr,dyr;
  Pixel u,v,l,r;
  if (p==q) {
    if (lambda->val[p] >= 0) /* seed pixel */      
      return(handicap->val[p]);
    else
      return(INT_MAX);
  }else{
    u.x = p%aimg->img->ncols;
    u.y = p/aimg->img->ncols;
    v.x = q%aimg->img->ncols;
    v.y = q/aimg->img->ncols;
    dx  = v.x - u.x; dy = v.y - u.y;
    d   = dx*dx + dy*dy;
    dxl = ((float)dy/d);
    dyl = (-(float)dx/d);
    dxr = (-(float)dy/d);
    dyr = ((float)dx/d);
    l.x = (int)((float)(u.x+v.x)/2.0 + dxl); 
    l.y = (int)((float)(u.y+v.y)/2.0 + dyl);
    if (!ValidPixel(aimg->img,l.x,l.y)) {l.x=u.x; l.y=u.y;}
    r.x = (int)((float)(u.x+v.x)/2.0 + dxr); 
    r.y = (int)((float)(u.y+v.y)/2.0 + dyr);
    if (!ValidPixel(aimg->img,r.x,r.y)) {r.x=u.x; r.y=u.y;}
    left  = l.x + aimg->img->tbrow[l.y];
    right = r.x + aimg->img->tbrow[r.y];
    return(aimg->cost->val[p] + 
	   255 - (aimg->img->val[left] - aimg->img->val[right]));
  }
}

int Fgeo(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q)
{
  Pixel u,v;
  int dx,dy;

  if (p==q) {
    if (lambda->val[p] >= 0) /* seed pixel */      
      return(0);
    else
      return(INT_MAX);
  }else {
    if (aimg->img->val[p]!=aimg->img->val[q])
      return(INT_MAX);
    u.x = p%aimg->img->ncols;
    u.y = p/aimg->img->ncols;
    v.x = q%aimg->img->ncols;
    v.y = q/aimg->img->ncols;
    dx  = abs(u.x-v.x);
    dy  = abs(u.y-v.y);
    return(aimg->cost->val[p]+5*MAX(dx,dy) + 2*MIN(dx,dy));
  }
}

AnnImg *IFTFIFO(Image *img, Image *lambda, Image *handicap, AdjRel *A, PathCost Pcost)
{
  Queue *Q=NULL;
  int i,p,q,n,cost;
  Pixel u,v;
  AnnImg *aimg=NULL;

  aimg = CreateAnnImg(img);
  n = aimg->img->ncols*aimg->img->nrows;
  Q = CreateQueue(QSIZE+1,n,aimg->cost->val);

  for (p=0; p < n; p++) {
    aimg->pred->val[p]  = NIL;
    aimg->label->val[p] = lambda->val[p]; 
    aimg->cost->val[p]  = Pcost(aimg,lambda,handicap,p,p);
    if (aimg->cost->val[p] < INT_MAX){
      InsertQueue(&Q,p);
    }
  }

  while(!EmptyQueue(Q)) {
    p=RemoveQueue(Q);
    u.x = p%aimg->img->ncols;
    u.y = p/aimg->img->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(aimg->img,v.x,v.y)){
	q = v.x + aimg->img->tbrow[v.y];
	if (aimg->cost->val[p] < aimg->cost->val[q]){
	  cost = Pcost(aimg,lambda,handicap,p,q); 
	  if (cost < aimg->cost->val[q]){
	    if (aimg->cost->val[q] == INT_MAX){
	      aimg->cost->val[q]  = cost;
	      InsertQueue(&Q,q);
	    }else
	      UpdateQueue(Q,q,cost);
	    aimg->pred->val[q]  = p;
	    aimg->label->val[q] = aimg->label->val[p];
	  }
	}
      }
    }
  }
  DestroyQueue(&Q);
  return(aimg);
}

AnnImg *IFTLIFO(Image *img, Image *lambda, Image *handicap, AdjRel *A, PathCost Pcost)
{
  Queue *Q=NULL;
  int i,p,q,n,cost;
  Pixel u,v;
  AnnImg *aimg=NULL;

  aimg = CreateAnnImg(img);
  n = aimg->img->ncols*aimg->img->nrows;
  Q = CreateQueue(QSIZE+1,n,aimg->cost->val);
  SetTieBreak(Q,LIFOBREAK);

  for (p=0; p < n; p++) {
    aimg->pred->val[p]  = NIL;
    aimg->label->val[p] = lambda->val[p]; 
    aimg->cost->val[p]  = Pcost(aimg,lambda,handicap,p,p);
    InsertQueue(&Q,p);    
  }

  while(!EmptyQueue(Q)) {
    p=RemoveQueue(Q);
    u.x = p%aimg->img->ncols;
    u.y = p/aimg->img->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if (ValidPixel(aimg->img,v.x,v.y)){
	q = v.x + aimg->img->tbrow[v.y];
	if (Q->L.elem[q].color == GRAY){ 
	  cost = Pcost(aimg,lambda,handicap,p,q);     
	  if (cost <= aimg->cost->val[q]){
	    UpdateQueue(Q,q,cost);
	    aimg->pred->val[q]  = p;
	    aimg->label->val[q] = aimg->label->val[p];
	  } 
	}
      }
    }
  }

  DestroyQueue(&Q);
  return(aimg);
}

