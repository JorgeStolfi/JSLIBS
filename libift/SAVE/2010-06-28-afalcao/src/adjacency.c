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

  This program is a collection of functions to create and destroy
  different types of adjacency relations.

*/

#include "adjacency.h"

AdjRel *CreateAdjRel(int n)
{
  AdjRel *A=NULL;

  A = (AdjRel *) calloc(1,sizeof(AdjRel));
  if (A != NULL){
    A->dx = AllocIntArray(n);
    A->dy = AllocIntArray(n);
    A->n  = n;
  } else {
    Error(MSG1,"CreateAdjRel");
  }

  return(A);
}

void DestroyAdjRel(AdjRel **A)
{
  AdjRel *aux;

  aux = *A;
  if (aux != NULL){
    if (aux->dx != NULL) free(aux->dx);
    if (aux->dy != NULL) free(aux->dy);
    free(aux);
    *A = NULL;
  }   
}

AdjRel *Circular(float r)
{
  AdjRel *A=NULL;
  int i,j,k,n,dx,dy,r0,r2,d,i0=0;
  float *da,*dr,aux;

  n=0;

  r0 = (int)r;
  r2  = (int)(r*r);
  for(dy=-r0;dy<=r0;dy++)
    for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)) <= r2)
	n++;
	
  A = CreateAdjRel(n);
  i=0;
  for(dy=-r0;dy<=r0;dy++)
    for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)) <= r2){
	A->dx[i]=dx;
	A->dy[i]=dy;
	if ((dx==0)&&(dy==0))
	  i0 = i;
	i++;
      }

  /* Set clockwise */
  
  da = AllocFloatArray(A->n);
  dr = AllocFloatArray(A->n);
  for (i=0; i < A->n; i++) {
    dx = A->dx[i];
    dy = A->dy[i];
    dr[i] = (float)sqrt((dx*dx) + (dy*dy));
    if (i != i0){ 
      da[i] = atan2(-dy,-dx)*180.0/PI;
      if (da[i] < 0.0)
	da[i] += 360.0;
    }
  }
  da[i0] = 0.0;
  dr[i0] = 0.0;
  
  /* place central pixel at first */
  
  aux    = da[i0];
  da[i0] = da[0];
  da[0]  = aux;
  aux    = dr[i0];
  dr[i0] = dr[0];
  dr[0]  = aux;
  d         = A->dx[i0];
  A->dx[i0] = A->dx[0];
  A->dx[0]  = d;
  d         = A->dy[i0];
  A->dy[i0] = A->dy[0];
  A->dy[0]  = d;

  /* sort by angle */
  
  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if (da[j] < da[k]){
	k = j;
      }
    aux   = da[i];
    da[i] = da[k];
    da[k] = aux;
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d   = A->dx[i];
    A->dx[i] = A->dx[k];
    A->dx[k] = d;
    d        = A->dy[i];
    A->dy[i] = A->dy[k];
    A->dy[k] = d;
  }

  /* sort by radius for each angle */
  
  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if ((dr[j] < dr[k])&&(da[j]==da[k])){
	k = j;
      }
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d        = A->dx[i];
    A->dx[i] = A->dx[k];
    A->dx[k] = d;
    d        = A->dy[i];
    A->dy[i] = A->dy[k];
    A->dy[k] = d;
  }

  free(dr);
  free(da);

  return(A);
}

AdjRel *RightSide(AdjRel *A)
{
  AdjRel *R=NULL;
  int i;
  float d;

  /* Let p -> q be an arc represented by the increments dx,dy. Its
     right side is given by the increments Dx = -dy/d + dx/2 and Dy =
     dx/d + dy/2, where d=sqrt(dx²+dy²). */

  R = CreateAdjRel(A->n);
  for (i=0; i < R->n; i++){
    d  = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
    if (d != 0){
      R->dx[i] = (int) (((float)A->dx[i]/2.0)-((float)A->dy[i]/d));
      R->dy[i] = (int)(((float)A->dx[i]/d)+((float)A->dy[i]/2.0));
    }
  }

  return(R);
}

AdjRel *LeftSide(AdjRel *A)
{
  AdjRel *L=NULL;
  int i;
  float d;

  /* Let p -> q be an arc represented by the increments dx,dy. Its
     left side is given by the increments Dx = dy/d + dx/2 and Dy =
     -dx/d + dy/2, where d=sqrt(dx²+dy²). */

  L = CreateAdjRel(A->n);
  for (i=0; i < L->n; i++){
    d  = sqrt(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
    if (d != 0){
      L->dx[i] = (int)(((float)A->dx[i]/2.0)+((float)A->dy[i]/d));
      L->dy[i] = (int)(((float)A->dy[i]/2.0)-((float)A->dx[i]/d));
    }
  }
  
  return(L);
}
