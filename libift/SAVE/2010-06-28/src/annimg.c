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
  -------------------------------------------------------------------------
  written by A.X. Falcão <afalcao@ic.unicamp.br>

  This program is a collection of functions to create, destroy, and
  manipulate an annotated image.
*/

#include "annimg.h"

AnnImg *CreateAnnImg(Image *img)
{
  AnnImg *aimg=NULL;

  aimg = (AnnImg *) calloc(1,sizeof(AnnImg));
  if (aimg == NULL)
    Error(MSG1,"CreateAnnImg");

  aimg->img   = img;
  aimg->cost  = CreateImage(img->ncols,img->nrows);    
  aimg->label = CreateImage(img->ncols,img->nrows);
  aimg->pred  = CreateImage(img->ncols,img->nrows);

  return(aimg);
}

void  DestroyAnnImg(AnnImg **aimg)
{
  AnnImg *aux;
  
  aux = *aimg;
  if (aux != NULL){
    DestroyImage(&(aux->cost));
    DestroyImage(&(aux->label));
    DestroyImage(&(aux->pred));
    free(aux);
    *aimg = NULL;
  }
}

Image *RootMap(Image *pred)
{
  Image *root=NULL;
  int p,n;

  root = CopyImage(pred);

  n = root->ncols*root->nrows;
  for (p=0; p < n; p++) 
    root->val[p] = Root(pred,p);
  return(root);
}

int Root(Image *pred, int p)
{
  if (pred->val[p]==NIL)
    return(p);
  else 
    return(Root(pred,pred->val[p]));    
}
