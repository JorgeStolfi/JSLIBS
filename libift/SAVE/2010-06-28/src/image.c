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

  This program is a collection of functions to create, destroy, and
  manipulate grayscale images.
*/

#include "image.h"

/* Internal functions */

int nc_fgets(char *s, int m, FILE *f);

Image *CreateImage(int ncols, int nrows)
{
  Image *img=NULL;
  int i;

  img = (Image *) calloc(1,sizeof(Image));
  if (img == NULL){
    Error(MSG1,"CreateImage");
  }

  img->val   = AllocIntArray(nrows*ncols);
  img->tbrow = AllocIntArray(nrows);

  img->tbrow[0]=0;
  for (i=1; i < nrows; i++)
    img->tbrow[i]=img->tbrow[i-1]+ncols;
  img->ncols = ncols;
  img->nrows = nrows;
 
 return(img);
}

void DestroyImage(Image **img)
{
  Image *aux;

  aux = *img;
  if(aux != NULL){
    if (aux->val != NULL)   free(aux->val); 
    if (aux->tbrow != NULL) free(aux->tbrow);
    free(aux);    
    *img = NULL;
  }
}

Image *ReadImage(char *filename)
{
  FILE *fp=NULL;
  unsigned char *value=NULL;
  char type[10];
  int  i,ncols,nrows,n;
  Image *img=NULL;
  char z[256];

  fp = fopen(filename,"r");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  fscanf(fp,"%s\n",type);
  if((strcmp(type,"P5")==0)){
    nc_fgets(z,255,fp);
    sscanf(z,"%d %d\n",&ncols,&nrows);
    n = ncols*nrows;
    nc_fgets(z,255,fp);
    sscanf(z,"%d\n",&i);
    value = (unsigned char *)calloc(n,sizeof(unsigned char));
    if (value != NULL){
      fread(value,sizeof(unsigned char),n,fp);
    }else{
      fprintf(stderr,"Insufficient memory in ReadImage\n");
      exit(-1);
    }
    fclose(fp);
    img = CreateImage(ncols,nrows);
    for (i=0; i < n; i++)
      img->val[i]=(int)value[i];
    free(value);
  }else{
    if((strcmp(type,"P2")==0)){
      nc_fgets(z,255,fp);
      sscanf(z,"%d %d\n",&ncols,&nrows);
      n = ncols*nrows;
      nc_fgets(z,255,fp);
      sscanf(z,"%d\n",&i);
      img = CreateImage(ncols,nrows);
      for (i=0; i < n; i++)
	fscanf(fp,"%d",&img->val[i]);
      fclose(fp);
    }else{
      fprintf(stderr,"Input image must be P2 or P5\n");
      exit(-1);
    }
  }

  return(img);
}

void WriteImage(Image *img,char *filename)
{
  FILE *fp;
  int i, n, Imax;

  fp = fopen(filename,"w");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  n    = img->ncols*img->nrows;
  if ((Imax=MaximumValue(img))==INT_MAX){
    Warning("Image with infinity values","WriteImage");
    Imax = INT_MIN;
    for (i=0; i < n; i++) 
      if ((img->val[i] > Imax)&&(img->val[i]!=INT_MAX))
	Imax = img->val[i];
    fprintf(fp,"P2\n");
    fprintf(fp,"%d %d\n",img->ncols,img->nrows);
    fprintf(fp,"%d\n",Imax+1);
  } else {
    fprintf(fp,"P2\n");
    fprintf(fp,"%d %d\n",img->ncols,img->nrows);
    if (Imax==0) Imax = 1;
    fprintf(fp,"%d\n",Imax);
  }
 
  for (i=0; i < n; i++) {
    if (img->val[i]==INT_MAX)
      fprintf(fp,"%d ",Imax+1);
    else
      fprintf(fp,"%d ",img->val[i]);
    if (((i+1)%17) == 0)
      fprintf(fp,"\n");
  }

  fclose(fp);

}

int MinimumValue(Image *img)
{
  int i,min,n;

  n = img->ncols*img->nrows;
  min = img->val[0];
  for (i=1; i < n; i++)
    if (img->val[i] < min)
      min = img->val[i];

  return(min);
}

int MaximumValue(Image *img)
{
  int i,max,n;

  max = img->val[0];
  n = img->ncols*img->nrows;
  for (i=1; i < n; i++)
    if (img->val[i] > max)
      max = img->val[i];

  return(max);
}

void SetImage(Image *img, int value)
{ 
  int i,n;
  n = img->ncols*img->nrows;
  for (i=0; i < n; i++){
    img->val[i]=value;
  }
}

bool ValidPixel(Image *img, int x, int y)
{
  if ((x >= 0)&&(x < img->ncols)&&
      (y >= 0)&&(y < img->nrows))
    return(true);
  else
    return(false);
}

Image *CopyImage(Image *img)
{
  Image *imgc;

  imgc = CreateImage(img->ncols,img->nrows);
  memcpy(imgc->val,img->val,img->ncols*img->nrows*sizeof(int));
  
  return(imgc);
}

Image *Threshold(Image *img, int lower, int higher)
{
  Image *bin=NULL;
  int p,n;

  bin = CreateImage(img->ncols,img->nrows);
  n = img->ncols*img->nrows;
  for (p=0; p < n; p++)
    if ((img->val[p] >= lower)&&(img->val[p] <= higher))
      bin->val[p]=1;
  return(bin);
}

Image *Dilate(Image *img, AdjRel *A)
{
  Image *dil=NULL;
  int p,q,max,i;
  Pixel u,v;

  dil = CreateImage(img->ncols,img->nrows);

  for (u.y=0; u.y < img->nrows; u.y++)
    for (u.x=0; u.x < img->ncols; u.x++){
      p   = u.x + img->tbrow[u.y];
      max = img->val[p];
      for (i=0; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(img,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  if (img->val[q] > max)
	    max = img->val[q];
	}
      }
      dil->val[p]=max;
    }
  return(dil);
}

Image *Erode(Image *img, AdjRel *A)
{
  Image *ero=NULL;
  int p,q,min,i;
  Pixel u,v;

  ero = CreateImage(img->ncols,img->nrows);

  for (u.y=0; u.y < img->nrows; u.y++)
    for (u.x=0; u.x < img->ncols; u.x++){
      p   = u.x + img->tbrow[u.y];
      min = img->val[p];
      for (i=0; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(img,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  if (img->val[q] < min)
	    min = img->val[q];
	}
      }
      ero->val[p]=min;
    }
  return(ero);
}

Image *Open(Image *img, AdjRel *A)
{
  Image *open=NULL,*ero=NULL;

  ero  = Erode(img,A);
  open = Dilate(ero,A);
  DestroyImage(&ero);

  return(open);
}

Image *Close(Image *img, AdjRel *A)
{
  Image *close=NULL,*dil=NULL;

  dil   = Dilate(img,A);
  close = Erode(dil,A);
  DestroyImage(&dil);

  return(close);
}

Image *MorphGrad(Image *img, AdjRel *A)
{
  Image *dil=NULL,*ero=NULL,*grad=NULL;

  dil  = Dilate(img,A);
  ero  = Erode(img,A);
  grad = Diff(dil,ero);

  DestroyImage(&dil);
  DestroyImage(&ero);

  return(grad);
}

Image  *Diff(Image *img1,Image *img2)
{
  Image *diff;
  int p,n;

  diff = CopyImage(img1);
  n = img1->ncols*img1->nrows;
  for (p=0; p < n; p++)
    diff->val[p] -= img2->val[p];
  return(diff);
}

Image *SQRT(Image *img)
{
  Image *sqrimg=NULL;
  int p,n;
  
  n = img->ncols*img->nrows;
  sqrimg = CreateImage(img->ncols,img->nrows);
  for (p=0; p < n; p++)
    sqrimg->val[p] = (int)sqrt((double)img->val[p]);
  return(sqrimg);
}

Image  *And(Image *img1,Image *img2)
{
  Image *and;
  int p,n;

  and = CreateImage(img1->ncols,img1->nrows);
  n   = img1->ncols*img1->nrows;
  for (p=0; p < n; p++)
    and->val[p] = MIN(img1->val[p],img2->val[p]);
  return(and);
}

Image  *Or(Image *img1,Image *img2)
{
  Image *or;
  int p,n;

  or = CreateImage(img1->ncols,img1->nrows);
  n   = img1->ncols*img1->nrows;
  for (p=0; p < n; p++)
    or->val[p] = MAX(img1->val[p],img2->val[p]);
  return(or);
}

Image  *Mult(Image *img1,Image *img2)
{
  Image *mult;
  int p,n;

  mult = CopyImage(img1);
  n = img1->ncols*img1->nrows;
  for (p=0; p < n; p++)
    mult->val[p] *= img2->val[p];
  return(mult);
}

Image *Add(Image *img, int value)
{
  Image *add;
  int Imax,p,n;

  add  = CopyImage(img);
  Imax = MaximumValue(img);
  n = img->ncols*img->nrows;
  for (p=0; p < n; p++)
    add->val[p] = MIN(add->val[p]+value,Imax);
  return(add);
}

Image *Sub(Image *img, int value)
{
  Image *sub;
  int p,n;

  sub = CopyImage(img);
  n = img->ncols*img->nrows;
  for (p=0; p < n; p++)
    sub->val[p] = MAX(sub->val[p]-value,0);
    
  return(sub);
}

Image *Complement(Image *img)
{
  Image *cimg=NULL;
  int p,n,Imax;
  
  n = img->ncols*img->nrows;
  cimg = CreateImage(img->ncols,img->nrows);
  Imax = MaximumValue(img);
  for (p=0; p < n; p++)
    cimg->val[p] = Imax - img->val[p];
  return(cimg);
}

Image *Abs(Image *img)
{
  Image *absimg=NULL;
  int p,n;
  
  n = img->ncols*img->nrows;
  absimg = CreateImage(img->ncols,img->nrows);
  for (p=0; p < n; p++)
    absimg->val[p] = abs(img->val[p]);
  return(absimg);
}

Image *GaussStretch(Image *img, float mean, float stdev)
{
  float *gauss=NULL,sq,var2;
  int i,Imax,n;
  Image *gimg=NULL;

  Imax  = MaximumValue(img);
  gauss = AllocFloatArray(Imax+1);
  var2  = 2*stdev*stdev;
  for (i=0; i < Imax+1; i++){
    sq  = ((float)i-mean)*((float)i-mean);
    gauss[i]=(float)(Imax*exp(-sq/var2));
  }
  n     = img->ncols*img->nrows;
  gimg     = CreateImage(img->ncols,img->nrows);
  for (i=0; i < n; i++){
    gimg->val[i] = (int)(gauss[img->val[i]]);
  }
  free(gauss);
  return(gimg);  
}

Image *DrawBorder(Image *img, Image *label, int value)
{
  Image *himg=NULL;
  int p,q,i;
  AdjRel *A=NULL;
  Pixel u,v;

  himg = CopyImage(img);
  A    = Circular(1.0);
  for (u.y=0; u.y < himg->nrows; u.y++){
    for (u.x=0; u.x < himg->ncols; u.x++){
      p = u.x + himg->tbrow[u.y];
      for (i=1; i < A->n; i++){
	v.x = u.x + A->dx[i];
	v.y = u.y + A->dy[i];
	if (ValidPixel(himg,v.x,v.y)){
	  q = v.x + himg->tbrow[v.y];
	  if (label->val[p] < label->val[q]){
	    himg->val[p] = value;
	    break;
	  }
	}
      }
    }
  }
  DestroyAdjRel(&A);

  return(himg);
}

void DrawPoint(Image *img, int x, int y, float radius, int value)
{
  Pixel v;
  int r0,r2,dx,dy,p;

  r0 = (int)radius;
  r2 = (int)(radius*radius);
  for(v.y = (y-r0), dy=-r0; v.y <= (y+r0); v.y++,dy++)
    for(v.x = (x-r0), dx=-r0; v.x <= (x+r0); v.x++,dx++)
      if (ValidPixel(img,v.x,v.y)&&(((dx*dx)+(dy*dy)) <= r2)){
	p = v.x + img->tbrow[v.y];
	img->val[p] = value;
      }
}

void DrawPath(Image *img,Image *pred, int dst,int value)
{
  int p;
  
  p = dst;
  while(pred->val[p]!=NIL){
    img->val[p]=value;
    p = pred->val[p];
  }
  img->val[p]=value;

}

/* jump # */
int nc_fgets(char *s, int m, FILE *f) {
  while(fgets(s,m,f)!=NULL)
    if (s[0]!='#') return 1;
  return 0;
}

