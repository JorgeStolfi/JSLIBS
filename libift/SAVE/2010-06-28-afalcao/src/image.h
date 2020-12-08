#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common.h"
#include "adjacency.h"

typedef struct _pixel {
  int x,y;
} Pixel;

typedef struct _image {
  int *val;
  int ncols,nrows;
  int *tbrow;
} Image;

typedef struct _kernel {
  float *val;
  int ncols,nrows;
} Kernel;


Image  *CreateImage(int ncols,int nrows);
void    DestroyImage(Image **img);
Image  *ReadImage(char *filename);
void    WriteImage(Image *img, char *filename);
int     MinimumValue(Image *img);
int     MaximumValue(Image *img);
void    SetImage(Image *img, int value);
bool    ValidPixel(Image *img, int x, int y);
Image  *CopyImage(Image *img);
Image  *Threshold(Image *img, int lower, int higher);
Image  *Diff(Image *img1,Image *img2); 
Image  *SQRT(Image *img);
Image  *And(Image *img1, Image *img2);
Image  *Or(Image *img1, Image *img2);
Image  *Mult(Image *img1,Image *img2);
Image  *Add(Image *img, int value);
Image  *Sub(Image *img, int value);
Image  *Complement(Image *img);
Image  *Abs(Image *img);
Image  *Dilate(Image *img, AdjRel *A);
Image  *Erode(Image *img, AdjRel *A);
Image  *Close(Image *img, AdjRel *A);
Image  *Open(Image *img, AdjRel *A);
Image  *MorphGrad(Image *img, AdjRel *A);
Image  *GaussStretch(Image *img, float mean, float stdev);
Image  *DrawBorder(Image *img, Image *label, int value);
void    DrawPoint(Image *img, int x, int y, float radius, int value);
void    DrawPath(Image *img,Image *pred, int dst,int value);

#endif





