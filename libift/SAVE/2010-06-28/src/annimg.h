#ifndef _ANNIMG_H_
#define _ANNIMG_H_

#include "image.h"
#include "adjacency.h"

typedef struct _annimg {
  Image *img;
  Image *cost;
  Image *label;
  Image *pred;
} AnnImg;

AnnImg *CreateAnnImg(Image *img);
void    DestroyAnnImg(AnnImg **aimg);
Image  *RootMap(Image *pred);
int     Root(Image *pred, int p);

#endif
