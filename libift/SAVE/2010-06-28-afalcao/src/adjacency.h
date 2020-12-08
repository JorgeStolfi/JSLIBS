#ifndef _ADJACENCY_H_
#define _ADJACENCY_H_

#include "common.h"

typedef struct _adjrel {
  int *dx;
  int *dy;
  int n;
} AdjRel;

AdjRel *CreateAdjRel(int n);
void    DestroyAdjRel(AdjRel **A);
AdjRel *Circular(float r);
AdjRel *RightSide(AdjRel *A);
AdjRel *LeftSide(AdjRel *A);


#endif
