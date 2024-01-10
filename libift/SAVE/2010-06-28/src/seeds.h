#ifndef _SEEDS_H_
#define _SEEDS_H_

#include "genift.h"

/* Lambda functions (lambda(p) >= 0 => p is a seed pixel.) */

Image  *LambdaInit(int ncols, int nrows); /* define empty seed set. */
Image  *LambdaTrivial(int ncols, int nrows);
Image  *LambdaRegMin(Image *img, AdjRel *A); 
Image  *LambdaFrame(int ncols, int nrows, int value);
Image  *LambdaContour(Image *bin);
Image  *LambdaContPixel(Image *bin, char label); /* label=1 uses
                                                    consecutive
                                                    integers, label=0
                                                    uses trivial
                                                    labeling function */
Image  *LambdaComp(Image *bin, AdjRel *A);
void    SetLambda(Image *lambda, int x, int y, int value);

/* auxiliary */

bool ValidContPoint(Image *bin, AdjRel *L, AdjRel *R, int p);

#endif
