/* fbox.h -- function boxes for branch and bound optimization */ 
/* Last edited on 2024-12-05 10:21:48 by stolfi */

#ifndef fbox_h
#define fbox_h

#include <stdint.h>

#include <ia.h>

typedef signed char sign;

typedef struct FBox  /* An axis-aligned multidimensional box. */
  { int32_t d;            /* Dimension */
    int32_t depth;        /* Subdivision depth */
    Interval fr;      /* Estimated range of function over this box. */
    Interval xr[1];   /* Coordinate ranges for this box (actually `d' of them). */
  } FBox;
  /* A function box (or f-box for short) {b} gives guaranteed bounds
    {b.fr} on the range of some function {f} within a {d}-dimensional
    axis-aligned box, the Cartesian product of the {d} intervals
    {b.xr[0..b.d-1]}. Specifically, is assumed that {f(x)} lies in
    {b.fr} for any argument point {x} such that {x[i] \in b.xr[i]} for
    all {i = 0..b.d-1}.
    
    The field {b.xr} actually has {d} elements, not 1,
    so the true size of {b} is not {sizeof(FBox)} but
    {sizeof(FBox) + (b.d-1)*sizeof(Interval)}. */

typedef sign (*FBoxCmp)(FBox *a, FBox *b); /* A box comparison function. */

FBox *fbox_make(int32_t d, int32_t depth, Interval *xr, Interval fr);
  /* Packages the given data as a box, making a copy of the 
    coordinate range list {xr}. */

void fbox_discard(FBox *b);
  /* Deallocate the box {b}. */
    
void fbox_print(FILE *f, FBox *b);
  /* Prints {b} to file {f}. */

#endif
