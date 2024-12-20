/* See bbopt.h */
/* Last edited on 2024-12-05 10:21:28 by stolfi */

#include <stdio.h>
#include <stdlib.h>

#include <affirm.h>
#include <ia.h>
#include <aa.h>

#include <fbox.h>
#include <fboxheap.h>
#include <fboxlist.h>

#include <bbopt.h>

#define MaxPieces 2

void bb_partition(FBox *b, Float *tol, Interval F(Interval *xr), FBox **c, int32_t *ncp);
  /* Partitions the box {b} into two or more sub-boxes {c[0..n-1]} and
    sets {*ncp} to the count {n} of those sub-boxes. The new boxes will 
    not share any storage with {b}, and their ranges are estimated with the
    {F} function.
    
    If the box {b} has all sides smaller than the corresponding
    tolerance, or cannot be split due to rounding errors, the
    procedure returns {*ncp = 0} and allocates no new storage.
    
    Currently the split is perpendicular to the axis {i}
    such that {b[i]/tol[i]} is largest. */
    
sign bb_L_compare(FBox *a, FBox *b);
  /* Compares boxes according to splitting priority. */
    
sign bb_Q_compare(FBox *a, FBox *b);
  /* Compares boxes according to lower bound {b.fr.lo} */

FBoxList bb_optimize
  ( int32_t d, 
    Interval F(Interval *xr),
    Interval *xr,
    Float *tol,
    Interval *fm,
    void report(Interval *xr, Interval fr, bool_t final)
  )
  { Float fmhi = PlusInfinity; /* Upper bound for global minimum */
    FBoxHeap *L = fbox_heap_new(100, bb_L_compare);
    FBoxHeap *Q = fbox_heap_new(100, bb_Q_compare);
    int32_t nc; 
    
    FBox *c[MaxPieces]; /* Return area for {partition()} */
    FBox *b = fbox_make(d, 0, xr, F(xr));
    fbox_heap_insert(L, b);

    fprintf(stderr, "initial box: \n");
    fbox_print(stderr, b);
    
    while(L->n > 0)
      { FBox *b = fbox_heap_pop(L);
        report(&(b->xr[0]), b->fr, FALSE);
        if (b->fr.lo > fmhi) 
          { /* ignore */
            fbox_discard(b);
          }
        else
          { int32_t i;
            bb_partition(b, tol, F, &(c[0]), &nc);
            affirm(nc <= MaxPieces, "oops, too many pieces");
            for(i = 0; i < nc; i++)
              { if (c[i]->fr.hi < fmhi) 
                  { fmhi = c[i]->fr.hi;
                    /* Flush bad boxes from output queue: */
                    while((Q->n > 0) && (Q->b[0]->fr.lo > fmhi))
                      { fbox_heap_pop(Q); }
                  }
              }
            if (nc == 0)
              { if (b->fr.lo <= fmhi) { fbox_heap_insert(Q, b); } }
            else
              { fbox_discard(b);
                for(i = 0; i < nc; i++)
                  { if (c[i]->fr.lo <= fmhi) { fbox_heap_insert(L, c[i]); } }
              }
          }
      }
    { Float fmlo = PlusInfinity; /* Lower bound for global minimum */
      FBoxList R = NULL;
      while (Q->n > 0)
        { FBox *b = fbox_heap_pop(Q); 
          report(&(b->xr[0]), b->fr, TRUE);
          if (b->fr.lo < fmlo) { fmlo = b->fr.lo; }
          R = fbox_cons(b, R);
        }
      (*fm) = (Interval){fmlo, fmhi};
      return R;
    }
  }
       
void bb_partition(FBox *b, Float *tol, Interval F(Interval *xr), FBox **c, int32_t *ncp)
  { int32_t i;
    int32_t imax = -1;
    Float wmax = -1.0, mmax = NAN;
    Interval *xr = &(b->xr[0]);
    for (i = 0; i < b->d; i++)
      { Interval xri = xr[i];
        Float wi = (xri.hi - xri.lo)/tol[i];
        Float mi = (xri.hi + xri.lo)/2;
        if ((wi > wmax) && (mi > xri.lo) && (mi < xri.hi))
          { wmax = wi; imax = i; mmax = mi; }
      }
    if (wmax <= 1.0)
      { /* Box cannot be split */
        (*ncp) = 0;
      }
    else
      { Interval zr = xr[imax];
        xr[imax].hi = mmax;
        c[0] = fbox_make(b->d, b->depth + 1, xr, F(xr));
        xr[imax].hi = zr.hi;
        
        xr[imax].lo = mmax;
        c[1] = fbox_make(b->d, b->depth + 1, xr, F(xr));
        xr[imax].lo = zr.lo;
        (*ncp) = 2;
      }
  }

sign bb_L_compare(FBox *a, FBox *b)
  { Float ma = (a->fr.lo + a->fr.hi)/2;
    Float mb = (b->fr.lo + b->fr.hi)/2;
    if (ma < mb) 
      { return -1; }
    else if (ma > mb)
      { return +1; }
    else
      { return 0; }
  }

sign bb_Q_compare(FBox *a, FBox *b)
  { Float ma = a->fr.lo;
    Float mb = b->fr.lo;
    if (ma > mb) 
      { return -1; }
    else if (ma > mb)
      { return +1; }
    else
      { return 0; }
  }
