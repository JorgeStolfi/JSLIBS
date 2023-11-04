/* See wt_median.h */
/* Last edited on 2023-11-04 00:13:19 by stolfi */

#define wt_median_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsmath.h>

#include <wt_median.h>

int32_t wt_median_index_set_update(int32_t nt, int32_t tx[], int32_t nx, int32_t ix, int32_t nw);
  /* Assumes that, on entry, {tx[0..nt-1]} is a set of {nt} consecutive
    indices in {0..nx-1}. On exit, {tx[0..nw-1]} will contain the {nw}
    consecutive indices centered at {ix} namely {ix-hw..ix+hw} where
    {hw=nw/2}.
    
    The input set size {nt} may be bigger or smaller than the window
    width {nw}, which must be odd.  In any case, the indices of the
    output set which are in the input set are merely copied, preserving
    their order.  The procedure returns the number of new indices, which
    will be appended at the end.  The range {ix-hw..ix+hw} must be a sub-range
    of {0..nx-1}. */

void wt_median_index_set_full_sort(int32_t nx, double x[], int32_t nw, int32_t tx[]);
void wt_median_index_set_insertion_sort(int32_t nx, double x[], int32_t nw, int32_t tx[]);
  /* These procedures sort the indices {tx[0..nw-1]}, a subset of {0..nx-1}, so that 
    the elements {x[tx[i]]} are in non-decreasing order when {i} ranges from 0 to {nw-1}.
    
    The insertion sort version may be faster if the list {tx[0..nw-1]}
    is already almost sorted. Otherwise the full sort version (which
    uses {qsort}) is probably faster. */


/* IMPLEMENTATIONS */

double wt_median(int32_t nx, double x[], int32_t ix, int32_t nw, double w[], bool_t interp, int32_t *nk_P, int32_t kx[])
  {
    bool_t debug = TRUE;
    
    demand((nw >= 1) && ((nw % 2) == 1), "invalid window width {nw}");
    int32_t hw = nw/2;
    demand((ix-hw >= 0) && (ix+hw < nx), "window spills outside {0..nx-1}");
    
    if (nw == 1)
      { /* Trivial median: */
        if (kx != NULL) { kx[0] = ix; (*nk_P) = 1; }
        return x[ix];
      }

    /* Get the set {tx[0..nw-1]} of indices of elems in window, sorted by {x} value: */
    int32_t *tx = NULL;
    { int32_t nt;
      if (kx == NULL)
        { nt = 0; tx = talloc(nw, int32_t); }
      else
        { demand(nk_P != NULL, "invalid {nk_P}");
          nt = (*nk_P);
          tx = kx;
        }
      int32_t nt_add = wt_median_index_set_update(nt, tx, nx, ix, nw);
      if (nt_add < 5)
        { wt_median_index_set_insertion_sort(nx, x, nw, tx); }
      else
        { wt_median_index_set_full_sort(nx, x, nw, tx); }
    }
    
    /* Find the highest {ia} such that the sum of weights in {0..ia} is at most 0.5: */
    int32_t ia = -1; double suma = 0;
    while (ia < nw-1)
      { double wa = w[hw + tx[ia+1]-ix];
        if (suma + wa > 0.5) { break; }
        ia++; suma += wa;
      }
    
    /* Find the lowest {ib} such that the sum of weights in {ib..nw-1} is at most 0.5: */
    int32_t ib = nw; double sumb = 0;
    while (ib >= 1)
      { double wb = w[hw + tx[ib-1]-ix];
        if (sumb + wb > 0.5) { break; }
        ib--; sumb += wb;
      }
    demand((suma <= 0.5) && (sumb <= 0.5), "weights add to more than 1.");
      
    if (ia > ib) 
      { /* The indices {ia} and {ib} should not have crossed. */
        /* Either there were roundoff errors or the weights did not quite add to 1. */
        /* Fix by force: */
        if (debug) { fprintf(stderr, "  ia = %d suma = %24.15e  ib = %d sumb = %24.15e", ia, suma, ib, sumb); }
        int32_t im = (ia + ib)/2; 
        if (im < 0) 
          { im = 0; }
        else if (im >= nw)
          { im = nw-1; }
        else
          { im = (ia + ib + (im % 2))/2; /* Round to even if {ia+ib} is odd. */ }
        ia = im; ib = im; 
        if (debug) { fprintf(stderr, " fixed to %d\n", ia); }
        suma = sumb = 0.5; /* Fake it, just in case. */
      }
      
    assert((0 <= ia) && (ia <= ib) && (ib < nw));
    double xm;
    if (ib == ia)
      { /* Perfect median: */
        xm = x[tx[ia]];
      }
    else if (ib == ia+1)
      { double xa = x[tx[ia]], wa = w[hw + tx[ib]-ix];
        double xb = x[tx[ib]], wb = w[hw + tx[ib]-ix];
        /* Median is in the range {[xa_xb]}: */
        assert(xa < xb);
        double da = fabs(sumb - suma + wa); /* What {|sumb-suma|} would be if {xa} is chosen. */ 
        double db = fabs(sumb - suma - wb); /* What {|sumb-suma|} would be if {xb} is chosen. */
        assert(db + da > 0);
        if (interp)
          { double f = da/(db + da);
            xm = (1-f)*xa + f*xb;
          }
        else
          { /* Median is either {ib} or {ia}: */
            if (da < db)
              { xm = xa; }
            else if (db < da) 
              { xm = xb; }
            else
              { /* Choose the one with even window index: */
                xm = ((ib % 2) == 0 ? xb : xa);
              }
          }
      }
    else
      { /* The median is between {ib} and {ia}: */
        assert(ib - ia >= 2);
        double xmb = x[tx[ib-1]], xma = x[tx[ia+1]];
        assert(xmb == xma);
        xm = xmb;
      }
      
    /* Return {nt} or cleanup: */
    if (tx == kx)
      { (*nk_P) = nw; }
    else
      { free(tx); }
    return xm;
  }
    
int32_t wt_median_index_set_update(int32_t nt, int32_t tx[], int32_t nx, int32_t ix, int32_t nw)
  { bool_t debug = TRUE;
  
    demand((nw >= 1) && ((nw % 2) == 1), "invalid window width {nw}");
    int32_t hw = nw/2;
    demand((ix-hw >= 0) && (ix+hw < nx), "window spills outside {0..nx-1}");
    
    demand((nt >= 0) && (nt <= nx), "invalid index set size {nt}");
    int32_t nt_add = -1; /* Indices that were added to {tx[0..nt-1]}. */
    if (nt == 0)
      { if (debug) { fprintf(stderr, "  index set was empty\n"); }
        for (int32_t k = -hw; k <= +hw; k++) { tx[nt] = ix + k; nt++; }
        nt_add = nw;
      }
    else
      { int32_t tx_min = INT32_MAX, tx_max = INT32_MIN;  /* Min and max valid indices in {tx[0..nt-1]}. */

        /* Remove from {tx[0..nt-1]} the indices that are outside the window: */
        /* Also find min and max indices {tx_min,tx_max} in remaining set: */
        int32_t mt = 0;
        for (int32_t jt = 0; jt < nt; jt++)
          { int32_t txj = tx[jt];
            demand((txj >= 0) && (txj < nx), "invalid index in {kx[..]}");
            if ((txj >= ix-hw) && (txj <= ix+hw))
              { tx[mt] = txj; mt++;
                if (txj < tx_min) { tx_min = txj; }
                if (txj > tx_max) { tx_max = txj; }
              }
          }
        if (debug) { fprintf(stderr, "  removed %d indices, left %d", nt - mt, mt); }
        if (mt != 0) 
          { if (debug) { fprintf(stderr, "  = {%d..%d}", tx_min, tx_max); }
            demand(mt == tx_max - tx_min + 1, "indices {kx[..]} are not consecutive");
          }
        if (debug) { fprintf(stderr, "\n"); }
        demand(mt <= nw, "there are repeated indices in {kx[..]}");
        nt = mt;
        
        /* Complete {tx[0..nt-1]} to the set of all window indices: */
        nt_add = nw - nt;
        if (nt == 0)
          { for (int32_t k = -hw; k <= +hw; k++) { tx[nt] = ix + k; nt++; } }
        else
          { /* Complete the set of window indices: */
            assert((tx_min >= ix-hw) && (tx_max <= ix+hw));
            if (debug && (tx_min > ix-hw)) { fprintf(stderr, "  adding {%d..%d}\n", ix-hw, tx_min-1); }
            while(tx_min > ix-hw) { tx_min--; tx[nt] = tx_min; nt++; }
            if (debug && (tx_max < ix+hw)) { fprintf(stderr, "  adding {%d..%d}\n", tx_max+1, ix+hw); }
            while(tx_max < ix+hw) { tx_max++; tx[nt] = tx_max; nt++; }
            assert(nt == nw);
          }
      }
    assert(nt_add >= 0);
    return nt_add;
  }
  
void wt_median_index_set_insertion_sort(int32_t nx, double x[], int32_t nw, int32_t tx[])
  { for (int32_t i = 1; i < nw; i++)
      { int32_t txi = tx[i];
        int32_t j = i;
        while ((j > 0) && (x[tx[j-1]] > x[txi])) { tx[j] = tx[j-1]; j--; }
        tx[j] = txi;
      }
  }
  
void wt_median_index_set_full_sort(int32_t nx, double x[], int32_t nw, int32_t tx[])
  {
    auto int32_t compx(const void *a, const void *b);
    
    qsort(tx, nw, sizeof(int32_t), compx);
    return;
    
    int32_t compx(const void *a, const void *b)
      { int32_t *ia = (int32_t *)a; assert((*ia) >= 0 && ((*ia) < nx));
        int32_t *ib = (int32_t *)b; assert((*ib) >= 0 && ((*ib) < nx));
        double xa = x[*ia], xb = x[*ib];
        if (xa < xb) 
          { return -1; }
        else if (xa > xb)
          { return +1; }
        else
          { return 0; }
      }
  }     
