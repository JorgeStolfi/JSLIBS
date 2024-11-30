/* See {baltuple.h}, */
/* Last edited on 2024-11-15 19:11:21 by stolfi */

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <baltuple.h>

void balt_first_tuple(int32_t n, int32_t d[], int32_t vmin, int32_t vmax, int32_t smin, int32_t smax)
  { int32_t k = 0;
    while (k < n)
      { /* Now {smin..smax} is the range for the sum {d[k..n-1]}: */
        /* Find range {dmin..dmax} for {d[k]}: */ 
        int32_t dmin = (int32_t)imax(vmin, smin);  /* Min valid value. */
        int32_t dmax = (int32_t)imin(vmax, smax);  /* Max valid value. */
        demand(dmin <= dmax, "no valid tuple");
        /* Pick first valid value {dk} in range: */
        int32_t dk = balt_ix_first(dmin, dmax);
        assert((dmin <= dk) && (dk <= dmin));
        /* Set {d[k]} to {dk}: */
        d[k] = dk;
        /* Update the sum range {smin..smax} for the sum {d[k+1..n-1]}: */
        smin = smin - dk; smax = smax - dk;
        /* Now {smin..smax} is the range for the sum {d[k+1..n-1]}: */
        k--;
       /* Now {smin..smax} is the range for the sum {d[k..n-1]}: */
       }     
  }

bool_t balt_next_tuple(int32_t n, int32_t d[], int32_t vmin, int32_t vmax, int32_t smin, int32_t smax)
  { 
    /* Update {smin,smax} to range of {d[n]}, it it existed: */
    for (uint32_t i = 0;  i < n; i++) { int32_t e = d[i];  smin = smin - e; smax = smax - e; }
    demand((smin <= 0) && (0 <= smax), "bad tuple sum");
    
    /* Find the last {d[k]} in {d[0..n-1]} that can be bumped: */
    int32_t k = n;
    while (TRUE)
      { /* Now {smin..smax} is the range of values for {d[k]} allowed by the sum bound. */
        if (k <= 0) { return TRUE; }
        /* Back up one entry, update {smin,smax}. */
        k--;
        smin = smin + (d[k] + vmax);
        smax = smax + (d[k] + vmin);
        /* Now {smin..smax} is the sum bound on {d[k]}. */
        /* Find the range {[dmin..dmax]} of valid cands to {d[k]}: */ 
        int32_t dmin = (int32_t)imax(vmin, smin);  /* Min valid value. */
        int32_t dmax = (int32_t)imin(vmax, smax);  /* Max valid value. */
        assert(dmin <= dmax);
        assert((dmin <= d[k]) && (d[k] <= dmax)); /* The current value must be valid. */
        /* Find the next valid value {dk} of coord {d[k]} in the search sequence: */
        int32_t dk = balt_ix_next(d[k], dmin, dmax);
        if ((dmin <= dk) && (dk <= dmax)) 
          { /* Bump {d[k]} to {dk}: */
            d[k] = dk;
            smin = smin - (dk + vmax);
            smax = smax - (dk + vmin);
            /* Reset all {d[k..n-1]} to their darliest values: */
            balt_first_tuple(n - k, &(d[k]), vmin, vmax, smin, smax);
            return FALSE;
          }
        else
          { /* Back up to the previous elememt. */ }
      }
  } 

int32_t balt_ix_first(int32_t vmin, int32_t vmax)
  {
    if (vmin > 0)
      { return vmin; }
    else if (vmax < 0)
      { return vmax; }
    else
      { return 0; }
  }

int32_t balt_ix_next(int32_t v, int32_t vmin, int32_t vmax)
  {
    if (vmin > 0)
      { return (v < vmin ? vmin : v + 1); }
    else if (vmax < 0)
      { return (v > vmax ? vmax : v - 1); }
    else
      { /* Range {vmin..vmax} contains zero. Try the next value from {SEQ}: */
        if (v < 0)
          { v = -v; }
        else 
          { v = -(v + 1); }
        if ((v >= vmin) && (v <= vmax)) { return v; }
        /* Try again with the next value: */
        if (v < 0)
          { v = -v; }
        else 
          { v = -(v + 1); }
        if ((v >= vmin) && (v <= vmax)) { return v; }
        /* Two successive failures; it means exhaustion: */
        return (vmax > -vmin ? 1 + vmax : 1 - vmin);
      }
  }

