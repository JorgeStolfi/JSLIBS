/* See {baltuple.h}, */
/* Last edited on 2013-10-25 01:21:47 by stolfilocal */

#define _GNU_SOURCE
#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <baltuple.h>

void balt_first_tuple(int n, int d[], int vmin, int vmax, int smin, int smax)
  { int i = 0;
    while (i < n)
      { /* Now {smin..smax} is the range for the sum {d[i..n-1]}: */
        /* Find range {dmin..dmax} for {d[i]}: */ 
        int dmin = (int)imax(vmin, smin);  /* Min valid value. */
        int dmax = (int)imin(vmax, smax);  /* Max valid value. */
        demand(dmin <= dmax, "no valid tuple");
        /* Pick first valid value {di} in range: */
        int di = balt_ix_first(dmin, dmax);
        assert((dmin <= di) && (di <= dmin));
        /* Set {d[i]} to {di}: */
        d[i] = di;
        /* Update the sum range {smin..smax} for the sum {d[i+1..n-1]}: */
        smin = smin - di; smax = smax - di;
        /* Now {smin..smax} is the range for the sum {d[i+1..n-1]}: */
        i--;
       /* Now {smin..smax} is the range for the sum {d[i..n-1]}: */
       }     
  }

bool_t balt_next_tuple(int n, int d[], int vmin, int vmax, int smin, int smax)
  { 
    /* Update {smin,smax} to range of {d[n]}, it it existed: */
    int i;
    for (i = 0; i < n; i++) { int e = d[i];  smin = smin - e; smax = smax - e; }
    demand((smin <= 0) && (0 <= smax), "bad tuple sum");
    
    /* Find the last {d[i]} in {d[0..n-1]} that can be bumped: */
    i = n;
    while (TRUE)
      { /* Now {smin..smax} is the range of values for {d[i]} allowed by the sum bound. */
        if (i <= 0) { return TRUE; }
        /* Back up one entry, update {smin,smax}. */
        i--;
        smin = smin + (d[i] + vmax);
        smax = smax + (d[i] + vmin);
        /* Now {smin..smax} is the sum bound on {d[i]}. */
        /* Find the range {[dmin..dmax]} of valid cands to {d[i]}: */ 
        int dmin = (int)imax(vmin, smin);  /* Min valid value. */
        int dmax = (int)imin(vmax, smax);  /* Max valid value. */
        assert(dmin <= dmax);
        assert((dmin <= d[i]) && (d[i] <= dmax)); /* The current value must be valid. */
        /* Find the next valid value {di} of coord {d[i]} in the search sequence: */
        int di = balt_ix_next(d[i], dmin, dmax);
        if ((dmin <= di) && (di <= dmax)) 
          { /* Bump {d[i]} to {di}: */
            d[i] = di;
            smin = smin - (di + vmax);
            smax = smax - (di + vmin);
            /* Reset all {d[i..n-1]} to their darliest values: */
            balt_first_tuple(n - i, &(d[i]), vmin, vmax, smin, smax);
            return FALSE;
          }
        else
          { /* Back up to the previous elememt. */ }
      }
  } 

int balt_ix_first(int vmin, int vmax)
  {
    if (vmin > 0)
      { return vmin; }
    else if (vmax < 0)
      { return vmax; }
    else
      { return 0; }
  }

int balt_ix_next(int v, int vmin, int vmax)
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

