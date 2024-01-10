/* See {msm_dyn.h} */
/* Last edited on 2008-01-11 19:07:32 by stolfi */

#define msm_dyn_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_dyn.h>
#include <msm_basic.h>
#include <msm_rung.h>

#include <vec.h>
#include <affirm.h>

#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

vec_typeimpl(msm_dyn_entry_vec_t,msm_dyn_entry_vec,msm_dyn_entry_t);
 
msm_dyn_tableau_t msm_dyn_tableau_new(void)
  { msm_dyn_tableau_t tb;
    tb.rMin = 0;
    tb.rMax = -1;
    tb.nr = 0;
    tb.sMin = int_vec_new(0);
    tb.sMax = int_vec_new(0);
    tb.ns = 0;
    tb.ev = msm_dyn_entry_vec_new(0);
    return tb;
  }
    
void msm_dyn_tableau_resize(msm_dyn_tableau_t *tb, int rMin, int rMax, int maxds)
  { demand(rMin <= rMax, "bad {r} range");
    tb->rMin = rMin;
    tb->rMax = rMax;
    tb->nr = rMax - rMin + 1;
    int_vec_expand(&(tb->sMin), tb->nr - 1);
    int_vec_expand(&(tb->sMax), tb->nr - 1);
    int dr;
    for (dr = 0; dr < tb->nr; dr++) { tb->sMin.el[dr] = +1; tb->sMax.el[dr] = -1; }
    demand(maxds % 2 == 0, "maxds should be even");
    tb->ns = maxds/2 + 1;
    msm_dyn_entry_vec_expand(&(tb->ev), tb->nr*tb->ns - 1);
  }
    
void msm_dyn_tableau_set_s_range(msm_dyn_tableau_t *tb, int r, int sMin, int sMax)
  { /* Get index {dr} into {tb->sMin,tb->sMax}: */
    demand((tb->rMin <= r) && (r <= tb->rMax), "invalid {r}");
    int dr = r - tb->rMin;
    /* Save the {s} range: */
    tb->sMin.el[dr] = sMin;
    tb->sMax.el[dr] = sMax;
    if (sMin <= sMax)
      { /* Check parity of {s} values: */
        demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
        demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
        /* Compute the number {ne} of entries for this {r}: */
        int ne = (sMax - sMin)/2 + 1;
        demand(ne <= tb->ns, "{s} range too large");
        /* Reset all tableau entries for this {r} value: */
        int ds;
        for (ds = 0; ds < ne; ds++)
          { int k = dr*tb->ns + ds;
            tb->ev.el[k] = (msm_dyn_entry_t){ -INF, msm_rung_none };
          }
      }
  }
    
void msm_dyn_tableau_get_s_range(msm_dyn_tableau_t *tb, int r, int *sMin, int *sMax)
  { if ((r < tb->rMin) || (r > tb->rMax))
      { /* {R} outside the range, return an empty interval: */
        (*sMin) = +1; (*sMax) = -1;
      }
    else
      { /* Return the stored {s}-interval: */
        (*sMin) = tb->sMin.el[r - tb->rMin];
        (*sMax) = tb->sMax.el[r - tb->rMin];
      }
  }

void msm_dyn_tableau_get_X_Y_ranges(msm_dyn_tableau_t *tb, int *xMin, int *xMax, int *yMin, int *yMax)
  {
    /* Start with empty intervals: */
    (*xMin) = +1; (*xMax) = -1;
    (*yMin) = +1; (*yMax) = -1;
    int r;
    for (r = tb->rMin; r <= tb->rMax; r++)
      { /* Get the {s} range corresponding to this {r} value: */
        int sMin, sMax;
        msm_dyn_tableau_get_s_range(tb, r, &sMin, &sMax);
        if (sMin <= sMax)
          { demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
            demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
            /* We need only consider the extremal cases: */
            int j;
            for (j = 0; j <= 1; j++)
              { int s = (j == 0 ? sMin : sMax);
                /* Convert {r,s} to {x,y}: */
                int x = (r + s)/2;
                int y = (r - s)/2;
                /* Expand X and Y ranges to include {(x,y)}: */
                if ((*xMin) > (*xMax)) 
                  { (*xMin) = (*xMax) = x; }
                else if (x < (*xMin))
                  { (*xMin) = x; }
                else if (x > (*xMax))
                  { (*xMax) = x; }
                if ((*yMin) > (*yMax)) 
                  { (*yMin) = (*yMax) = y; }
                else if (y < (*yMin))
                  { (*yMin) = y; }
                else if (y > (*yMax))
                  { (*yMax) = y; }
              }
          }
      }
  }

msm_dyn_entry_t *msm_dyn_tableau_get_entry_address(msm_dyn_tableau_t *tb, msm_rung_t g)
  { /* Convert {g} to the R,S coordinate system: */
    int r = g.c[0] + g.c[1];
    int s = g.c[0] - g.c[1];
    /* Find index {dr} into {sMin,sMax}: */
    if ((r < tb->rMin) || (r > tb->rMax)) { return NULL; }
    int dr = r - tb->rMin;
    /* Get relative S-coordinate {ds}: */
    int sMin = tb->sMin.el[dr];
    int sMax = tb->sMax.el[dr];
    if ((s < sMin) || (s > sMax)) { return NULL; }
    assert((s - sMin) % 2 == 0);
    int ds = (s - sMin)/2;
    int k = dr*tb->ns + ds;
    return &(tb->ev.el[k]);
  }

void msm_dyn_tableau_free(msm_dyn_tableau_t *tb)
  { free(tb->sMin.el);
    free(tb->sMax.el);
    free(tb->ev.el);
  }

