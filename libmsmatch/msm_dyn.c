/* See {msm_dyn.h} */
/* Last edited on 2022-10-20 07:58:57 by stolfi */

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
    tb.sMin = int32_vec_new(0);
    tb.sMax = int32_vec_new(0);
    tb.ns = 0;
    tb.ev = msm_dyn_entry_vec_new(0);
    return tb;
  }
    
void msm_dyn_tableau_resize(msm_dyn_tableau_t *tb, int32_t rMin, int32_t rMax, int32_t maxds)
   { 
    if(rMin > rMax)
      { fprintf(stderr, "Houston, we have a problem...\n");
        fprintf(stderr, "\nrMin = %d, rMax = %d\n", rMin, rMax);
      }
    demand(rMin <= rMax, "bad {r} range");
    tb->rMin = rMin;
    tb->rMax = rMax;
    tb->nr = rMax - rMin + 1;
    int32_vec_expand(&(tb->sMin), tb->nr - 1);
    int32_vec_expand(&(tb->sMax), tb->nr - 1);
    int32_t dr;
    for (dr = 0; dr < tb->nr; dr++) { tb->sMin.e[dr] = +1; tb->sMax.e[dr] = -1; }
    demand(maxds % 2 == 0, "maxds should be even");
    tb->ns = maxds/2 + 1;
    msm_dyn_entry_vec_expand(&(tb->ev), tb->nr*tb->ns - 1);
  }
    
void msm_dyn_tableau_set_s_range(msm_dyn_tableau_t *tb, int32_t r, int32_t sMin, int32_t sMax)
  { /* Get index {dr} into {tb->sMin,tb->sMax}: */
    demand((tb->rMin <= r) && (r <= tb->rMax), "invalid {r}");
    int32_t dr = r - tb->rMin;
    /* Save the {s} range: */
    tb->sMin.e[dr] = sMin;
    tb->sMax.e[dr] = sMax;
    if (sMin <= sMax)
      { /* Check parity of {s} values: */
        demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
        demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
        /* Compute the number {ne} of entries for this {r}: */
        int32_t ne = (sMax - sMin)/2 + 1;
        demand(ne <= tb->ns, "{s} range too large");
        /* Reset all tableau entries for this {r} value: */
        int32_t ds;
        for (ds = 0; ds < ne; ds++)
          { int32_t k = dr*tb->ns + ds;
            tb->ev.e[k] = (msm_dyn_entry_t){ -INF, msm_rung_none };
          }
      }
  }
    
void msm_dyn_tableau_get_s_range(msm_dyn_tableau_t *tb, int32_t r, int32_t *sMin, int32_t *sMax)
  { if ((r < tb->rMin) || (r > tb->rMax))
      { /* {R} outside the range, return an empty interval: */
        (*sMin) = +1; (*sMax) = -1;
      }
    else
      { /* Return the stored {s}-interval: */
        (*sMin) = tb->sMin.e[r - tb->rMin];
        (*sMax) = tb->sMax.e[r - tb->rMin];
      }
  }

void msm_dyn_tableau_get_i0_i1_ranges(msm_dyn_tableau_t *tb, int32_t *i0Min, int32_t *i0Max, int32_t *i1Min, int32_t *i1Max)
  {
    /* Start with empty intervals: */
    (*i0Min) = +1; (*i0Max) = -1;
    (*i1Min) = +1; (*i1Max) = -1;
    int32_t r;
    for (r = tb->rMin; r <= tb->rMax; r++)
      { /* Get the {s} range corresponding to this {r} value: */
        int32_t sMin, sMax;
        msm_dyn_tableau_get_s_range(tb, r, &sMin, &sMax);
        if (sMin <= sMax)
          { demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
            demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
            /* We need only consider the extremal cases: */
            int32_t j;
            for (j = 0; j <= 1; j++)
              { int32_t s = (j == 0 ? sMin : sMax);
                /* Convert {r,s} to {i0,i1}: */
                int32_t i0 = (r + s)/2;
                int32_t i1 = (r - s)/2;
                /* Expand X and Y ranges to include {(i0,i1)}: */
                if ((*i0Min) > (*i0Max)) 
                  { (*i0Min) = (*i0Max) = i0; }
                else if (i0 < (*i0Min))
                  { (*i0Min) = i0; }
                else if (i0 > (*i0Max))
                  { (*i0Max) = i0; }
                if ((*i1Min) > (*i1Max)) 
                  { (*i1Min) = (*i1Max) = i1; }
                else if (i1 < (*i1Min))
                  { (*i1Min) = i1; }
                else if (i1 > (*i1Max))
                  { (*i1Max) = i1; }
              }
          }
      }
  }

msm_dyn_entry_t *msm_dyn_tableau_get_entry_address(msm_dyn_tableau_t *tb, msm_rung_t g)
  { /* Convert {g} to the R,S coordinate system: */
    int32_t r = g.c[0] + g.c[1];
    int32_t s = g.c[0] - g.c[1];
    /* Find index {dr} into {sMin,sMax}: */
    if ((r < tb->rMin) || (r > tb->rMax)) { return NULL; }
    int32_t dr = r - tb->rMin;
    /* Get relative S-coordinate {ds}: */
    int32_t sMin = tb->sMin.e[dr];
    int32_t sMax = tb->sMax.e[dr];
    if ((s < sMin) || (s > sMax)) { return NULL; }
    assert((s - sMin) % 2 == 0);
    int32_t ds = (s - sMin)/2;
    int32_t k = dr*tb->ns + ds;
    return &(tb->ev.e[k]);
  }

void msm_dyn_tableau_free(msm_dyn_tableau_t *tb)
  { free(tb->sMin.e);
    free(tb->sMax.e);
    free(tb->ev.e);
  }

