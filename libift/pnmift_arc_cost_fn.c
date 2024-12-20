/* pnmift_arc_cost_fn.c - implementation of pnmift_arc_cost_fn.h */
/* Last edited on 2024-12-05 10:29:10 by stolfi */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <affirm.h>
#include <bool.h>
#include <frgb.h>
#include <frgb_ops.h>
  
#include <ift.h>
#include <pnmift_arc_cost_fn.h>

pnmift_arc_cost_fn_t *pnmift_arc_cost_fn_from_name(char *name)
  {
    if (strcmp(name, "ediff_rgb") == 0)
      { return &pnmift_arc_cost_fn_ediff_rgb; }
    else if (strcmp(name, "ediff_yuv") == 0)
      { return &pnmift_arc_cost_fn_ediff_yuv; }
    else if (strcmp(name, "lum") == 0)
      { return &pnmift_arc_cost_fn_lum; }
    else 
      { demand(FALSE, "bad function name"); }
    return NULL;
  }

pnmift_arc_cost_t pnmift_arc_cost_fn_ediff_rgb(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns)
  {
    if (chns == 1)
      { return fabs(q.c[0] - p.c[0]); }
    else if (chns == 3)
      { double sum_d2 = 0;
        int chn;
        for (chn = 0; chn < chns; chn++)
          { double d = ((double)p.c[chn]) - ((double)q.c[chn]);
            sum_d2 += d*d;
          }
        return sqrt(sum_d2/chns);
      }
    else
      { demand(FALSE, "bad channel count"); }
  }

pnmift_arc_cost_t pnmift_arc_cost_fn_ediff_yuv(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns)
  {
    frgb_YUV_to_yuv(&p, 0.02);
    frgb_YUV_to_yuv(&q, 0.02);
    return pnmift_arc_cost_fn_ediff_rgb(p, ra, q, chns);
  }

pnmift_arc_cost_t pnmift_arc_cost_fn_lum(frgb_t p, ift_rel_arc_t *ra, frgb_t q, int chns)
  {
    return frgb_get_Y(&q);
  }
