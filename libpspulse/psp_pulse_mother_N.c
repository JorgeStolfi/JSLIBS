/* See psp_pulse_mother_N.h */
/* Last edited on 2011-09-18 13:51:19 by stolfilocal */

#include <math.h>

#include <affirm.h>

#include <psp_basic.h>
#include <psp_pulse.h>
#include <psp_pulse_mother_N.h>

/* INTERNAL PROTOTYPES */

void psp_pulse_mother_N_main_bezier
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix,
    double *bz
  );
  
/* IMPLEMENTATIONS */

unsigned int psp_pulse_mother_N_count (psp_cont_t c, psp_degree_t g)
  { affirm(g >= 0, "negative degree");
    affirm(c >= -1, "bad continuity class");
    if (g <= 2*c)
      { return 0; }
    else
      { return g-c; }
  }

psp_grid_size_t psp_pulse_mother_N_supp_count
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix
  )
  { affirm(c >= -1, "bad contin for N-pulse");
    affirm(g >= 0, "negative degree for N-pulse");
    affirm(g >= 2*c + 1, "bad degree for N-pulse");
    affirm((pix >= 0) && (pix <= g-c-1), "invalid index");
    return (pix > c ? 1 : 2);
  }

void psp_pulse_mother_N_to_bezier
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t x,
    double *bz
  )
  { affirm(c >= -1, "bad contin for N-pulse");
    affirm(g >= 0, "negative degree for N-pulse");
    affirm(g >= 2*c + 1, "bad degree for N-pulse");
    affirm(g <= psp_pulse_MAX_DEGREE, "degree too high for N-pulse");
    affirm((pix >= 0) && (pix <= g-c-1), "invalid index");
    int j;
    for (j = 0; j <= g; j++) { bz[j] = 0; }
    if (pix > c)
      { /* A Bernstein-Bézier element with single-cell support: */
        if (x == 0) { bz[pix] = 1; }
      }
    else
      { /* A two-interval pulse: */
        int r, s;
        if (x == 0)
          { int xip = g-c+pix;
            bz[xip] = 1.0;
            for (r = xip+1; r <= g; r++) { bz[r] = bz[r-1]/2; }
            for (s = pix; s > 0; s--)
              { bz[xip] /= 2;
                for (r = xip+1; r <= g; r++) { bz[r] = (bz[r] + bz[r-1])/2; }
              }
          }
        else if (x == 1) 
          { bz[pix] = 1.0;
            for (r = pix-1; r >= 0; r--) { bz[r] = bz[r+1]/2; }
            for (s = c-pix; s > 0; s--)
              { bz[pix] = bz[pix]/2;
                for (r = pix-1; r >= 0; r--) { bz[r] = (bz[r] + bz[r+1])/2; }
              }
          }
     }
  }
  
/* SPECIALIZED PROCS FOR MAIN PATCH: */

void psp_pulse_mother_N_main_bezier
  ( psp_cont_t c, 
    psp_degree_t g, 
    psp_pulse_mother_index_t pix,
    double *bz
  )
  { int r, s;
    bz[pix] = 1.0;
    for (r = pix-1; r >= 0; r--) { bz[r] = bz[r+1]/2; }
    for (s = pix+1; s <= c; s++)
      { bz[pix] = bz[pix]/2;
        for (r = pix-1; r >= 0; r--) { bz[r] = (bz[r] + bz[r+1])/2; }
      }
  }
