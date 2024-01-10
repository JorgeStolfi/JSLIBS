/* See udg_pulse_N.h */
/* Last edited on 2009-02-10 09:02:19 by stolfi */

#include <math.h>

#include <bz_basic.h>
#include <affirm.h>

#include <udg_pulse.h>
#include <udg_pulse_N.h>

/* INTERNAL PROTOTYPES */

void udg_pulse_N_mother_main_bezier
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix,
    double *bz
  );
  
/* IMPLEMENTATIONS */

unsigned int udg_pulse_N_num_mothers (udg_cont_t c, udg_degree_t g)
  { affirm(g >= 0, "negative degree");
    affirm(c >= -1, "bad continuity class");
    if (g <= 2*c)
      { return 0; }
    else
      { return g-c; }
  }

udg_grid_size_t udg_pulse_N_mother_supp_count
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix
  )
  { affirm(c >= -1, "bad contin for N-pulse");
    affirm(g >= 0, "negative degree for N-pulse");
    affirm(g >= 2*c + 1, "bad degree for N-pulse");
    affirm((pix >= 0) && (pix <= g-c-1), "invalid index");
    return (pix > c ? 1 : 2);
  }

void udg_pulse_N_mother_to_bezier
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t x,
    double *bz
  )
  { affirm(c >= -1, "bad contin for N-pulse");
    affirm(g >= 0, "negative degree for N-pulse");
    affirm(g >= 2*c + 1, "bad degree for N-pulse");
    affirm(g <= udg_pulse_MAX_DEGREE, "degree too high for N-pulse");
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

void udg_pulse_N_mother_main_bezier
  ( udg_cont_t c, 
    udg_degree_t g, 
    udg_pulse_mother_index_t pix,
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
