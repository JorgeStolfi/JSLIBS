/* See psp_pulse_mother_H.h */
/* Last edited on 2009-08-23 18:24:34 by stolfi */

#include <math.h>

#include <bz_basic.h>
#include <affirm.h>

#include <psp_pulse.h>
#include <psp_pulse_mother_H.h>

unsigned int psp_pulse_mother_H_count (psp_cont_t c, psp_degree_t g)
  { affirm(g >= 0, "negative degree");
    affirm(c >= -1, "bad continuity class");
    if ((c < 0) || (g != 2*c + 1))
      { return 0; }
    else
      { return c+1; }
  }

psp_grid_size_t psp_pulse_mother_H_supp_count
  ( psp_cont_t c, 
    psp_degree_t g,
    psp_pulse_mother_index_t pix
  )
  { affirm(c >= 0, "bad contin for H-pulse");
    affirm(g == 2*c + 1, "bad degree for H-pulse");
    affirm((pix >= 0) && (pix <= c), "invalid H-pulse index");
    return 2;
  }

void psp_pulse_mother_H_to_bezier
  ( psp_cont_t c, 
    psp_degree_t g,
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t x,
    double *bz
  )
  { 
    affirm(c >= 0, "bad contin for H-pulse");
    affirm(g == 2*c + 1, "bad degree for H-pulse");
    affirm((pix >= 0) && (pix <= c), "invalid H-pulse index");
    int r, s;
    int xip = c - pix;
    if (x == 0)
      { /* Coeffs {0..c} are all zero: */
        for (r = 0; r <= c; r++) { bz[r] = 0; }
        /* Coeffs {c+1..g} are obtained by {m}-fold backwards summation of (-1)^m: */
        for (r = g; r >= g-xip; r--) { bz[r] = 1; }
        for(s = xip+1; s <= c; s++)
          { double sum = bz[g];
            bz[g] = 0;
            for (r = g-1; r > g-s; r--) { double t = bz[r]; bz[r] = -sum; sum += t; }
            bz[g-s] = -sum;
          }
      }
    else if (x == 1)
      { /* Coeffs {0..c} are obtained by {m}-fold summation of 1,...1: */
        for (r = 0; r <= xip; r++) { bz[r] = 1; }
        for(s = xip+1; s <= c; s++)
          { double sum = bz[0];
            bz[0] = 0;
            for (r = 1; r < s; r++) { double t = bz[r]; bz[r] = sum; sum += t; }
            bz[s] = sum;
          }
        /* Coeffs {c+1..g} are all zero: */
        for (r = g; r > c; r--) { bz[r] = 0; }
      }
    else
      { /* Coefs are all zero: */
        for (r = 0; r <= g; r++) { bz[r] = 0; }
      }
  }

