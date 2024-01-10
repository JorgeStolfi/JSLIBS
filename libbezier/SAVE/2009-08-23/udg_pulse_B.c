/* See udg_pulse_B.h */
/* Last edited on 2009-02-09 18:24:20 by stolfi */

#include <math.h>

#include <bz_basic.h>
#include <affirm.h>

#include <udg_pulse.h>
#include <udg_pulse_B.h>

unsigned int udg_pulse_B_num_mothers (udg_cont_t c, udg_degree_t g)
  { affirm(g >= 0, "negative degree");
    affirm(c >= -1, "bad continuity class");
    if (g != c + 1)
      { return 0; }
    else
      { return 1; }
  }

udg_grid_size_t udg_pulse_B_mother_supp_count
  ( udg_cont_t c, 
    udg_degree_t g,
    udg_pulse_mother_index_t pix
  )
  { affirm(c >= -1, "bad contin for B-pulse");
    affirm(g == c + 1, "bad degree for B-pulse");
    affirm(pix == 0, "invalid B-pulse index");
    return g + 1;
  }

void udg_pulse_B_mother_to_bezier
  ( udg_cont_t c, 
    udg_degree_t g,
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t x,
    double *bz
  )
  { 
    affirm(c >= -1, "bad contin for B-pulse");
    affirm(g == c + 1, "bad degree for B-pulse");
    affirm(pix == 0, "invalid B-pulse index");
    if ((x < 0) || (x > g)) 
      { int r;
        for (r = 0; r < g; r++) { bz[r] = 0; }
        return;
      }
    /* Should give a general formula: */
    if (c == -1)
      { /* Only one interval: */
        bz[0] = +1;
      }
    else if (c == 0)
      { if (x == 0)
          { bz[0] = 00; bz[1] = +1; }
        else
          { bz[0] = +1; bz[1] = 00; }
      }
    else if (c == 1)
      { if (x == 0)
          { bz[0] = 00; bz[1] = 00; bz[2] = +1; }
        else if (x == 1)
          { bz[0] = +1; bz[1] = +2; bz[2] = +1; }
        else 
          { bz[0] = +1; bz[1] = 00; bz[2] = 00; }
      }
    else if (c == 2)
      { if (x == 0)
          { bz[0] = 00; bz[1] = 00; bz[2] = 00; bz[3] = +1; }
        else if (x == 1)
          { bz[0] = +1; bz[1] = +2; bz[2] = +4; bz[3] = +4; }
        else if (x == 2)
          { bz[0] = +4; bz[1] = +4; bz[2] = +2; bz[3] = +1; }
        else 
          { bz[0] = +1; bz[1] = 00; bz[2] = 00; bz[3] = 00; }
      }
    else 
      { demand(FALSE, "cont order not implemented"); }
  }
