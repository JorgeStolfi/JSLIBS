/* See {ball.h} */
/* Last edited on 2024-11-23 06:22:24 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <ball.h>

double ball_vol(uint32_t d)
  {
    if (d == 0)
      { return 1; }
    else if (d == 1)
      { return 2; }
    else
      { return M_PI/(((double)d)/2) * ball_vol(d-2); }
  }

double ball_cap_vol_frac_pos(uint32_t d, double z)
  {
    /* fprintf(stderr, "    vol_frac_pos(%3d,%.10f)", d, z); */
    double f;
    if (z < -1.0) 
      { f = 0.0; }
    else if (z > 1.0) 
      { f = 1.0; }
    else if (d == 0)
      { /* The zero-dimensional sphere has just two points: */
        f = (z == -1.0 ? 0.25 : (z == 1 ? 0.75 : 0.50));
      }
    else if (z == -1.0) 
      { f = 0.0; }
    else if (z == 1.0) 
      { f = 1.0; }
    else if (d == 1)
      { /* The one-dimensional sphere is the interval [-1 _ +1]: */
        f = (z + 1)/2;
      }
    else
      { double w = asin(z);
        f = 0.5 + ball_zone_vol_frac_ang(d,w);
      }
    /* fprintf(stderr, " = %.10f\n", f); */
    return f;
  }

double ball_zone_vol_frac_ang(uint32_t d, double w)
  {
    /* fprintf(stderr, "    vol_frac_ang(%3d,%.10f)", d, w); */
    if (w + M_PI/2 < +1.0e-14)
      { return -0.5; }
    else if (w - M_PI/2 > -1.0e-14)
      { return +0.5; }
    else if (d == 0)
      { return 0; }
    else if (d == 1)
      { return 0.5*sin(w); }
    else if (d == 2)
      { return (w + 0.5*sin(2*w))/M_PI; }
    else
      { int32_t p = (d % 2);
        double cw = cos(w);

        double K = (p == 0 ? cw : 1);
        double f = K;
        int32_t i;
        for (i = 3 - p; i < d; i += 2)
          { K = K * cw*cw*((double)i-1)/((double) i); 
            f += K; 
          }
        f = f*sin(w);
        if (p == 0) { f = 2*(f + w)/M_PI; }
        f = 0.5*f;
        /* fprintf(stderr, " = %.10f\n", f); */
        return f;
      }
  }
