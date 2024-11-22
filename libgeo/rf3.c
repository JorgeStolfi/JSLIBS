/* See rf3.h. */
/* Last edited on 2024-11-20 13:55:39 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <rf3.h>

#include <affirm.h>
#include <jsrandom.h>

#define N 3

rf3_t rf3_add (rf3_t* const a, rf3_t* const b)
  { rf3_t r;
    r.c[0] = (float)(((double)a->c[0]) + b->c[0]);
    r.c[1] = (float)(((double)a->c[1]) + b->c[1]);
    r.c[2] = (float)(((double)a->c[2]) + b->c[2]);
    return r;
  }

rf3_t rf3_sub (rf3_t* const a, rf3_t* const b)
  { rf3_t r;
    r.c[0] = (float)(((double)a->c[0]) - b->c[0]);
    r.c[1] = (float)(((double)a->c[1]) - b->c[1]);
    r.c[2] = (float)(((double)a->c[2]) - b->c[2]);
    return r;
  }

rf3_t rf3_scale (double s, rf3_t* const a)
  { rf3_t r;
    r.c[0] = (float)(s * a->c[0]);
    r.c[1] = (float)(s * a->c[1]);
    r.c[2] = (float)(s * a->c[2]);
    return r;
  }

rf3_t rf3_mix (double s, rf3_t* const a, double t, rf3_t* const b)
  { rf3_t r;
    r.c[0] = (float)(s * a->c[0] + t * b->c[0]);
    r.c[1] = (float)(s * a->c[1] + t * b->c[1]);
    r.c[2] = (float)(s * a->c[2] + t * b->c[2]);
    return r;
  }

float rf3_max (rf3_t* const a)
  { return fmaxf(a->c[0], fmaxf(a->c[1], a->c[2])); }

rf3_t rf3_rot_axis (rf3_t* const a, uint32_t i, uint32_t j, double ang)
  { demand((i >= 0) && (i < N), "rf3_rot_axis: bad index {i}");
    demand((j >= 0) && (j < N), "rf3_rot_axis: bad index {j}");
    demand(i != j, "r4_rot_axis: axes not distinct");
    double c = cos(ang);
    double s = sin(ang);
    double x = + c*a->c[i] - s*a->c[j];
    double y = + s*a->c[i] + c*a->c[j];
    rf3_t r = (*a);
    r.c[i] = (float)x;
    r.c[j] = (float)y;
    return r;
  }

rf3_t rf3_rot_gen (rf3_t* const a, rf3_t* const d, double ang)
  { fprintf(stderr, "!! {rf3_rot_gen} not implemented yet!\n");
    return (*a);
  }

double rf3_norm (rf3_t* const a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    return sqrt(a0*a0 + a1*a1 + a2*a2);
  }

rf3_t rf3_throw_cube (void)
  { rf3_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(2.0 * frandom() - 1.0); }
    return r;
  }
