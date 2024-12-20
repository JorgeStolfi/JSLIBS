/* See {frgb_inksep.h}. */
/* Last edited on 2024-12-05 01:02:32 by stolfi */ 

#include <stdint.h>
#include <math.h>

#include <affirm.h>
#include <r4x4.h>
#include <r4.h>
#include <bool.h>
#include <frgb.h>

#include <frgb_inksep.h>

#define CHNS frgb_inksep_CHNS
#define LAYS frgb_inksep_LAYS

void frgb_inksep_shrink_color_set(r4x4_t *mix, double shrink)
  { if (shrink != 0)
      { /* Shrink the colors towards their average: */
        double s = shrink, t = 1 - s;
        for(int32_t c = 0; c < CHNS; c ++)
          { double bj = 0;
            for (int32_t i = 0;  i < LAYS; i++) { bj += mix->c[i][c]; }
            bj /= LAYS;
            for (int32_t i = 0;  i < LAYS; i++) { mix->c[i][c] = t*mix->c[i][c] + s*bj; }
          }
      }
  }

void frgb_inksep_separate_layers(frgb_t fv, r4x4_t *sep, r4_t *m)
  { /* Convert color {fv[0..chns-1]} into RGBW homogeneous quadruple: */
    r4_t y;
    for (int32_t c = 0;  c < LAYS; c++) { y.c[c] = fv.c[c]; }
    y.c[LAYS-1] = 1.0;
    
    /* Perform convex color decomposition: */
    r4x4_map_row(&y, sep, m);
  }
  
void frgb_inksep_clip_to_simplex(r4_t *m, int32_t *nbadP)
  { /* Compute {z} = {x} clipped to simplex: */
    int32_t nbad = 0; /* Number of negative barycoords. */
    int32_t ngud = 0; /* Number of non-negative barycoords. */
    int32_t kord[LAYS]; /* First {nbad} elems are bad, rest are good. */
    for (int32_t i = 0;  i < LAYS; i++) 
      { if (m->c[i] < 0.0) 
          { kord[nbad] = i; nbad++; }
        else
          { kord[LAYS-1-ngud] = i; ngud++; }
      }
      
    if (nbad > 0)
      { /* Must clip the {z} tuple. */
        /* Should clip to the nearest point in RGB distance, but */
        /* It is easier to clip in barycentric distance. */

        /* Compute sum of negative elems, set them to 0: */
        double sbad = 0.0;
        for (int32_t k = 0;  k < nbad; k++)
          { int32_t i = kord[k]; sbad += m->c[i]; m->c[i] = 0.0; }

        /* Distribute {sbad} over other elems: */
        double sgud = 1.0 - sbad;
        for (int32_t k = 0;  k < ngud; k++)
          { int32_t i = kord[LAYS-1-k]; 
            m->c[i] += sbad*m->c[i]/sgud;
          }
      }
    (*nbadP) = nbad;
  }

void frgb_inksep_reveal_colors(r4_t *m)
  { double t = 1.0;
    for (int32_t i = LAYS-1; i > 0; i--)
      { double v = m->c[i];
        if (t < v) { t = v; }
        m->c[i] = (t == 0.0 ? 0.5 : v / t);
        t -= v; if (t < 0.0) { t = 0.0; }
      }
    m->c[0] = 1.0;
  }

void frgb_inksep_compute_separations(frgb_t bgColor, r4_t *m, r4x4_t *mix, frgb_t gv[])
  { for (int32_t i = 0;  i < LAYS; i++)
      { double s = m->c[i], t = 1 - s;
        for (int32_t c = 0;  c < CHNS; c++)
          { gv[i].c[c] = (float)(s*mix->c[i][c] + t*bgColor.c[c]); }
      }
  }

frgb_t frgb_inksep_remap_color(r4_t *m, r4x4_t *syn)
  {
    r4_t sh;
    r4x4_map_row(m, syn, &sh);
    frgb_t sv;
    for (int32_t c = 0;  c < CHNS; c++) { sv.c[c] = (float)(sh.c[c]/sh.c[3]); }
    return sv;
  }
