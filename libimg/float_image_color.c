/* See {float_image_color.h}. */
/* Last edited on 2024-12-04 23:26:53 by stolfi */

#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_color.h>
#include <frgb.h>
#include <frgb_ops.h>
   
frgb_t fic_get_frgb_pixel(float_image_t *A, int32_t cR, int32_t cG, int32_t cB, int32_t x, int32_t y)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    demand((cR >= 0) && (cR < NC), "bad R channel");
    if ((cG < 0) || (cG >= NC)) { cG = cR; }
    if ((cB < 0) || (cB >= NC)) { cB = cR; }
    demand((x >= 0) && (x < NX), "bad x");
    demand((y >= 0) && (y < NY), "bad y");
    frgb_t p;
    float *sp = float_image_get_sample_address(A, 0, x, y);
    p.c[0] = sp[cR*A->st[0]];
    p.c[1] = sp[cG*A->st[0]];
    p.c[2] = sp[cB*A->st[0]];
    return p;
  }

void fic_set_frgb_pixel(float_image_t *A, int32_t cR, int32_t cG, int32_t cB, int32_t x, int32_t y, frgb_t *p)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    demand((cR >= 0) && (cR < NC), "bad R channel");
    demand((cG >= 0) && (cG < NC), "bad G channel");
    demand((cB >= 0) && (cB < NC), "bad B channel");
    demand((x >= 0) && (x < NX), "bad x");
    demand((y >= 0) && (y < NY), "bad y");
    float *sp = float_image_get_sample_address(A, 0, x, y);
    sp[cR*A->st[0]] = p->c[0];
    sp[cG*A->st[0]] = p->c[1];
    sp[cB*A->st[0]] = p->c[2];
  }

void fic_normalize_colors(float_image_t *A, int32_t cR, int32_t cG, int32_t cB)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    if ((NX == 0) || (NY == 0)) { /* Nothing to do: */ return; }
    demand((cR >= 0) && (cR < NC), "bad R channel");
    demand((cG >= 0) && (cG < NC), "bad G channel");
    demand((cB >= 0) && (cB < NC), "bad B channel");
    
    float vmax, vmin;
    int32_t c, x, y;
    vmax = -INF;
    vmin = +INF;
    for (y = 0; y < NY; y++)
      for (x = 0; x < NX; x++)
        { frgb_t p = fic_get_frgb_pixel(A, cR, cG, cB, x, y);
          for (c = 0; c < 3; c++)
            { if (p.c[c] > vmax) vmax = p.c[c];
              if (p.c[c] < vmin) vmin = p.c[c];
            }
        }
    
    double a, b;
    if (vmin < vmax)
      { /* Map {[vmin _ vmax]} to {[0 _ 1]}: */
        a = 1.0/(vmax - vmin); b = 0;
      }
    else
      { /* Set all pixels to 0.5 gray: */ 
        a = 0.0; b = 0.5;
      }

    for (y = 0; y < NY; y++)
      for (x = 0; x < NX; x++)
        { frgb_t p = fic_get_frgb_pixel(A, cR, cG, cB, x, y);
          p.c[0] = (float)(a * (p.c[0] - vmin) + b);
          p.c[1] = (float)(a * (p.c[1] - vmin) + b);
          p.c[2] = (float)(a * (p.c[2] - vmin) + b);
          fic_set_frgb_pixel(A, cR, cG, cB, x, y, &p);
        }
  }

