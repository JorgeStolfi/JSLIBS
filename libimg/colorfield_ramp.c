/* see colorfield_ramp.h 
** Last edited on 2013-10-21 00:12:42 by stolfilocal
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>
 
#include <colorfield_ramp.h>

cfld_ramp_args_t *cfld_ramp_parse_general(argparser_t *pp)
  {
    cfld_ramp_args_t *rfa = (cfld_ramp_args_t *)malloc(sizeof(cfld_ramp_args_t));
    int32_t k;
    for (k = 0; k < 3; k++)
      { int32_t Hk = (int32_t)argparser_get_next_int(pp, -INT32_MAX, INT32_MAX);
        int32_t Vk = (int32_t)argparser_get_next_int(pp, -INT32_MAX, INT32_MAX);
        frgb_t colork = frgb_parse_color(pp);
        rfa->p[k] =  (cfld_int_pair_t){{Hk, Vk}};
        rfa->color[k] = colork;
      }
      
    return rfa;
  }

cfld_ramp_params_t *cfld_ramp_compute_params
  ( cfld_ramp_args_t *rfa, 
    frgb_adjuster_t *adjust,
    int32_t logarithmic
  )
  { cfld_ramp_params_t *rfp = (cfld_ramp_params_t *)malloc(sizeof(cfld_ramp_params_t));
    frgb_t fv[3];
    cfld_int_pair_t *p = &(rfa->p[0]);
    int32_t k;
    /* Note: user-input colors are all RGB, even when {o->gray} is true. */
    for (k = 0; k < 3; k++)
      { fv[k] = adjust(&(rfa->color[k]), p[k].c[0], p[k].c[1]);
        if (logarithmic) { frgb_log_scale(&(fv[k]), frgb_CHANNELS); }
      }
      
    /* Compute coefficients of linear field: */
    { /* The point matrix: */
      double m00 = 1.0, m01 = p[0].c[0], m02 = p[0].c[1];
      double m10 = 1.0, m11 = p[1].c[0], m12 = p[1].c[1];
      double m20 = 1.0, m21 = p[2].c[0], m22 = p[2].c[1];
      
      /* The color matrix: */
      double f00 = fv[0].c[0], f01 = fv[0].c[1], f02 = fv[0].c[2];
      double f10 = fv[1].c[0], f11 = fv[1].c[1], f12 = fv[1].c[2];
      double f20 = fv[2].c[0], f21 = fv[2].c[1], f22 = fv[2].c[2];
      
      /* Pivoting {m00,m10,m20}: */
      m10 -= m00; m11 -= m01; m12 -= m02; f10 -= f00; f11 -= f01; f12 -= f02;
      m20 -= m00; m21 -= m01; m22 -= m02; f20 -= f00; f21 -= f01; f22 -= f02;
      
      /* Pivoting {m11,m21}: */
      { double d = sqrt(m11*m11 + m21*m21);
        double c = m11/d, s = -m21/d;
        double t;
        t = c*m21 + s*m11; m11 = c*m11 - s*m21; m21 = 0.0 /*t*/; 
        t = c*m22 + s*m12; m12 = c*m12 - s*m22; m22 = t; 
        
        t = c*f20 + s*f10; f10 = c*f10 - s*f20; f20 = t; 
        t = c*f21 + s*f11; f11 = c*f11 - s*f21; f21 = t; 
        t = c*f22 + s*f12; f12 = c*f12 - s*f22; f22 = t; 
      }
      
      /* Check for degenerate triangle: */
      assert(fabs(m22) >= 1.0e-8); 
      
      /* Normalize {m22}: */
      { double t = m22;
        m22 = 1.0; /* {m22/t} */
        f20 /= t; f21 /= t; f22 /= t;
      }

      /* Clear rest of column 2: */
      { double t = m12;
        m12 = 0.0; /* {m12-t*m22} */
        f10 -= t*f20;
        f11 -= t*f21;
        f12 -= t*f22;
      }
      
      { double t = m02;
        m02 = 0.0; /* {m02-t*m22} */
        f00 -= t*f20;
        f01 -= t*f21;
        f02 -= t*f22;
      }

      /* Check for degenerate triangle: */
      assert(fabs(m11) >= 1.0e-8);
      
      /* Normalize {m11} */
      { double t = m11;
        m11 = 1.0; /* {m11/t} */
        f10 /= t; f11 /= t; f12 /= t;
      }
      
      /* Clear rest of column 1: */
      { double t = m01;
        m01 = 0.0; /* {m01-t*m11} */
        f00 -= t*f10;
        f01 -= t*f11;
        f02 -= t*f12;
      }
      
      /* m00 is already 1.0 */ 
      assert(m00 == 1.0);
      
      rfp->forg = (frgb_t){{ (float)f00, (float)f01, (float)f02 }};
      rfp->fcol = (frgb_t){{ (float)f10, (float)f11, (float)f12 }};
      rfp->frow = (frgb_t){{ (float)f20, (float)f21, (float)f22 }};
    }

    return rfp;
  }

void cfld_ramp_eval
  ( cfld_ramp_params_t *rfp,
    int32_t logarithmic, 
    int32_t col, 
    int32_t row,
    frgb_t *fv,
    int32_t chns
  )
  {
    float *forg = rfp->forg.c;
    float *fcol = rfp->fcol.c;
    float *frow = rfp->frow.c;
    int32_t i;
    for (i = 0; i < chns; i++)
      { fv->c[i] = (float)(forg[i] + fcol[i]*(double)col + frow[i]*(double)row); 
        if (logarithmic) { fv->c[i] = (float)exp(fv->c[i]); }
      }
  }
