/* See pst_fit_sphere.h */
/* Last edited on 2025-03-01 19:28:11 by stolfi */ 

#define _GNU_SOURCE
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <argparser.h>
#include <r2.h>
#include <rn.h> 
#include <affirm.h> 

#include <pst_fit_sphere_and_lights.h>
#include <pst_fit_light.h>
#include <pst_normal_map.h>
#include <pst_geom.h>
#include <pst_light.h>
#include <pst_shading.h>
#include <pst_argparser.h>
#include <pst_basic.h>

/* INTERNAL PROTOTYPES */

void pst_fit_sphere_adjust_bounds(double *mid, double *adj, double low);
  /* Adjusts the arguments {*mid} and {*adj} so as to clip the range
    {[mid-adj _ mid+adj]} against the range {[low _ +INF]}. 
    In particular, if {mid+adj} is less than {low}, sets {mid} to {low}
    and {adj} to 0. */ 

#define Pr fprintf
#define Er stderr

/* IMPLEMENTATIONS */

double pst_fit_sphere_and_lights
  ( pst_map_vec_t *IMGV,      /* Photos of a spherical object. */  
    ellipse_crs_t *E,       /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,          /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,          /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,          /* Maximum ± adjustment allowed in {str} coordinates. */
    double adjStep,         /* Ideal adjustment step for geometry parameters. */
    pst_light_vec_t *lhtv,  /* (IN/OUT) light model. */
    int adjustDir,          /* Index of lamp for direction adjustment, or -1. */
    double dirStep,         /* Max lamp direction change per iteration (radians). */ 
    double weightBias,      /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative,     /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,      /* Minimum Z of normal to consider in light fitting. */
    int iterLight,          /* Max geometry-fitting iterations. */
    double tolLight,        /* Convergence criterion for geometry fitting. */
    float_image_t *NRM,     /* (OUT) Normal map of sphere. */ 
    pst_map_vec_t *SYNV       /* (OUT) Synthetic photos of the fitted sphere. */
  )
  { bool_t debug = TRUE;
  
    /* Get and check all dimensions: */
    int NF = IMGV->ne;
    demand(lhtv->ne == NF, "bad output light field vector");
    demand(SYNV->ne == NF, "bad output photo vector");

    int NC, NX, NY;
    float_image_get_size(NRM, &NC, &NX, &NY);
    demand(NC == 4, "bad normal map");
    
    int i;
    for (i = 0; i < NF; i++)
      { float_image_check_size(IMGV->e[i], 1, NX, NY, "mismatched input images");
        float_image_check_size(SYNV->e[i], 1, NX, NY, "mismatched output sinthetic images");
      }
      
    /* Grab parameters in local variables {centerX,centerY,radius,stretchX,stretchY}: */
    double centerX = E->ctr.c[0];
    double centerY = E->ctr.c[1];
    double radius = E->rad;
    double stretchX = E->str.c[0];
    double stretchY = E->str.c[1];
      
    /* Adjust {radius}, {radAdj} so as to avoid negative radii: */
    double radMin = 0.2; /* Minimum sphere radius in pixels. */
    pst_fit_sphere_adjust_bounds(&radius, &radAdj, radMin);

    /* We use a slow and dumb method, namely vary each variable parameter by a 
      fixed step, and record the best combination. */
      
    /* Compute number of values to try on each side of central value of each parameter. */

    auto int compute_steps(double adj, double step);
      /* Computes the number of {K} of trial values for a geometry
        parameter which has to be adjusted between {val-adj} and
        {val+adj}, for some center valeu {val}, with granularity at most
        {step}. The count is always odd and positive; if {adj} is zero,
        the result is 1, else it is at least 3. */  

    int compute_steps(double adj, double step)
      { int K = 1 + 2*(step == 0 ? 0 : (int)ceil(adj/step));
        assert((K > 0) && (K % 2 == 1)); 
        return K;
      }

    int KC = compute_steps(ctrAdj,adjStep); /* Trials in {ctr} X and Y. */
    int KR = compute_steps(radAdj,adjStep); /* Trials in {rad}. */
    int KS = compute_steps(strAdj,adjStep); /* Trials in {str} length. */

    auto double try_val(double mid, double adj, int K, int i); 
      /* Computes the trial parameter value number {i} (counting from
        0) given the parameter's middle value {mid} and the amount of
        adjustment {adj} on each side of {mid}. Trial value 0 is
        {mid}, then alternates between {mid+k*step} and
        {mid-k*step}. */
      
    auto double try_val_log(double mid, double adj, int K, int i); 
      /* Similar to {tri_val} but uses logscale instead of linear scale. */
      
    auto double try_val_rel(int K, int i); 
      /* Returns the relative displacement to be used by {try_val}
        (or {try_val_log}, in log scale), ranging in {[-1 _ +1]} . */
    
    double try_val(double mid, double adj, int K, int i)
      { if (adj == 0) { return mid; }
        double r = try_val_rel(K, i);
        return mid + r*adj;
      }
    
    double try_val_log(double mid, double adj, int K, int i)
      { if ((mid == 0) || (adj == 0)) { return mid; }
        double r = try_val_rel(K, i);
        return mid * pow((mid+adj)/mid, r);
      }
    
    double try_val_rel(int K, int i)
      { int H = (K - 1)/2;     /* Number of steps on each side. */
        if (H <= 0) { return 0; }
        int m = (i + 1)/2;     /* Layer index, from 0. */
        int s = 2*(i % 2) - 1; /* Direction, +1 or -1. */
        return m*s/H;
      }
    
    /* Now vary the parameters, save best in {*ctr}, {*rad}, {*str}: */
    double d2Best = +INF;
    int icx, icy, ir, isx, isy;
    ellipse_crs_t try;
    for (icy = 0; icy < KC; icy++)
      { try.ctr.c[1] = try_val(centerY, ctrAdj, KC, icy);
        for (icx = 0; icx < KC; icx++)
          { try.ctr.c[0] = try_val(centerX, ctrAdj, KC, icx);
            for (ir = 0; ir < KR; ir++)
              { try.rad = try_val_log(radius, radAdj, KR, ir);
                for (isy = 0; isy < KS; isy++)
                  { try.str.c[1] = try_val(stretchY, strAdj, KS, isy);
                    for (isx = 0; isx < KS; isx++)
                      { try.str.c[0] = try_val(stretchX, strAdj, KS, isx);
                        if (debug)
                          { double *ctr = try.ctr.c; double *str = try.str.c;
                            Pr(Er, "  tentative geometry:\n");
                            Pr(Er, "    radius  = %8.5f\n", try.rad);
                            Pr(Er, "    center  = ( %8.5f %8.5f )\n", ctr[0], ctr[1]);
                            Pr(Er, "    stretch = ( %+8.5f %+8.5f )\n", str[0], str[1]);
                          }

                        double d2 = pst_fit_sphere_evaluate_geometry
                          ( IMGV, 
                            &try, 
                            lhtv, 
                            adjustDir, dirStep, weightBias, nonNegative, 
                            minNormalZ, iterLight, tolLight, 
                            NRM, SYNV
                          );
                        if (debug) { Pr(Er, "  avg diff sqr = %8.5f\n", d2); }
                        if (debug) { Pr(Er, "\n"); }
                        if (d2 < d2Best) { d2Best = d2; (*E) = try; }
                      }
                  }
              }
          }
      }
        
    /* Recompute best {NRM}, {lmpv}, and {SYNV}: */
    d2Best = pst_fit_sphere_evaluate_geometry
      ( IMGV, 
        E,
        lhtv, 
        adjustDir, dirStep, weightBias, nonNegative, 
        minNormalZ, iterLight, tolLight, 
        NRM, SYNV
      );
    return d2Best;
  }
  
double pst_fit_sphere_evaluate_geometry
  ( pst_map_vec_t *IMGV,      /* Actual photos of a spherical object. */
    ellipse_crs_t *E,       /* Geometric parameters of sphere's projection. */
    pst_light_vec_t *lhtv,  /* (IN/OUT) light model. */
    int adjustDir,          /* Index of lamp for direction adjustment, or -1. */
    double dirStep,         /* Max lamp direction change per iteration (radians). */ 
    double weightBias,      /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative,     /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,      /* Minimum Z of normal to consider in light fitting. */
    int iterLight,          /* Max iterations for light field fitting. */
    double tolLight,        /* Stopping tolerance for light field fitting. */
    float_image_t *NRM,     /* (OUT) Normal map of best-fit sphere. */
    pst_map_vec_t *SYNV       /* (OUT) Synthetic photos of the fitted sphere. */
  )
  {
    bool_t debug = FALSE;
    
    if (debug)
      { Pr(Er, "entering %s\n", __FUNCTION__);
        Pr(Er, "E = "); ellipse_crs_print(Er, E, "%8.2f"); Pr(Er, "\n");
      }
    
    int NF = lhtv->ne; /* Number of light fields. */
    demand(SYNV->ne == NF, "light field count mismatch");
    double d2tot = 0.0;
    
    /* Create the normal map of an ideal sphere with geometry {E}: */
    hr3_pmap_t M = pst_geom_sphere_view_map(C, E);
    pst_normal_map_from_proc(pst_geom_sphere_compute_normal, 2, &xym_to_uvm, &uvw_to_xyz, NRM);
    
    int i;
    for (i = 0; i < NF; i++)
      { float_image_t *IMG = IMGV->e[i];
        /* float_image_t *MSK = MSKV->e[i]; */
        /* Adjust simple light field model: */
        pst_light_t *lht = &(lhtv->e[i]); 
        bool_t pwrAdjust = TRUE;
        bool_t ambAdjust = TRUE;
        if (adjustDir >= 0)
          { /* Adjust lamp direction: */
            if (debug) { Pr(Er, "  adjusting direction of lamp %d\n", adjustDir); }
            pst_lamp_t *src = lht->lmpv.e[adjustDir]; 
            pst_fit_light_single_iterative
              ( IMG, NRM, 
                lht, 
                src, dirStep, pwrAdjust, ambAdjust,
                weightBias, nonNegative, minNormalZ, 
                iterLight, tolLight
              );
          }

        /* Finish off with a global intensity optimization: */
        pst_fit_light_multi(IMG, NRM, lht, weightBias, nonNegative, minNormalZ);

        /* Create synthetic photo: */
        float_image_t *SYN = SYNV->e[i];
        float_image_fill(SYN, 0.0);
        pst_shading_add_diffuse(NRM, NULL, lht, SYN);
        /* Compare actual and synthetic photos: */
        double d2 = pst_fit_sphere_diff_sqr(NRM, IMG, SYN);
        if (debug) { Pr(Er, "    photo %d avg diff sqr = %8.5f\n", i, d2); }
        /* Accumulate the difference: */
        d2tot += d2;
      }

    if (debug)
      { Pr(Er, "leaving %s\n", __FUNCTION__);
        Pr(Er, "  total diff sqr = %8.5f\n", d2tot);
      }

    return d2tot / NF;
  }

double pst_fit_sphere_diff_sqr
  ( float_image_t *NRM,
    float_image_t *AIMG, 
    float_image_t *BIMG
  )
  { /* Get/check sizes: */
    int NC, NX, NY;
    float_image_get_size(AIMG, &NC, &NX, &NY);
    float_image_check_size(BIMG, NC, NX, NY, "bad {BIMG} image");
    float_image_check_size(NRM, 4, NX, NY, "bad normal map");
    /* Compute sum of squared differences: */
    double d2sum = 0.0; /* Sum of squared differences between valid pixels. */
    double np = 0;      /* Number of valid pixels. */
    int c, x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            if (r3_L_inf_norm(&nrm) != 0)
              { /* Valid pixel: */
                for (c = 0; c < NC; c++) 
                  { double a = float_image_get_sample(AIMG, c, x, y);
                    double b = float_image_get_sample(BIMG, c, x, y);
                    double d = a - b;
                    d2sum += d*d;
                    np++;
                  }
              }
          }
      }
    /* Return average over valid pixels: */
    return d2sum/np;
  }
    
void pst_fit_sphere_adjust_bounds(double *mid, double *adj, double low)
  { double md = (*mid), ad = (*adj); 
    if (md + ad <= low)
      { md = low; ad = 0.0; }
    else if (md - ad < low) 
      { md = (low + md + ad)/2; ad = md - low; }
    (*mid) = md; (*adj) = ad;
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
**
**   Copyright © 2006 by the State University of Campinas (UNICAMP).
**
** Created on apr/2006 by Jorge Stolfi, IC-UNICAMP.
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
