/* See pst_virtual_gauge.h */
/* Last edited on 2025-01-04 04:52:35 by stolfi */ 

#include <stdint.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <r3.h> 
#include <affirm.h> 
#include <ellipse_crs.h> 
#include <ellipse_crs_args.h> 

// #include <pst_normal_map.h>
#include <pst_basic.h>
// #include <pst_geom.h>
#include <pst_shading.h>
#include <pst_camera.h>

#include <pst_virtual_gauge.h>

void pst_virtual_gauge_color_in_pixel
  ( int32_t x,
    int32_t y,
    ellipse_crs_t *E, 
    r3x3_t *VM, 
    frgb_t *albedo,
    pst_shading_func_t* shading,
    frgb_t *val_P,
    r3_t *vnr_P,
    double *opac_P
  );
  /* Computes the mean apparent color (per-channel radiance) {val}, in
    the pixel with bottom left corner {(x,y)}, of the synthetic image of
    a virtual gauge described by {E,VM,albedo}, and the {shading}
    function.  Also computes the mean normal direction {vnr} of the 
    gauge's surface in that pixel. Also computes the fraction {opac} of the pixel that is
    covered by the gauge's projection.  These results are returned
    in {*val_P}, {*vnr_P}, and {opac_P}.
    
    The color {val} is the weighted average of the color {valp(p)} at
    several sampoints (sub-sampling points) {p} on the image in and
    around the pixel {x,y}. The opacity is the fraction of those
    sampoints that fall inside the gauge's projection.
    
    The normal {vnr} is the weighted average of the normal {vnrp(p)}
    at those sampoints, normalized to unit length. If the opacity 
    os zero, meaning that no sampoint hit the gauge's projection,
    the mean {vnr} is set to {(0,0,0)}.
    
    The point color {valp(p)} takes into account the surface
    normal direction {vnrp(p)}, the viewing angle, the lighting conditions, and
    the surface albedo. Specifically, it will be the apparent
    intensity (radiance, brightness) value returned by the function
    {shading(anrm)} times {albedo}, channel by channel. 
    
    The gauge's projection on the image is assumed to be the ellipse
    {E}. The procedure computes the true normal {anrm} at the point of
    the gauge or background visible at point {p} of the image. It
    assumes that the horizontal distance from the gauge to the /optical
    axis/ may not be negligible, so that the view direction {view} (from
    the gauge's surface towards the camera's optical center) may be
    significantly tilted relative to the vertical. However, it assumes
    that the dustance from the gauge to the /camera center/ (virtual
    observer's position) is large enough, compared to the radius of the
    gauge, for {view} to be assumed constant over the entire gauge; so
    that an entire hemisphere of the latter is visible. The view matrix
    {VM} will be used to correct the computed surface normal to account
    for the oblique view. */

/* IMPLEMENTATIONS */

void pst_virtual_gauge_paint
  ( pst_virtual_gauge_data_t *gd,
    pst_shading_func_t* shading,
    float_image_t *img,
    float_image_t *nrm
  )
  { int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    demand((NC == 3) || (NC == 4), "image must have 3 or 4 channels");
    if (nrm != NULL) { float_image_check_size(nrm, 3, NX, NY); }
    
    /* Get the normal correction matrix: */
    /* !!! Should compute {p} from the point {p} and camera position {obs} !!! */
    r3x3_t VM = pst_camera_normal_correction_matrix(&(gd->view));
    
    for (int32_t x = 0; x < NX; x++)
      { for (int32_t y = 0; y < NY; y++)
          { frgb_t val;
            double opac;
            r3_t vnr;
            pst_virtual_gauge_color_in_pixel(x, y, &(gd->E), &VM, &(gd->albedo), shading, &val, &vnr, &opac);
            assert(isfinite(opac));
            if (opac > 0)
              { for (int32_t c = 0; c < NC; c++)
                  { float *smp = float_image_get_sample_address(img, c, x, y);
                    double smpc;
                    if (c < 3)
                      { double valc = val.c[c];
                        assert(isfinite(valc));
                        smpc = (opac == 1 ? valc : ((1-opac)*(*smp) + opac*valc));
                      }
                    else 
                      { assert(c == 3);
                        smpc = opac;
                      }
                    (*smp) = (float)smpc;
                  }
              }
            if (nrm != NULL) 
              { for (int32_t c = 0; c < 3; c++)
                  { float_image_set_sample(nrm, c, x, y, (float)(vnr.c[c])); }
              }
          } 
      }
  }

void pst_virtual_gauge_color_in_pixel
  ( int32_t x,
    int32_t y,
    ellipse_crs_t *E, 
    r3x3_t *VM, 
    frgb_t *albedo,
    pst_shading_func_t* shading,
    frgb_t *val_P,
    r3_t *vnr_P,
    double *opac_P
  )
  { /* !!! Should use Hann window weights extending outside the pixel. !!! */
    int32_t HS = 2;  /* Pixel subsamples are numbered {-HS..+HS} along each axis. */
    /* int32_t NS = (2*HS+1)*(2*HS+1); */  /* Total number of pixel sub-samples. */
    int32_t NCC = 3; /* Count of color channels (excluding opacity). */
    double sum_hwv[NCC]; /* Sum of {hit*weight*color} */
    double sum_hwn[3];   /* Sum of {hit*weight*normal}. */
    double sum_hw = 0;   /* Sum of {hit*weight} */
    double sum_w = 0;    /* Sum of {weight}. */
    uint32_t num_h = 0;  /* Count of hits. */
    for (int32_t c = 0; c < NCC; c++) { sum_hwv[c] = 0; }
    for (int32_t kx = -HS; kx <= +HS; kx++)
      { double px = (double)x + 0.5 + (double)(kx)/(HS+1.0);
        for (int32_t ky = -HS; ky <= +HS; ky++)
          { double py = (double)y + 0.5 + 0.5*(double)(ky)/(HS+1.0);
            r2_t p = (r2_t){{ px, py }}; /* Subsampling point. */
            double w = 1.0; /* Sampoint weight -- for now. */
            sum_w += w;
            /* Subtract the center and undo the stretching: */
            r2_t q = ellipse_crs_relative_coords(E, &p);
            
            double r2 = r2_norm_sqr(&q);
            if (r2 <= 1.0)
              { /* Compute the view-relative normal {vnrm} to the sphere: */
                r3_t vnrm = (r3_t){{ q.c[0], q.c[1], sqrt(1 - r2) }};
                /* Compute the absolute normal {anrm}: */
                r3_t anrm; r3x3_map_row(&vnrm, VM, &anrm);
                frgb_t shd = shading(&anrm);
                double pval[NCC];
                bool_t ok = TRUE;
                for (int32_t c = 0; c < NCC; c++) 
                  { pval[c] = shd.c[c]* albedo->c[c];
                    if (! isfinite(pval[c])) { ok = FALSE; }
                  }
                if (ok)
                  { /* Sampoint is valid: */
                    num_h++;
                    sum_hw += w;
                    for (int32_t c = 0; c < NCC; c++) { sum_hwv[c] += w*pval[c]; }
                    for (int32_t c = 0; c < 3; c++) { sum_hwn[c] += w*anrm.c[c]; }
                  }
              }
          }
      }
    assert(sum_w > 0);
    /* double opac = ((double)num_h)/((double)NS); */
    double opac = sum_hw/sum_w;
    frgb_t val;
    r3_t vnr;
    if (sum_hw > 0)
      { for (int32_t c = 0; c < NCC; c++) { val.c[c] = (float)(sum_hwv[c] / sum_hw); } 
        for (int32_t c = 0; c < 3; c++) {vnr.c[c] = sum_hwn[c]/sum_hw; }
        (void)r3_dir(&vnr, &vnr);
      }
    else
      { val = (frgb_t){{ NAN, NAN,NAN }}; 
        vnr = (r3_t){{ 0,0,0 }};
      }
    (*val_P) = val;
    (*vnr_P) = vnr;
    (*opac_P) = opac;
  }

void pst_virtual_gauge_args_parse(argparser_t *pp, pst_virtual_gauge_data_t *gd)
  { 
    ellipse_crs_args_parse(pp, &(gd->E), NULL, NULL, NULL);
    
    if (argparser_keyword_present_next(pp, "view"))
      { gd->view.c[0] = argparser_get_next_double(pp, -100000.0, +100000.0);
        gd->view.c[1] = argparser_get_next_double(pp, -100000.0, +100000.0);
        gd->view.c[2] = argparser_get_next_double(pp, -100000.0, +100000.0);
        /* Normalize to unit length: */
        (void)r3_dir(&(gd->view), &(gd->view));
      }
    else
      { gd->view = (r3_t){{ NAN, NAN, NAN }}; } 
      
    if (argparser_keyword_present_next(pp, "albedo"))
      { gd->albedo.c[0] = (float)argparser_get_next_double(pp, 0.000, 1000.0);
        int c = 1;
        while ((c < 3) && argparser_next_is_number(pp))
          { gd->albedo.c[c] = (float)argparser_get_next_double(pp, 0.000, 1000.0);
            c++;
          }
        /* Replicate last value through all remaining channels: */
        while (c < 3) { gd->albedo.c[c] = gd->albedo.c[c-1]; c++; }
      }
    else
      { gd->albedo = (frgb_t){{ NAN, NAN, NAN }}; }
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
**
**   Copyright © 2004 by the Fluminense Federal University (UFF).
**
** Created on jul/2005 by Rafael Saracchini, IC-UFF.
** Modified by Jorge Stolfi, mar/2006.
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
