/* See pst_Kodak_Q13.h */
/* Last edited on 2012-07-22 11:21:05 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_average.h>
#include <r2.h> 
#include <r3.h> 
#include <r3x3.h> 
#include <rn.h> 
#include <gauss_elim.h> 
#include <qmin_simplex.h> 
#include <affirm.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_Kodak_Q13.h>

/* INTERNAL PARAMETERS */

#define pst_Kodak_Q13_scale_strip_height (pst_Kodak_Q13_patch_height - 2*pst_Kodak_Q13_strip_gap)
  /* Y extent of strips in gray scale (mm). */

#define pst_Kodak_Q13_s_strip_bot_y (0.0 + pst_Kodak_Q13_strip_gap)
#define pst_Kodak_Q13_s_strip_top_y (pst_Kodak_Q13_s_strip_bot_y + pst_Kodak_Q13_scale_strip_height)
  /* Bottom and top Y coordinates (mm) of the scale strip. */

#define pst_Kodak_Q13_r0_strip_top_y (0.0 - pst_Kodak_Q13_strip_gap)
#define pst_Kodak_Q13_r0_strip_bot_y (pst_Kodak_Q13_r0_strip_top_y - pst_Kodak_Q13_ref_strip_height)
  /* Bottom and top Y coordinates (mm) of the bottom reference strip (white frame). */

#define pst_Kodak_Q13_r1_strip_bot_y (pst_Kodak_Q13_scale_height + pst_Kodak_Q13_strip_gap)
#define pst_Kodak_Q13_r1_strip_top_y (pst_Kodak_Q13_r1_strip_bot_y + pst_Kodak_Q13_ref_strip_height)
  /* Bottom and top Y coordinates (mm) of the top reference strip (gray background). */

/* INTERNAL PROTOTYPES */

float_image_t *pst_Kodak_Q13_extract_patches
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double step,          /* Output pixel size (mm). */
    double ylo,           /* Bottom Y of patches to extract (mm). */
    double yhi            /* Top Y of patches to extract (mm). */
  );
  /* Extracts a row of {pst_Kodak_Q13_num_steps} patches between Y
    coordinates {ylo} and {yhi} in the chart coordinate system.
    
    The extracted patches have the same integer width (in pixels).
    Each pixel corresponds to a square with side {step} (in mm)
    on the chart. 
    
    All extracted patches (including the first and last ones) are
    equally spaced in the X direction by
    {pst_Kodak_Q13_mid_patch_width} (mm). Their centers are vertically aligned
    with the centers of typical patches of the gray-scale, except for
    the first and last patches.
    
    In principle, the height {dy} of each extracted patch (in chart
    coordinates) is {yhi - ylo}, and the width {dx} is
    {pst_Kodak_Q13_mid_patch_width} minus {pst_Kodak_Q13_patch_trim_x}
    on each side, both rounded down to integer multiples of {step}. The
    procedure fails if {step} is less than {dx} or {dy}. */

/* IMPLEMENTATIONS */

void pst_Kodak_Q13_get_matrix(r2_t *tl, r2_t *tr, r2_t *bl, r2_t *br, r3x3_t *P)
  {
    bool_t debug = TRUE;

    auto void debug_matrix(char *head, r3x3_t *M, double dx, double dy, char *foot);
      /* Prints the matrix {M} and its effect on the corners of a {dx × dy} rectangle. */
    
    /* Find the weights of the points {tl}, {bl}, {br}: */
    r3_t w;
    { r3x3_t Q;
      Q.c[0][0] = -1.0; Q.c[0][1] = -bl->c[0]; Q.c[0][2] = -bl->c[1];
      Q.c[1][0] = +1.0; Q.c[1][1] = +br->c[0]; Q.c[1][2] = +br->c[1];
      Q.c[2][0] = +1.0; Q.c[2][1] = +tl->c[0]; Q.c[2][2] = +tl->c[1];
      r3x3_inv(&Q, &Q);
      w.c[0] = 1.0;  w.c[1] = tr->c[0];  w.c[2] = tr->c[1];
      r3x3_map_row(&w, &Q, &w);
      if (debug)
        { fprintf(stderr, "  weights:");
          fprintf(stderr, "  bl = %11.5f", w.c[0]);
          fprintf(stderr, "  br = %11.5f", w.c[1]);
          fprintf(stderr, "  tl = %11.5f", w.c[2]);
          fprintf(stderr, "\n");
        }
    }
    
    /* Build the matrix {M} that maps the unit square {U} to the four given points: */
    P->c[0][0] = w.c[0]; P->c[0][1] = w.c[0]*bl->c[0]; P->c[0][2] = w.c[0]*bl->c[1];
    P->c[1][0] = w.c[1]; P->c[1][1] = w.c[1]*br->c[0]; P->c[1][2] = w.c[1]*br->c[1];
    P->c[2][0] = w.c[2]; P->c[2][1] = w.c[2]*tl->c[0]; P->c[2][2] = w.c[2]*tl->c[1];
    
    int i, j;
    for (i = 1; i < 3; i++)
      { for (j = 0; j < 3; j++)
          { P->c[i][j] -= P->c[0][j]; }
      }
      
    if (debug) { debug_matrix("unit square:", P, 1.0, 1.0, ""); }

    /* Combine with the matrix that maps the actual Q-13 rectangle to {U}: */
    for (i = 0; i < 3; i++)
      { P->c[1][i] /= pst_Kodak_Q13_total_width;
        P->c[2][i] /= pst_Kodak_Q13_total_height;
      }

    if (debug) { debug_matrix("Q-13 chart:", P, pst_Kodak_Q13_total_width, pst_Kodak_Q13_total_height, ""); }

    void debug_matrix(char *head, r3x3_t *M, double dx, double dy, char *foot)
      { if (head != NULL) { fprintf(stderr, "  %s\n", head); }
        int x, y;
        for (y = 0; y < 2; y++)
          for (x = 0; x < 2; x++)
            { r3_t p = (r3_t){{ 1.0, x*dx, y*dy }};
              r3_t q;
              r3x3_map_row(&p, M, &q);
              r2_t t;
              t.c[0] = q.c[1]/q.c[0]; 
              t.c[1] = q.c[2]/q.c[0];
              rn_gen_print(stderr, 3, p.c, "%+7.2f", "  [ ", " ", " ]");
              rn_gen_print(stderr, 3, q.c, "%+7.2f", " = [ ", " ", " ]");
              rn_gen_print(stderr, 2, t.c, "%+7.2f", " = ( ", " ", " )");
              fprintf(stderr, "\n");
            }
         if (foot != NULL) { fprintf(stderr, "  %s\n", foot); }
      }
  }

float_image_t *pst_Kodak_Q13_extract_chart
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Width of chart image (in pixels). */
  )
  { 
    /* Get input image dimensions (samples): */
    int NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    
    /* Chart dimensions (mm): */
    double mmX = pst_Kodak_Q13_total_width;   /* Chart width in mm. */
    double mmY = pst_Kodak_Q13_total_height;  /* Chart height in mm. */

    /* Output image dimensions (pixels): */
    int ONX = (int)floor(mmX/pixelSize + 0.5);
    int ONY = (int)floor(mmY/pixelSize + 0.5);
    
    /* Left and right edges of chart (mm): */
    double xsz = ONX*pixelSize;
    double xlo = (pst_Kodak_Q13_total_width - xsz)/2;
    double xhi = xlo + xsz;
    
    /* Top and bottom edges of chart (mm): */
    double ysz = ONY*pixelSize;
    double ylo = (pst_Kodak_Q13_total_height - ysz)/2;
    double yhi = ylo + ysz;
    
    /* Create and fill the chart image: */
    float_image_t *omg = float_image_new(NC, ONX, ONY);

    ix_reduction_t red = ix_reduction_SINGLE;
    float undef = 0.5;
    bool_t avg = TRUE;
    int order = 1;
    float_image_transform_copy_persp_rectangle(img, red, xlo, xhi, ylo, yhi, P, undef, avg, order, 0, 0, ONX, ONY, NULL, omg);

    return omg;
  }

float_image_t *pst_Kodak_Q13_extract_gray_scale
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  )
  { 
    double ylo = pst_Kodak_Q13_patch_trim_y;
    double yhi = pst_Kodak_Q13_scale_height - pst_Kodak_Q13_patch_trim_y;
    return pst_Kodak_Q13_extract_patches(img, P, pixelSize, ylo, yhi);
  }
  
float_image_t *pst_Kodak_Q13_extract_body_strip
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  )
  { 
    double ylo = pst_Kodak_Q13_scale_height + pst_Kodak_Q13_ref_strip_gap;
    double yhi = ylo + pst_Kodak_Q13_ref_strip_height;
    return pst_Kodak_Q13_extract_patches(img, P, pixelSize, ylo, yhi);
  }
  
float_image_t *pst_Kodak_Q13_extract_frame_strip
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  )
  { 
    double yhi = 0 - pst_Kodak_Q13_ref_strip_gap;
    double ylo = yhi - pst_Kodak_Q13_ref_strip_height;
    return pst_Kodak_Q13_extract_patches(img, P, pixelSize, ylo, yhi);
  }

float_image_t *pst_Kodak_Q13_extract_patches
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize,          /* Output pixel size (mm). */
    double ylo,           /* Bottom Y of patches to extract (mm). */
    double yhi            /* Top Y of patches to extract (mm). */
  )
  {
    /* Get input image dimensions (samples): */
    int NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    /* Compute extracted patch dimensions {mmX_pt,mmY_pt} in mm: */
    double mmX_pt = pst_Kodak_Q13_mid_patch_width - 2*pst_Kodak_Q13_patch_trim_x;
    double mmY_pt = yhi - ylo;
    /* Compute extracted patch dimensions {SNX_pt,SNY_pt} in pixels: */
    int SNX_pt = (int)floor(mmX_pt/pixelSize + 0.00001);
    int SNY_pt = (int)floor(mmY_pt/pixelSize + 0.00001);
    demand(SNX_pt > 0, "pixel size {pixelSize} is larger than patch width");
    demand(SNY_pt > 0, "pixel size {pixelSize} is larger than patch height");
    /* Recompute ideal extracted patch Y dimensions{mmX_pt,mmY_pt} in mm: */
    mmX_pt = SNX_pt*pixelSize;
    mmY_pt = SNY_pt*pixelSize;
    /* Recompute ideal extracted patch Y range in mm: */
    double yCenter = (ylo + yhi)/2;
    ylo = yCenter - mmY_pt/2;
    yhi = yCenter + mmY_pt/2;
    /* Compute output image dimensions (pixels): */
    int nSteps = pst_Kodak_Q13_num_steps; 
    int SNX = SNX_pt * nSteps;
    int SNY = SNY_pt;
    /* Create output image: */
    float_image_t *omg = float_image_new(NC, SNX, SNY);
    /* Compute X displacement from chart edge to regularized gray scale: */
    double xShift = (pst_Kodak_Q13_total_width - nSteps*pst_Kodak_Q13_mid_patch_width)/2;
    int i;
    for (i = 0; i < nSteps; i++)
      { /* Compute center X coordinate of patch {i}: */
        double xCenter = xShift + (i + 0.5)*pst_Kodak_Q13_mid_patch_width;
        /* Compute left and right edges {xlo,xhi} of patch {i}: */
        double xlo = xCenter - mmX_pt/2;
        double xhi = xCenter + mmX_pt/2;
        /* Extract patch: */
        ix_reduction_t red = ix_reduction_SINGLE;
        float undef = 0.5;
        bool_t avg = TRUE;
        int order = 1;
        float_image_transform_copy_persp_rectangle
          ( img, red, xlo, xhi, ylo, yhi, P, undef, avg, order, i*SNX_pt, 0, SNX_pt, SNY_pt, NULL, omg );
      }
        
    return omg;
  }

double pst_Kodak_Q13_patch_albedo(int i)
  { 
    return pow(0.1, i*0.10 + 0.05);
  }
  
