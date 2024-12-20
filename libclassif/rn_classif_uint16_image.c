/* See rn_classif.h. */
/* Last edited on 2024-12-05 10:24:36 by stolfi */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <filefmt.h>
#include <r2.h>
#include <uint16_image.h>
#include <jsrandom.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <rn_classif.h>
#include <rn_classif_float_image.h>
#include <rn_classif_uint16_image.h>

uint16_image_t *rn_classif_uint16_image
  ( int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1,
    int NC2, 
    rn_classif_labeler_t *lab2,
    int NCD, 
    rn_classif_dataset_t *D, 
    double HD[],
    int classD[], 
    r2_t *ictr, 
    double HV, 
    int NXY, 
    int NSUB, 
    double sigma
  )
  { 
    demand(NA == 2, "problem is not two-dimensional");
    
    float_image_t *fim = float_image_new(3, NXY, NXY);
    
    if ((lab1 != NULL) || (lab2 != NULL))
      { /* Paint the two labelers: */
        frgb_t *cmap = rn_classif_pick_class_colors(NC1);
        rn_classif_float_image_paint_labeler(fim, NA, NC1, lab1, NC2, lab2, ictr, HV, NSUB, sigma, cmap);
        free(cmap);
      }
    else
      { /* paint the whole image white: */
        float_image_fill(fim, 1.0);
      }
      
    /* Paint the dots over it: */
    if (D != NULL)
      { int NSUBD = 3; /* Always antialias the dots. */
        /* Choose site colors: */
        frgb_t *cmap = NULL;
        if (classD != NULL)
          { /* Use class-colored dots: */
            cmap = rn_classif_pick_class_colors(NCD);
            /* Class 0 is black,not white, if background is white/gray: */
            if (lab1 == NULL) { cmap[0] =  (frgb_t){{ 0,0,0 }}; }
          }
        /* Choose style for handicap disks: */
        bool_t fillH = ((HD != NULL) && (lab1 == NULL));
        double penH = 0.75;
        /* Choose style for dots: */
        bool_t fillP = TRUE;
        double radP = 0.0;
        double penP = 0.0;
        if ((classD != NULL) && (cmap != NULL))
          { /* Use class-colored dots: */
            if (lab1 != NULL)
              { /* Colored background, needs black border: */
                radP = 2.00; 
                penP = 0.75;
              }
            else
              { /* White/gray background; no border needed: */
                radP = 1.75; 
                penP = 0.0;
              }
          }
        else
          { /* Classes not provided, use smaller black dots: */
            radP = 1.25;
            penP = 0.00;
          }
        rn_classif_float_image_paint_dataset(fim, NA, NCD, D, HD, classD, ictr, HV, fillH, penH, cmap, fillP, radP, penP, NSUBD);
        free(cmap);
      }
      
    int c;
    for (c = 0; c < 3; c++) { float_image_apply_gamma(fim, c, 0.450, 0.0327); }
    bool_t isMask = FALSE; /* Assume smooth distr of pixel values. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, 3, NULL, NULL, NULL, 255, TRUE, FALSE);
    float_image_free(fim);
    return pim;
  }
