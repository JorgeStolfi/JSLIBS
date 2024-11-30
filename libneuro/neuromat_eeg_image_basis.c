/* See {neuromat_eeg_image_basis.h}. */
/* Last edited on 2024-11-11 07:38:33 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <rmxn.h>
#include <lsq.h>
#include <float_image.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_image.h>
#include <neuromat_eeg_image_basis.h>
                
float_image_t **neuromat_eeg_image_basis_make
  ( int32_t ne, 
    neuromat_eeg_func_basis_eval_t eval, 
    float_image_t *msk, 
    int32_t msub,
    r3_t *srad, 
    r2_t *ictr, 
    r2_t *irad
  )
  {
    bool_t verbose = TRUE;
    
    demand(msk->sz[0] == 1, "mask should be monochromatic");
    int32_t NX = (int32_t)msk->sz[1];
    int32_t NY = (int32_t)msk->sz[2];
    
    /* Allocate the basis images: */
    float_image_t **bas = notnull(malloc(ne*sizeof(float_image_t*)), "no mem");
    for (uint32_t ie = 0;  ie < ne; ie++) { bas[ie] = float_image_new (1, NX, NY); }
    
    /* Now loop on pixels, paint every image: */
    double *bval = rn_alloc(ne); /* Element values at a certain subsampling point. */
   double *bsum = rn_alloc(ne); /* Summed/averaged values in a pixel. */
    for (uint32_t iy = 0;  iy < NY; iy++)
      { for (uint32_t ix = 0;  ix < NX; ix++)
          { bool_t debugpx = (ix == NX/2) & (iy == NY/2);
            /* Get mask weight of pixel: */
            double mxy = float_image_get_sample(msk, 0, ix, iy);
            demand((mxy >= 0) && (mxy <= 1.0), "invalid mask value");
            /* Set {bval[ie]} to value of element {ie} at center of pixel {ix,iy}: */
            for (uint32_t ie = 0;  ie < ne; ie++) { bval[ie] = 0; }
            if (mxy > 0)
              { /* Sample {msub} by {msub} points {p} in pixel {ix,iy}:  */
                for (uint32_t ie = 0;  ie < ne; ie++) { bsum[ie] = 0; }
                for (uint32_t dy = 0;  dy < msub; dy++) 
                  { for (uint32_t dx = 0;  dx < msub; dx++) 
                      { /* Compute the subsampling point {p3D} on the unit sphere: */
                        r2_t qxy = (r2_t) {{ ix + (dx + 0.5)/msub,  iy + (dy + 0.5)/msub }}; /* Pt in image domain. */
                        r2_t pxy = neuromat_eeg_geom_disk_from_ellipse(&qxy, ictr, irad); /* Corresp point in unit-disk schematic. */
                        r3_t p3D = neuromat_eeg_geom_3D_from_2D(&pxy, NULL, srad); /* Corresp point on unit sphere. */
              
                        if (verbose && debugpx){ r3_gen_print(stderr, &p3D, "%+8.5f", "  p3D = ( ", " ", " )\n"); }
                        eval(ne, bval, &p3D);

                        /* Accumulate: */
                        for (uint32_t ie = 0;  ie < ne; ie++) { bsum[ie] += bval[ie]; }
                      }
                  }
                /* Compute pixel average:: */
                for (uint32_t ie = 0;  ie < ne; ie++) { bval[ie] = bsum[ie]/(msub*msub); }
              }
            /* Now set pixels to {bval[0..ne-1]}: */ 
            for (uint32_t ie = 0;  ie < ne; ie++) 
              { float fval = (float)(bval[ie]);
                if (verbose && debugpx){ fprintf(stderr, "    bval[%3d] = %+8.5f fval = %+8.5f\n", ie, bval[ie], fval); }
                float_image_set_sample(bas[ie], 0, ix, iy, fval);
              }
          }
        if (verbose) { fprintf(stderr, "."); if ((iy == NY-1) || (((iy + 1) % 50) == 0)) { fprintf(stderr, "\n"); } }
      }
      
    free(bval);
    return bas;
  }
