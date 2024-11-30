/* Last edited on 2023-11-26 07:09:50 by stolfi */ 
/* Test of the {ppv_array_blur.h} module. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <uint16_image.h>
#include <uint16_image_write_png.h>

#include <ppv_array.h>
#include <ppv_image.h>

#include <ppv_array_blur.h>

void pbt_do_test(ppv_dim_t d, int32_t radius, int32_t stride);
ppv_array_t *pbt_make_array(ppv_dim_t d);

int32_t main (int32_t argn, char **argv);

int32_t main (int32_t argn, char **argv)
  {
    ppv_dim_t dmax = 4;
    int32_t rmax  = 3;

    for (ppv_dim_t d = 0; d <= dmax; d++)
      { for (uint32_t radius = 1;  radius < rmax; radius++)
          { for (uint32_t i = 0;  i < 3; i++)
              { int32_t stride;
                if (i == 0)
                  { stride = 1; }
                else if (i == 1)
                  { stride = radius+1; }
                else 
                  { /* Look for an allowed value of {stride} other than 0 or {radius}: */
                    stride = -1;
                    for (int32_t s = 2; s <= radius; s++)
                      { if ((radius+1) % s == 0) { stride = s; } }
                  }
                if (stride >= 1) { pbt_do_test(d, radius, stride); } 
              }
          }
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
     
void pbt_do_test(ppv_dim_t d, int32_t radius, int32_t stride)
  { fprintf(stderr, "--- testing {ppv_array_blur}");
    fprintf(stderr, " dimension %d radius = %d stride = %d ---\n", d, radius, stride);
    
    /* Create the weight table: */
    int32_t szw = 2*radius + 1;
    double wt[szw]; /* Hahn weight table, sums to 1. */
    int32_t stride_nat;
    wt_table_hann_fill(szw, 0.0, wt, &stride_nat);
    demand(stride_nat != 0, "weight table is not partition of constant");
    demand(stride_nat % stride == 0, "given {stride} is wrong");
    wt_table_normalize_sum(szw, wt);
    wt_table_print(stderr, "hann", szw, wt, stride);
    
    /* Create the array and fill it with a suitable test pattern: */
    ppv_array_t *A = pbt_make_array(d);
    
    /* Choose the number of bits per sample of {G}: */
    int32_t nw = (int32_t)ipow(szw, d); /* Samples in neighborhood. */
    int32_t lognw = 0;
    while ((1 << lognw) < nw) { lognw++; }
    int32_t bps_min = (int32_t)imax(1, lognw*A->bps);
    int32_t bps_max = (int32_t)imin(ppv_MAX_BPS, 2*bps_min + 1); 
    assert(bps_min <= bps_max);
    ppv_nbits_t bpsG = (ppv_nbits_t)int32_abrandom(bps_min, bps_max);
    fprintf(stderr, "nw = %d bpsG = %d\n", nw, bpsG);
    ppv_sample_t maxsmpG = ppv_max_sample(bpsG);

    /* Blur array: */
    /* ??? Should test with custom {quantize,floatize}. ??? */
    ppv_array_t *G = ppv_array_blur(A, NULL, maxsmpG, radius, wt, stride, NULL);

    /* Write out as image if possible: */
    char *fpref = jsprintf("out/J_d%d_r%02d_s%02d", d, radius, stride); 
    if ((G->d == 2) && (bpsG <= 16))
      { /* Write {G} as a greyscale image: */
        uint16_image_t *J = ppv_image_from_array(G);
        char *fname = txtcat(fpref, ".png");
        uint16_image_write_png_named(fname, J, 1.0, TRUE);
        free(fname);
      }
    else if ((G->d == 3) && (bpsG <= 16))
      { /* Write {G} as a stack of greyscale images: */
        int32_t nc = (int32_t)G->size[2]; /* Number of slices. */
        for (uint32_t ic = 0;  ic < nc; ic++)
          { ppv_array_t *F = ppv_slice(G, 2, ic);
            uint16_image_t *J = ppv_image_from_array(F);
            char *fname = jsprintf("%s_%05d.png", fpref, ic);
            uint16_image_write_png_named(fname, J, 1.0, TRUE);
            free(fname);
          }
      }
    free(fpref);
    
    return;
  }
        
ppv_array_t *pbt_make_array(ppv_dim_t d)
  { 
    ppv_size_t sz[d];
    ppv_choose_test_size(d, 1000000, sz);
    ppv_sample_t maxsmp = (ppv_sample_t)abrandom(1, 8);
    ppv_array_t *A = ppv_array_new(d, sz, maxsmp);
    ppv_throw_balls(A);
    return A;
  }

