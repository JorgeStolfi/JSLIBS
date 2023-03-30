/* Last edited on 2023-03-18 10:57:41 by stolfi */ 
/* Test of the {ppv_image.h} module. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <uint16_image.h>
#include <uint16_image_write_png.h>

#include <ppv_array.h>
#include <ppv_image.h>

void pit_do_test(ppv_dim_t d, int32_t chns, uint16_t maxsmp);
void pit_fill_array(ppv_array_t *A);

int32_t main (int32_t argn, char **argv)
  {
    ppv_dim_t dmin = 2; /* Min number of axes. */
    ppv_dim_t dmax = 3; /* Max number of axes. */

    for (ppv_dim_t d = dmin; d <= dmax; d++)
      { int32_t chmin = (d == 2 ? 1 : 2);  /* Max count of color channels. */
        int32_t chmax = (d == 2 ? 1 : 4);  /* Max count of color channels. */
        for (int32_t chns = chmin; chns <= chmax; chns++)
          { for (int32_t kb = 0; kb < 4; kb++) 
              { uint16_t maxsmp = (uint16_t)(ipow(5,kb+1)-4);
                 pit_do_test(d, chns, maxsmp);
              }
          }
      }
    fprintf(stderr, "done.\n");
    return 0;
  }

void pit_do_test(ppv_dim_t d, int32_t chns, uint16_t maxsmp)
  { 
    fprintf(stderr, "=== testing with d = %d chns = %d maxsmp = %d ===\n", d, chns, maxsmp);
    
    demand(maxsmp < (1 << 16), "invalid max sample value"); 
    
    /* Chose image dimensions and {maxsmp}: */
    int32_t npix = 300000; /* Target pixel count. */
    int32_t cols = (int32_t)ceil(2.0*sqrt(npix));
    int32_t rows = (int32_t)ceil(0.5*sqrt(npix));

    /* Choose array parameters: */
    ppv_size_t sz[d];
    sz[0] = cols;
    sz[1] = rows;
    if (d == 3) { sz[2] = chns; }

    /* Create array and fill with test pattern: */
    ppv_array_t *A = ppv_array_new(d, sz, maxsmp);
    ppv_print_descriptor (stderr, "A: ", A, "\n");
    pit_fill_array(A);
    
    /* Convert to image and write it out: */
    uint16_image_t *J = ppv_image_from_array(A);
    char *fname = NULL;
    asprintf(&fname, "out/J_d%d_chns%d_maxsmp%05d.png", d, chns, maxsmp);
    uint16_image_write_png_named(fname, J, 1.0, TRUE);
    free(fname);

    /* Convert back and compare: */
    ppv_array_t *B = ppv_image_to_array(J);
    ppv_print_descriptor (stderr, "B: ", B, "\n");
    assert(B->d == A->d);
    assert(B->size[0] == A->size[0]);
    assert(B->size[1] == A->size[1]);
    if (d == 3) assert(B->size[2] == A->size[2]);
    auto bool_t comp(const ppv_index_t *ix, ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(comp, FALSE, A, B, NULL);
    return;

    bool_t comp(const ppv_index_t *ix, ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t smpA = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        ppv_sample_t smpB = ppv_get_sample_at_pos(B->el, B->bps, B->bpw, pB);
        demand(smpA == smpB, "samples differ");
        return FALSE;
      }
  }
     
void pit_fill_array(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    ppv_sample_t maxsmp = A->maxsmp;
    if (d == 2)
      { /* Paint random grey balls: */
        ppv_throw_balls(A);
      }
    else if (d == 3)
      { /* Paint random balls in one channel and replicate: */
        int32_t chns = (int32_t)A->size[2];
        /* Take slice {C} of channel 0: */
        ppv_array_t *C = ppv_slice(A, 2, 0); 
        ppv_print_descriptor (stderr, "C: ", C, "\n");
        assert(C->d == 2);
    
        /* Fill channel 0 with balls: */
        ppv_throw_balls(C);
        /* Fill channels with versions of channel 0, leaving 0 for last: */
        for (int32_t ch = chns-1; ch >= 0; ch--)
          { /* Get channel {ch} as a 2D array {D}: */
            ppv_array_t *D = ppv_slice(A, 2, ch); 
            ppv_print_descriptor (stderr, "D: ", D, "");
            
            /* Pick a random mapping of {C}-samples to {D}-samples: */
            ppv_sample_t vmap[maxsmp+1];
            for (ppv_sample_t smp = 0; smp <= maxsmp; smp++)
              { vmap[smp] = (ppv_sample_t)uint32_abrandom(1, maxsmp); }
            
            /* Copy channel {C} (channel 0) to {D} (channel {ch}) with sample map: */
            auto bool_t chdup(const ix_index_t ix[], ix_pos_t pC, ix_pos_t pD, ix_pos_t pX);
            ppv_enum(chdup, FALSE, C, D, NULL);
            
            bool_t chdup(const ix_index_t ix[], ix_pos_t pC, ix_pos_t pD, ix_pos_t pX)
              { ppv_sample_t smpC = ppv_get_sample_at_pos(C->el, C->bps, C->bpw, pC);
                ppv_sample_t smpD = (smpC == 0 ? 0 : vmap[smpC]);
                ppv_set_sample_at_pos(D->el, D->bps, D->bpw, pD, smpD);
                return FALSE;
              }
          }
      }
    else
      { assert(FALSE); }
    return;
  }

