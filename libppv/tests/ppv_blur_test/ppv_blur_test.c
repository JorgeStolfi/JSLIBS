/* Last edited on 2021-07-03 09:33:06 by jstolfi */ 
/* Test of the {ppv_blur.h} module. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <wt_table.h>

#include <ppv_array.h>
#include <ppv_blur.h>

int main (int argn, char **argv)
  {
    ppv_dim_t dmax = 4;
    ppv_size_t rmax  = 2;

    for (ppv_dim_t d = 0; d <= dmax; d++)
      { for (ppv_size_t radius = 0; radius < rmax; radius++)
          { for (int32_t i = 0; i < 3; i++)
              { int32_t hstep;
                if (i == 0)
                  { hstep = 0; }
                else if (i == 1)
                  { hstep = radius; }
                else 
                  { /* Look for an allowed value of {hstep} other than 0 or {radius}: */
                    hstep = -1;
                    for (int32_t hh = 1; hh < radius; hh++)
                      { if ((2*radius+1) % (2*hh + 1) == 0) { hstep = hh; } }
                  }
                if (hstep >= 0) { pbt_do_test(d, radius, (ppv_size_t)hstep); } 
              }
          }
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
     
void pbt_do_test(ppv_dim_t d, ppv_size_t radius, ppv_size_t hstep)
  { fprintf(stderr, "--- testing {ppv_array_blur}");
    fprintf(stderr, " dimension %d radius = %d hstep = %d ---\n", d, radius, hstep);
    
    ppv_size_t t = 2*radius +1;
    double wt[t]; /* Hahn weight table, sums to 1. */
    wt_table_fill_hann((int)t, wt);
    
    ppv_array_t *A = pbt_make_array(d);
    
    /* Choose the number of bits per sample of {G}: */
    int32_t nw = ipow(t, d); /* Samples in neighborhood. */
    int32_t lognw = 0;
    while ((1 << lognw) < nw) { lognw++; }
    int32_t bps_min = lognw*A->bps;
    int32_t bps_max = imin(ppv_BPS_MAX, 2*bps_min + 1); 
    assert(bps_min <= bps_max);
    ppv_nbts_t bps = int32_abrandom(bps_min, bps_max);

    /* Blur array: */
    /* ??? Should test with custom {quantize,floatize}. ??? */
    ppv_array_t *G = ppv_array_blur(A, NULL, bps, radius, wt, hstep, NULL);

    /* ??? Should check and/or write out the array {G}. */
    return;
  }
        
ppv_array_t *pbt_make_array(ppv_dim_t d)
  { 
    ppv_nbits_t bps = (ppv_nbits_t)abrandom(1, 3);
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    ppv_size_t sz[d];
    ppv_choose_test_size(d, 1000000, sz);
    ppv_array_t *A = ppv_array_new(d, sz, bps, bpw);
    
    assert(ppv_blur_dimension(b) == d);
    ppv_sample_count_t nvtot = ppv_blur_voxel_count(b);
    fprintf(stderr, "{ppv_blur_voxel_count} says %ld voxels\n", nvtot);
    
    fprintf(stderr, "index ranges:");
    ppv_index_t ixlo[d], ixhi[d];
    ppv_blur_index_ranges(b, ixlo, ixhi);
    for (ppv_axis_t ax = 0; ax < d; ax++) 
      { fprintf(stderr, " {%+ld .. %+ld}", ixlo[ax], ixhi[ax]);
        demand(ixlo[ax] == -ixhi[ax], "ball brush index ranges are not symmetric");
        demand((rad-1 < ixhi[ax]) && (ixhi[ax] <= rad), "ball brush index range is bad");
      }
    fprintf(stderr, "\n");
    
    /* Check and count brush voxels through {ppv_blur_enum}: */
    ppv_sample_count_t nvtot_enum = 0; /* Voxels found by {ppv_blur_enum}. */
    
    auto bool_t ckvox(const ppv_index_t ix[]);
      /* Checks whether the center of voxel {ix} is at distance {rad} or less
        from the center of voxel {{0,0,..0}}. Also increments {nvtot_enum}.
        Always returns {FALSE}. */
    
    bool_t stop = ppv_blur_enum(ckvox, b);
    assert(! stop);
    fprintf(stderr, "{ppv_blur_enum} visited %ld voxels\n", nvtot_enum);
    demand(nvtot_enum == nvtot, "counts do not match");
    
    /* Check and count voxels that should be in brush by hand: */
    ppv_sample_count_t nvtry_hand = 0; /* Voxels tried by hand. */
    ppv_sample_count_t nvtot_hand = 0; /* Voxels found by hand. */
    ppv_index_t irad = (ppv_index_t)ceil(rad + 1.0); /* An upper bound on brush voxel indices. */
    fprintf(stderr, "checking all voxels with indices in {%+ld..%+ld}\n", -irad, +irad);
    ppv_index_t iradt[dmax];  /* Max abs index along each axis. */
    for (ppv_axis_t ax = 0; ax < dmax; ax++) { iradt[ax] = (ax < d ? irad : 0); }
    ppv_index_t ixt[dmax];
    assert(dmax == 4);
    for (ixt[0] = -iradt[0]; ixt[0] <= +iradt[0]; ixt[0]++)
      for (ixt[1] = -iradt[1]; ixt[1] <= +iradt[1]; ixt[1]++)
        for (ixt[2] = -iradt[2]; ixt[2] <= +iradt[2]; ixt[2]++)
          for (ixt[3] = -iradt[3]; ixt[3] <= +iradt[3]; ixt[3]++)
            { nvtry_hand++;
              double x[d];
              for (ppv_axis_t ax = 0; ax < d; ax++) { x[ax] = (double)ixt[ax]; }
              if (rn_norm(d, x) <= rad) { nvtot_hand++; }
            }
    fprintf(stderr, "hand enumeration tried %ld voxels found %ld in ball\n", nvtry_hand, nvtot_hand);
    demand(nvtot_hand == nvtot, "counts do not match");
    return;
    
    /* Internal procedures: */
    
    bool_t ckvox(const ppv_index_t ix[])
      { double x[d];
        for (ppv_axis_t ax = 0; ax < d; ax++) { x[ax] = (double)ix[ax]; }
        demand(rn_norm(d, x) <= rad, "ball voxel center is not in geometric ball");
        nvtot_enum++;
        return FALSE;
      }
  }

void test_blur_enum_stop(ppv_dim_t d, double rad)
  {
    fprintf(stderr, "--- testing {ppv_blur_enum} early stop for dimension %d radius = %20.16f ---\n", d, rad);
    
    ppv_dim_t dmax = 4;
    assert(d <= dmax);
    ppv_blur_t *b = ppv_blur_make_ball(d, rad);
    assert(ppv_blur_dimension(b) == d);
    ppv_sample_count_t nvtot = ppv_blur_voxel_count(b);
    assert(nvtot > 0);
    ppv_sample_count_t nvstop = (nvtot+1)/2;
    fprintf(stderr, "{ppv_blur_enum} should stop after visiting %ld voxels\n", nvstop);
    
    /* Check and count brush voxels through {ppv_blur_enum}: */
    ppv_sample_count_t nvtot_enum = 0; /* Voxels found by {ppv_blur_enum}. */
    
    auto bool_t ckvox(const ppv_index_t ix[]);
      /* Increments {nvtot_enum}, then returns {TRUE} iff {nvtot_enum >= nvstop}. */
    
    bool_t stop = ppv_blur_enum(ckvox, b);
    assert(stop);
    fprintf(stderr, "{ppv_blur_enum} stopped after visiting %ld voxels\n", nvtot_enum);
    demand(nvtot_enum == nvstop, "counts do not match");
    return;
    
    /* Internal procedures: */
    
    bool_t ckvox(const ppv_index_t ix[])
      { nvtot_enum++;
        return nvtot_enum >= nvstop;
      }
  }
