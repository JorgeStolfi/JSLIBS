/* Last edited on 2023-03-18 10:57:22 by stolfi */ 
/* Test of the {ppv_brush.h} module. */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <rn.h>

#include <ppv_array.h>
#include <ppv_brush.h>

#define bug(MSG) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fputs(MSG, stderr); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

void test_brush_make_ball(ppv_dim_t d, double rad);
void test_brush_enum_stop(ppv_dim_t d, double rad);

int32_t main (int32_t argn, char **argv)
  {
    ppv_dim_t dmax = 4;
    double rad = 3.1416; /* A random radius for test */
    for (ppv_dim_t d = 0; d <= dmax; d++)
      { test_brush_make_ball(d, rad);
        test_brush_make_ball(d, 0.0);
        test_brush_enum_stop(d, rad);
        test_brush_enum_stop(d, rad);
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
        
void test_brush_make_ball(ppv_dim_t d, double rad)
  { fprintf(stderr, "--- testing {ppv_brush_make_ball} for dimension %d radius = %20.16f ---\n", d, rad);
    
    ppv_dim_t dmax = 4;
    assert(d <= dmax);
    ppv_brush_t *b = ppv_brush_make_ball(d, rad);
    
    assert(ppv_brush_dimension(b) == d);
    ppv_sample_count_t nvtot = ppv_brush_voxel_count(b);
    fprintf(stderr, "{ppv_brush_voxel_count} says %ld voxels\n", nvtot);
    
    fprintf(stderr, "index ranges:");
    ppv_index_t ixlo[d], ixhi[d];
    ppv_brush_index_ranges(b, ixlo, ixhi);
    for (ppv_axis_t ax = 0; ax < d; ax++) 
      { fprintf(stderr, " {%+ld .. %+ld}", ixlo[ax], ixhi[ax]);
        demand(ixlo[ax] == -ixhi[ax], "ball brush index ranges are not symmetric");
        demand((rad-1 < ixhi[ax]) && (ixhi[ax] <= rad), "ball brush index range is bad");
      }
    fprintf(stderr, "\n");
    
    /* Check and count brush voxels through {ppv_brush_enum}: */
    ppv_sample_count_t nvtot_enum = 0; /* Voxels found by {ppv_brush_enum}. */
    
    auto bool_t ckvox(const ppv_index_t ix[]);
      /* Checks whether the center of voxel {ix} is at distance {rad} or less
        from the center of voxel {{0,0,..0}}. Also increments {nvtot_enum}.
        Always returns {FALSE}. */
    
    bool_t stop = ppv_brush_enum(ckvox, b);
    assert(! stop);
    fprintf(stderr, "{ppv_brush_enum} visited %ld voxels\n", nvtot_enum);
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

void test_brush_enum_stop(ppv_dim_t d, double rad)
  {
    fprintf(stderr, "--- testing {ppv_brush_enum} early stop for dimension %d radius = %20.16f ---\n", d, rad);
    
    ppv_dim_t dmax = 4;
    assert(d <= dmax);
    ppv_brush_t *b = ppv_brush_make_ball(d, rad);
    assert(ppv_brush_dimension(b) == d);
    ppv_sample_count_t nvtot = ppv_brush_voxel_count(b);
    assert(nvtot > 0);
    ppv_sample_count_t nvstop = (nvtot+1)/2;
    fprintf(stderr, "{ppv_brush_enum} should stop after visiting %ld voxels\n", nvstop);
    
    /* Check and count brush voxels through {ppv_brush_enum}: */
    ppv_sample_count_t nvtot_enum = 0; /* Voxels found by {ppv_brush_enum}. */
    
    auto bool_t ckvox(const ppv_index_t ix[]);
      /* Increments {nvtot_enum}, then returns {TRUE} iff {nvtot_enum >= nvstop}. */
    
    bool_t stop = ppv_brush_enum(ckvox, b);
    assert(stop);
    fprintf(stderr, "{ppv_brush_enum} stopped after visiting %ld voxels\n", nvtot_enum);
    demand(nvtot_enum == nvstop, "counts do not match");
    return;
    
    /* Internal procedures: */
    
    bool_t ckvox(const ppv_index_t ix[])
      { nvtot_enum++;
        return nvtot_enum >= nvstop;
      }
  }
