#define PROG_NAME "test_kdtom_split"
#define PROG_DESC "Test {kdtom_split.h} functions"
#define PROG_VERS "1.0"

#define tkds_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-16 21:46:06 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <ppv_array.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_test.h>

/* INTERNAL PROTOTYPES */

#define smpFMT ppv_sample_t_FMT
#define ixFMT ppv_index_t_FMT

void tkds_do_tests(ppv_dim_t d, ppv_sample_t maxsmp);

kdtom_split_t *tkds_test_make
  ( ppv_dim_t d, 
    ppv_sample_t maxsmp, 
    ppv_index_t ixlo[], 
    ppv_size_t size[],
    ppv_axis_t ax,
    ppv_size_t size0
  );
  /* Creates a {kdtom_split_t} node {T} with given attributes.
    The low sub-domains {T.DK0} will have size {size0}
    in the direction {ax}. Requires {size0} to be in {1..size[ax]-1}. */
  
void tkds_test_translate(kdtom_split_t *T);
  /* Tests {kdtom_translate(t,dx)} with some displacement {dx}. */

void tkds_test_clip_core(kdtom_split_t *T, ppv_axis_t ax, ppv_size_t size0);
  /* Tests {kdtom_clip_core(T,ixlo,size)} with various clip boxes.
    
    For each axis {r}, the low and high indices of the clip box along
    that axis may be just below {T.DK}, at the low end of {T.DK},
    somewhere inside {T.DK}, at the upper end of {T.DK}, and just above
    that end.
    
    When {r} is {ax}, additional values for the low and high indices
    of the clip box may be tried: just below {ixcut = T.h.ixlo[ax] + size0}, at
    {ixcut}, and just above {ixcut}.
    
    The combinations tried will include at least one empty box.
    
    In any case, along axes other than {r}, the clipbox will strictly contain 
    the core domain {T.DK}. */

kdtom_const_t *tkds_make_test_const
  ( ppv_dim_t d, 
    ppv_sample_t maxsmp, 
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[],
    ppv_sample_t smp
  );
  /* Creates a {kdtom_const_t} node {R} whose domain is usually overlapping
    with the box with low corner {ixlo[0..d-1]} and size {size[0..d-1]].
    The parameters {R.h.d}, {R.h.maxsmp}, {R.h.fill},{R.smp} will be as specified. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    int32_t nms = 33;
    ppv_sample_t *ms = kdtom_test_pick_max_samples(nms);
    
    for (ppv_dim_t d = 1; d <= 6; d++)
      { for (int32_t k = 0; k < nms; k++)
          { tkds_do_tests(d, ms[k]); }
      }
      
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void tkds_do_tests(ppv_dim_t d, ppv_sample_t maxsmp)
  { 
    fprintf(stderr, "=== testing d = %u maxsmp = " smpFMT " ===\n", d, maxsmp);
    
    /* Choose the core domain {ixlo,size}: */
    ppv_size_t size[d];
    ppv_sample_count_t npos = 10000;
    ppv_choose_test_size(d, npos, size);
    ppv_index_t ixlo[d];
    for (ppv_axis_t k = 0; k < d;  k++) 
      { ixlo[k] = 18 + 10*k;
        assert(size[k] > 0);
      }
      
    /* Choose the split axis {ax} and size of low sub-domain {size0}: */
    ppv_axis_t ax = (ppv_axis_t)uint32_abrandom(0, d-1);
    if (size[ax] < 6) { size[ax] = 6; }
    ppv_size_t size0 = (ppv_size_t)uint64_abrandom(1, size[ax]-1);

    kdtom_split_t *T = tkds_test_make(d, maxsmp, ixlo, size, ax, size0);
    kdtom_print_node(stderr, 0, "T", (kdtom_t *)T, TRUE);
    assert(T->h.kind == kdtom_kind_SPLIT); 
    assert(T->h.d == d); 
    assert(T->h.maxsmp == maxsmp);

    tkds_test_translate(T);
    
    tkds_test_clip_core(T, ax, size0);
      
    return;
  }
    
void tkds_test_translate(kdtom_split_t *T)
  {
    kdtom_test_translate((kdtom_t *)T);
    return;
  }

void tkds_test_clip_core(kdtom_split_t *T, ppv_axis_t ax, ppv_size_t size0)
  {
    ppv_dim_t d = T->h.d;

    /* Try clipping along all axes {r}: */
    for (ppv_axis_t r = 0; r < d;  r++) 
      { 
        auto void do_test_clip_core(ppv_index_t ixlo_r, ppv_size_t size_r);
          /* Tests clipping to range {ixlo_r .. ixlo_r + size_r - 1} on axis {r}. */
          
        if (r == ax)
          { kdtom_test_enum_ranges_split(T->h.ixlo[r], T->h.size[r], size0, do_test_clip_core); }
        else
          { kdtom_test_enum_ranges_single(T->h.ixlo[r], T->h.size[r], do_test_clip_core); }

        continue;

        void do_test_clip_core(ppv_index_t ixlo_r, ppv_size_t size_r)
          { kdtom_test_clip_core((kdtom_t *)T, r, ixlo_r, size_r);  }
      }
    return;
  }

kdtom_split_t *tkds_test_make
  ( ppv_dim_t d, 
    ppv_sample_t maxsmp, 
    ppv_index_t ixlo[], 
    ppv_size_t size[],
    ppv_axis_t ax, 
    ppv_size_t size0
  )
  {
    fprintf(stderr, "--- testing {kdtom_split_make} d = %u maxsmp = " smpFMT " ---\n", d, maxsmp);
    ixbox_print(stderr, d, "core = [", ixlo, size, " ]\n");
     
    fprintf(stderr, "ax = %u size0 = " ppv_size_t_FMT "\n", ax, size0);
    assert((size0 > 0) && (size0 < size[ax]));
    ppv_size_t size1 = size[ax] - size0;
    
    /* Choose the fill sample values: */
    ppv_sample_t fill = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));
    ppv_sample_t fill0 = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));
    ppv_sample_t fill1 = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));
    if (uint32_abrandom(0, 2) == 0) { fill0 = fill; }
    if (uint32_abrandom(0, 2) == 0) { fill1 = fill; }

    fprintf(stderr, "fill = " smpFMT, fill);
    assert(fill <= maxsmp);
    fprintf(stderr, " fill0 = " smpFMT, fill0);
    assert(fill0 <= maxsmp);
    fprintf(stderr, " fill1 = " smpFMT "\n", fill1);
    assert(fill1 <= maxsmp);
    
    /* Vector to translate the two sub-nodes to the origin: */
    ppv_index_t dx[d];
    for (ppv_axis_t k = 0; k < d;  k++) { dx[k] = -ixlo[k]; }

    ppv_index_t iax = ixlo[ax];  /* Save given {ixlo[ax]}. */
    ppv_size_t sax = size[ax];   /* Save given {size[ax]}. */

    /* Create the two sub-nodes {R0,R1} around {T.DK0} and {T,DK1}: */
    
    ixlo[ax] = iax;
    size[ax] = size0;
    kdtom_const_t *R0 = tkds_make_test_const(d, maxsmp, fill0, ixlo, size, maxsmp);

    ixlo[ax] = iax + size0;
    size[ax] = size1;
    kdtom_const_t *R1 = tkds_make_test_const(d, maxsmp, fill1, ixlo, size, maxsmp);

    /* Create the clipped versions {T0,T1} shifted down to the origin: */
    ixlo[ax] = iax;
    size[ax] = size0;
    kdtom_t *T0 = kdtom_clip_core((kdtom_t *)R0, ixlo, size);
    dx[ax] = -ixlo[ax];
    kdtom_translate(T0, dx);
    
    ixlo[ax] = iax + size0;
    size[ax] = size1;
    kdtom_t *T1 = kdtom_clip_core((kdtom_t *)R1, ixlo, size);
    dx[ax] = -ixlo[ax];
    kdtom_translate(T1, dx);

    /* Restore the given {ixlo,size}: */
    ixlo[ax] = iax;
    size[ax] = sax;

    /* Build the SPLIT node: */
    kdtom_split_t *T = kdtom_split_make(d, maxsmp, fill, ixlo, size, ax, T0, size0, T1, size1);
    assert(T->h.kind == kdtom_kind_SPLIT); 
    assert(T->h.d == d); 
    demand(T->h.maxsmp == maxsmp, "wrong {maxsmp}");
    demand(T->h.fill == fill, "wrong {fill}");
    for (ppv_axis_t k = 0; k < d;  k++) 
      { demand(T->h.ixlo[k] == ixlo[k], "wrong {ixlo}");
        demand(T->h.size[k] == size[k], "wrong {size}");
      }
    demand(T->size0 == size0, "wrong {size0}");
    
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
      /* Requires {ix} to be inside {T.DK}, and returns the expected {T.V[ix]}. */

    kdtom_test_get_sample((kdtom_t *)T, getsmp);

    return T;

    auto ppv_sample_t getsmp(ppv_index_t ix[])
      { assert(ixbox_has(d, ix, ixlo, size));
        ppv_sample_t smp;
        if (ix[ax] < ixlo[ax] + size0)
          { smp = kdtom_get_sample((kdtom_t *)R0, ix); }
        else
          { smp = kdtom_get_sample((kdtom_t *)R1, ix); }
        return smp; 
      }
  }
    
kdtom_const_t *tkds_make_test_const
  ( ppv_dim_t d, 
    ppv_sample_t maxsmp, 
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[],
    ppv_sample_t smp
  )
  {
    ppv_index_t ixlo_R[d];
    ppv_size_t size_R[d];
    for (ppv_axis_t k = 0; k < d;  k++) 
      { ppv_index_t szm = size[k]/2;
        ppv_index_t ixlo_min = ixlo[k] - 3;
        ppv_index_t ixlo_max = ixlo[k] + szm - 1;
        ppv_index_t ixlok = (ppv_index_t)int64_abrandom(ixlo_min, ixlo_max);
        ppv_index_t ixhi_min = ixlo[k] + szm;
        ppv_index_t ixhi_max = ixlo[k] + size[k] + 2;
        ppv_index_t ixhik = (ppv_index_t)int64_abrandom(ixhi_min, ixhi_max);
        ixlo_R[k] = ixlok;
        size_R[k] = ixhik - ixlok + 1;
      }
    
    kdtom_const_t *R = kdtom_const_make(d, maxsmp, fill, ixlo_R, size_R, smp);
    return R;
  }
