#define PROG_NAME "test_kdtom_const"
#define PROG_DESC "Test {kdtom_const.h} functions"
#define PROG_VERS "1.0"

#define tkdc_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-18 00:19:37 by jstolfi */

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
#include <ixbox.h>

#include <kdtom.h>
#include <kdtom_const.h>
#include <kdtom_test.h>

/* INTERNAL PROTOTYPES */

#define smpFMT ppv_sample_t_FMT

void tkdc_do_tests(ppv_dim_t d, ppv_sample_t maxsmp);

kdtom_const_t *tkdc_test_make(ppv_dim_t d, ppv_sample_t maxsmp, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Creates a {kdtom_const_t} node with given attributes. */

void tkdc_test_translate(kdtom_const_t *T);
  /* Tests {kdtom_translate(t,dx)} with some displacement {dx}. */

void tkdc_test_clip_core(kdtom_const_t *T);
  /* Tests {kdtom_clip_core(T,ixlo,size)} with various clip boxes.
    
    For each axis {ax}, the low and high indices of the clip box along
    that axis may be just below {T.DK}, at the low end of {T.DK},
    somewhere inside {T.DK}, at the upper end of {T.DK}, and just above
    that end.  The combinations tried will include one empty box.
    
    In any case, along axes other than {ax}, the clipbox will strictly contain 
    the core domain {T.DK}. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    int32_t nms = 33;
    ppv_sample_t *ms = kdtom_test_pick_max_samples(nms);
    
    for (ppv_dim_t d = 1; d <= 6; d++)
      { for (int32_t k = 0; k < nms; k++)
          { tkdc_do_tests(d, ms[k]); }
      }
      
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void tkdc_do_tests(ppv_dim_t d, ppv_sample_t maxsmp)
  { 
    fprintf(stderr, "=== testing d = %u maxsmp = " smpFMT " ===\n", d, maxsmp);
    
    /* Choose the core domain {ixlo,size}: */
    ppv_size_t size[d];
    ppv_sample_count_t npos = 10000;
    ppv_choose_test_size(d, npos, size);
    ppv_index_t ixlo[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ixlo[k] = 18 + 10*k; }

    kdtom_const_t *T = tkdc_test_make(d, maxsmp, ixlo, size);
    kdtom_print_node(stderr, 0, "T", (kdtom_t *)T, TRUE);
    assert(T->h.kind == kdtom_kind_CONST); 
    assert(T->h.d == d); 
    assert(T->h.maxsmp == maxsmp);
    
    tkdc_test_translate(T);
    tkdc_test_clip_core(T);
      
    return;
  }
    
void tkdc_test_translate(kdtom_const_t *T)
  {
    kdtom_test_translate((kdtom_t *)T);
    return;
  }

void tkdc_test_clip_core(kdtom_const_t *T)
  {
    ppv_dim_t d = T->h.d;
    
    /* Try clipping along all axes {r}: */
    for (ppv_axis_t r = 0; r < d;  r++) 
      { 
        auto void do_test_clip_core(ppv_index_t ixlo_r, ppv_size_t size_r);
          /* Tests clipping to range {ixlo_r .. ixlo_r + size_r - 1} on axis {r}. */
        kdtom_test_enum_ranges_single(T->h.ixlo[r], T->h.size[r], do_test_clip_core);

        continue;
        
        void do_test_clip_core(ppv_index_t ixlo_r, ppv_size_t size_r)
          { assert(T->h.d == d);
            kdtom_test_clip_core((kdtom_t *)T, r, ixlo_r, size_r);
          }
      }
    return;
  }

kdtom_const_t *tkdc_test_make(ppv_dim_t d, ppv_sample_t maxsmp, ppv_index_t ixlo[], ppv_size_t size[])
  {
    fprintf(stderr, "--- testing {kdtom_const_make} d = %u maxsmp = " smpFMT " ---\n", d, maxsmp);
    ixbox_print(stderr, d, "core = [", ixlo, size, " ]\n");
    
    ppv_sample_t fill = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));

    ppv_sample_t smp = (ppv_sample_t)(maxsmp <= 1 ? 0 : uint64_abrandom(1, maxsmp-1));
    if ((maxsmp > 0) && (smp == fill)) { smp++; }
 
    fprintf(stderr, "fill = " smpFMT " smp = " smpFMT "\n", fill, smp);
    assert(fill <= maxsmp);
    assert(smp <= maxsmp);
    if (maxsmp > 0) { assert(smp != fill); }

    kdtom_const_t *T = kdtom_const_make(d, maxsmp, fill, ixlo, size, smp);
    assert(T->h.kind == kdtom_kind_CONST); 
 
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
      /* Requires {ix} to be inside {T.DK}, and returns {T.V[ix]}. */

    kdtom_test_get_sample((kdtom_t *)T, getsmp);
    return T;

    auto ppv_sample_t getsmp(ppv_index_t ix[])
      { assert(ixbox_has(d, ix, ixlo, size));
        return smp;
      }

  }
    
