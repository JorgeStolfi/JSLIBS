#define PROG_NAME "test_kdtom_const"
#define PROG_DESC "Test {kdtom_const.h} functions"
#define PROG_VERS "1.0"

#define tkdc_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-13 02:39:18 by jstolfi */

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
#include <kdtom_const.h>
#include <kdtom_test.h>

/* INTERNAL PROTOTYPES */

void tkdc_do_tests(ppv_dim_t d, ppv_sample_t maxsmp);

kdtom_const_t *tkdc_test_make(ppv_dim_t d, ppv_sample_t maxsmp, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Creates a {kdtom_const_t} node with given attributes. */

void tkdc_test_translate(kdtom_const_t *T);
  /* Tests {kdtom_translate(t,dx)} with some displacement {dx}. */

void tkdc_test_clip(kdtom_const_t *T, int32_t loclip, int32_t hiclip);
  /* Tests {kdtom_clip(T,ixlo,size)} with various clip boxes.
    See {kdtom_test_clip} for the meaning of {loclip,hiclip}
    and the clipping boxes that are used.  */

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
    fprintf(stderr, "Testing d = %u maxsmp = " ppv_sample_t_FMT " ...\n", d, maxsmp);
    
    /* Choose the core domain {ixlo,size}: */
    ppv_size_t size[d];
    ppv_sample_count_t npos = 10000;
    ppv_choose_test_size(d, npos, size);
    ppv_index_t ixlo[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ixlo[k] = 18 + 10*k; }

    kdtom_const_t *T = tkdc_test_make(d, maxsmp, ixlo, size);
    
    tkdc_test_translate(T);
    
    for (int32_t loclip = -2; loclip <= +2; loclip++)
      { for (int32_t hiclip = -2; hiclip <= +2; hiclip++)
          { if ((loclip == -1) || (loclip <= hiclip))
              { tkdc_test_clip(T, loclip, hiclip); }
          }
      }
      
    return;
  }
    
void tkdc_test_translate(kdtom_const_t *T)
  {
    kdtom_test_translate((kdtom_t *)T);
    return;
  }

void tkdc_test_clip(kdtom_const_t *T, int32_t loclip, int32_t hiclip)
  {
    kdtom_test_clip((kdtom_t *)T, loclip, hiclip);
    return;
  }

kdtom_const_t *tkdc_test_make(ppv_dim_t d, ppv_sample_t maxsmp, ppv_index_t ixlo[], ppv_size_t size[])
  {
    fprintf(stderr, "Testing {kdtom_const_make} d = %u maxsmp = " ppv_sample_t_FMT " ...\n", d, maxsmp);
    fprintf(stderr, "ixlo =");
    for (ppv_axis_t k = 0; k < d;  k++) { fprintf(stderr, " " ppv_index_t_FMT, ixlo[k]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "size =");
    for (ppv_axis_t k = 0; k < d;  k++) { fprintf(stderr, " " ppv_size_t_FMT, size[k]); }
    fprintf(stderr, "\n");
    
    ppv_sample_t fill = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));

    ppv_sample_t smp = (ppv_sample_t)(maxsmp <= 1 ? 0 : uint64_abrandom(1, maxsmp-1));
    if ((maxsmp> 0) && (smp == fill)) { smp++; }
 
    fprintf(stderr, "fill = " ppv_sample_t_FMT " smp = " ppv_sample_t_FMT "\n", fill, smp);
    assert(fill <= maxsmp);
    assert(smp <= maxsmp);
    if (maxsmp > 0) { assert(smp != fill); }

    kdtom_const_t *T = kdtom_const_make(d, maxsmp, fill, ixlo, size, smp);
    assert(T->h.kind == kdtom_kind_CONST); 
 
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
      /* Requires {ix} to be inside {T.DK}, and returns {T.V[ix]}. */

    kdtom_test_get_sample(T, getsmp);
    return T;

    auto ppv_sample_t getsmp(ppv_index_t ix[])
      { assert(kdtom_index_is_in_box(d, ix, ixlo, size));
        return smp;
      }

  }
    
