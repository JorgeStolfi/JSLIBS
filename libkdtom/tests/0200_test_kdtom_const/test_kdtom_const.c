#define PROG_NAME "test_kdtom_const"
#define PROG_DESC "Test {kdtom_const.h} functions"
#define PROG_VERS "1.0"

#define tkdc_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-12 22:32:08 by jstolfi */

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
  
void tkdc_test_get_sample(kdtom_const_t *T, ppv_index_t ixlo[], ppv_size_t size[]);
  /* Tests {kdtom_get_sample(T,ix)} with varous indices {ix} just inside and 
    just outside the core domain. */

void tkdc_test_translate(kdtom_const_t *T);
  /* Tests {kdtom_translate(t,dx)} with some displacement {dx}. */

void tkdc_test_clip(kdtom_const_t *T, sign_t loclip, sign_t hiclip);
  /* Tests {kdtom_clip(T,ixlo,size)} with various clip boxes.
    
    For each axis {ax}, as {loclip} ranges in {-2..+2}, the low index of
    the clip box will be just below {T.DK}, at the low end of {T.DK},
    inside {T.DK}, at the upper end of {T.DK}, and just above that end.
    
    The high index of the clip box along axis {ax} will be similarly
    selected as {hiclip} ranges in {-2..+2}.  Note that some combinations
    of {loclip} and {hiclip} may result in an empty box.
    
    Along axes other than {ax}, the clipbox will strictly contain 
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
    fprintf(stderr, "Testing d = %u maxsmp = " ppv_sample_t_FMT " ...\n", d, maxsmp);
    
    /* Choose the core domain {ixlo,size}: */
    ppv_size_t size[d];
    ppv_sample_count_t npos = 10000;
    ppv_choose_test_size(d, npos, size);
    ppv_index_t ixlo[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ixlo[k] = 18 + 10*k; }

    kdtom_const_t *T = tkdc_test_make(d, maxsmp, ixlo, size);

    tkdc_test_get_sample(T, ixlo, size);
    
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
    fprintf(stderr, "Testing {kdtom_translate} ...\n");
    ppv_dim_t d = T->h.d;
    
    bool_t empty = kdtom_has_empty_core((kdtom_t *)T);
    ppv_index_t ixlo_old[d]; /* Saved original {T.h.ixlo}. */
    ppv_index_t dx[d];       /* Displacement to apply. */
    ppv_index_t ixlo_new[d];  /* Expected {T.h.ixlo} after translation. */
    ppv_size_t size[d];      /* Size vector (should not change) */
    for (ppv_axis_t k = 0; k < d;  k++) 
      { size[k] = T->h.size[k];
        ixlo_old[k] = T->h.ixlo[k];
        dx[k] = 400 -10*k;
        ixlo_new[k] = (empty ? 0 : ixlo_old[k]  + dx[k]);
      }

    kdtom_translate((kdtom_t *)T, dx);
    
    /* Check: */
    for (ppv_axis_t k = 0; k < d; k++) 
      { demand(T->h.ixlo[k] == ixlo_new[k], "{kdtom_translate} bad {T.ixlo}"); 
        demand(T->h.size[k] == size[k], "{kdtom_translate} bad {T.size}"); 
      }
    
    tkdc_test_get_sample(T, ixlo_new, size);

    return;
  }

void tkdc_test_get_sample(kdtom_const_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    ppv_dim_t d = T->h.d;
    ppv_sample_t fill =T->h.fill;
    ppv_sample_t smp = T->smp;
    bool_t empty = kdtom_has_empty_core((kdtom_t *)T);

    ppv_index_t ix[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ix[k] = T->h.ixlo[k]; }

    ppv_sample_t val;
    for (ppv_axis_t k = 0; k < d;  k++)
      { val = kdtom_get_sample((kdtom_t *)T, ix);
        assert(val == (empty ? fill : smp));
        ix[k]--;
        val = kdtom_get_sample((kdtom_t *)T, ix);
        assert(val == fill);
        ix[k]++;
        ix[k] += size[k];
        val = kdtom_get_sample((kdtom_t *)T, ix);
        assert(val == fill);
        ix[k]--;
        val = kdtom_get_sample((kdtom_t *)T, ix);
        assert(val == (empty ? fill : smp));
        ix[k]++;
        ix[k] -= size[k];
      }
    return;
  }
  
void tkdc_test_clip(kdtom_const_t *T, sign_t loclip, sign_t hiclip)
  {
    fprintf(stderr, "Testing {kdtom_clip} loclip = %+d hiclip =%+d ...\n", loclip, hiclip);
    
    ppv_dim_t d = T->h.d;
    
    kdtom_test_show_box(stderr, d, "original box = [ ", T->h.ixlo, T->h.size, " ]\n");

    ppv_index_t ixlo_box[d]; /* {ixlo} value of clip box. */
    ppv_size_t size_box[d];  /* Size vector of clip box. */
    ppv_index_t ixlo_new[d]; /* Expected {T.h.ixlo} after clipping. */
    ppv_size_t size_new[d];  /* Expected {T.h.size} after clipping. */
    
    for (ppv_axis_t ax = 0; ax < d;  ax++) 
      { /* Get current core index range {ixlok_old..ixhik_old} on axis {ax}: */
        ppv_index_t ixlok_old = T->h.ixlo[ax];
        ppv_index_t ixhik_old = ixlok_old + T->h.size[ax] - 1;
                
        /* Select the clipping box index range {ixlok_box..ixhik_box} on axis {ax}: */
        ppv_index_t ixa[5] = { ixlok_old - 2, ixlok_old, ixlok_old + 1, ixhik_old, ixhik_old + 1};
        ppv_index_t ixb[5] = { ixlok_old - 1, ixlok_old, ixhik_old - 1, ixhik_old, ixhik_old + 2};
        ppv_index_t ixlok_box = ixa[loclip+2];
        ppv_index_t ixhik_box = ixb[hiclip+2];
        bool_t empty_box = (ixlok_box > ixhik_box);
        if (empty_box) { ixlok_box = 0; ixhik_box = -1; }
        ixlo_box[ax] = ixlok_box; 
        size_box[ax] = ixhik_box - ixlok_box + 1;
        
        /* Set index range along other axes as non-clipping: */
        for (ppv_axis_t k = 0; k < d;  k++)
          { if (k != ax)
              { ixlo_box[k] = (empty_box ? 0 : T->h.ixlo[k] - 1);
                size_box[k] = (empty_box ? 0 : T->h.size[k] + 2);
              }
          }
        kdtom_test_show_box(stderr, d, "clipping box = [ ", ixlo_box, size_box, " ]\n");
        
        /* Compute the expected clipped core box {ixlo_new[ax],size_new[ax]}: */
        kdtom_intersect_boxes(d, T->h.ixlo, T->h.size, ixlo_box, size_box, ixlo_new, size_new); 
        kdtom_test_show_box(stderr, d, "expected box = [ ", ixlo_new, size_new, " ]\n");

        kdtom_t *S = kdtom_clip((kdtom_t *)T, ixlo_box, size_box);
        kdtom_test_show_box(stderr, d, "clipped box =  [ ", S->ixlo, S->size, " ]\n");

        demand(S->d == T->h.d, "{kdtom_clip} bug: {S.d}");
        demand(S->maxsmp == T->h.maxsmp, "{kdtom_clip} bug: {S.maxsmp}");
        demand(S->fill == T->h.fill, "{kdtom_clip} bug: {S.fill}");
        for (ppv_axis_t k = 0; k < d;  k++)
          { demand(S->ixlo[k] == ixlo_new[k], "{kdtom_clip} bug: {S.ixlo}");
            demand(S->size[k] == size_new[k], "{kdtom_clip} bug: {S.size}");
          }
        
        demand(S->kind == kdtom_kind_CONST, "{kdtom_clip} bug: {S.kind}");
        tkdc_test_get_sample((kdtom_const_t *)S, ixlo_new, size_new);
        fprintf(stderr, "\n");
      }
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
    
    return T;
  }
    
