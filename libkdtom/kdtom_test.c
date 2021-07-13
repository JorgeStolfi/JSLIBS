/* See {kdtom_test.h}. */
/* Last edited on 2021-07-13 02:57:39 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsrandom.h>
#include <ppv_array.h>
#include <jsmath.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_array.h>
#include <kdtom_const.h>

#include <kdtom_test.h>

void kdtom_test_get_sample(kdtom_t *T, kdtom_test_get_sample_proc_t *getsmp)
  {
    ppv_dim_t d = T->d;
    ppv_sample_t fill = T->fill;
    bool_t empty = kdtom_has_empty_core((kdtom_t *)T);

    ppv_index_t ix[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ix[k] = T->ixlo[k]; }

    ppv_sample_t val_T; /* Sample value fetched from {T}. */
    ppv_sample_t val_X; /* Sample value obtained from {getsmp}. */
    if (empty)
      { val_T = kdtom_get_sample((kdtom_t *)T, ix); 
        assert(val_T == fill);
        for (ppv_axis_t k = 0; k < d;  k++)
          { ix[k]--;
            val_T = kdtom_get_sample((kdtom_t *)T, ix);
            assert(val_T == fill);
            ix[k]++;

            ix[k]++;
            val_T = kdtom_get_sample((kdtom_t *)T, ix);
            assert(val_T == fill);
            ix[k]--;
          }
      }
    else 
      { val_T = kdtom_get_sample((kdtom_t *)T, ix); val_X = getsmp(ix);
        assert(val_T == val_X);

        for (ppv_axis_t k = 0; k < d;  k++)
          { ix[k]--;
            val_T = kdtom_get_sample((kdtom_t *)T, ix);
            assert(val_T == fill);
            ix[k]++;

            ix[k] += T->size[k];

            val_T = kdtom_get_sample((kdtom_t *)T, ix);
            assert(val_T == fill);

            ix[k]--;
            val_T = kdtom_get_sample((kdtom_t *)T, ix); val_X = getsmp(ix);
            assert(val_T == val_X);
            ix[k]++;

            ix[k] -= T->size[k];
          }
      }
    return;
  }

void kdtom_test_clip(kdtom_t *T, int32_t loclip, int32_t hiclip)
  {
    fprintf(stderr, "Testing {kdtom_clip} loclip = %+d hiclip =%+d ...\n", loclip, hiclip);
    
    ppv_dim_t d = T->d;
    
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
      /* Requires {ix} to be inside {T.DK}, and returns {T.V[ix]}. */

    kdtom_test_show_box(stderr, d, "original box = [ ", T->ixlo, T->size, " ]\n");

    ppv_index_t ixlo_box[d]; /* {ixlo} value of clip box. */
    ppv_size_t size_box[d];  /* Size vector of clip box. */
    ppv_index_t ixlo_new[d]; /* Expected {T.h.ixlo} after clipping. */
    ppv_size_t size_new[d];  /* Expected {T.h.size} after clipping. */
    
    for (ppv_axis_t ax = 0; ax < d;  ax++) 
      { /* Get current core index range {ixlok_old..ixhik_old} on axis {ax}: */
        ppv_index_t ixlok_old = T->ixlo[ax];
        ppv_index_t ixhik_old = ixlok_old + T->size[ax] - 1;
                
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
              { ixlo_box[k] = (empty_box ? 0 : T->ixlo[k] - 1);
                size_box[k] = (empty_box ? 0 : T->size[k] + 2);
              }
          }
        kdtom_test_show_box(stderr, d, "clipping box = [ ", ixlo_box, size_box, " ]\n");
        
        /* Compute the expected clipped core box {ixlo_new[ax],size_new[ax]}: */
        kdtom_intersect_boxes(d, T->ixlo, T->size, ixlo_box, size_box, ixlo_new, size_new); 
        kdtom_test_show_box(stderr, d, "expected box = [ ", ixlo_new, size_new, " ]\n");

        kdtom_t *S = kdtom_clip((kdtom_t *)T, ixlo_box, size_box);
        kdtom_test_show_box(stderr, d, "clipped box =  [ ", S->ixlo, S->size, " ]\n");

        demand(S->d == T->d, "{kdtom_clip} bug: {S.d}");
        demand(S->maxsmp == T->maxsmp, "{kdtom_clip} bug: {S.maxsmp}");
        demand(S->fill == T->fill, "{kdtom_clip} bug: {S.fill}");
        for (ppv_axis_t k = 0; k < d;  k++)
          { demand(S->ixlo[k] == ixlo_new[k], "{kdtom_clip} bug: {S.ixlo}");
            demand(S->size[k] == size_new[k], "{kdtom_clip} bug: {S.size}");
          }
        
        kdtom_test_get_sample(S, getsmp);
        fprintf(stderr, "\n");
        free(S);
      }

    return;
    
    auto ppv_sample_t getsmp(ppv_index_t ix[])
      { 
        assert(kdtom_index_is_in_box(d, ix, T->ixlo, T->size));
        return kdtom_get_sample(T, ix);
      }
  }

void kdtom_test_translate(kdtom_t *T)
  {
    fprintf(stderr, "Testing {kdtom_translate} ...\n");
    ppv_dim_t d = T->d;
    kdtom_test_show_box(stderr, d, "T core = [ ", T->ixlo, T->size, " ]\n");
    
    bool_t empty = kdtom_has_empty_core((kdtom_t *)T);
    ppv_index_t dx[d];       /* Displacement to apply. */
    fprintf(stderr, "displ = [");
    for (ppv_axis_t k = 0; k < d;  k++) 
      { dx[k] = 400 -10*k;
        fprintf(stderr, " " ppv_index_t_FMT, dx[k]);
      }
    fprintf(stderr, " ]\n");

    kdtom_t *S = kdtom_clone(T);
    
    kdtom_translate((kdtom_t *)S, dx);
    kdtom_test_show_box(stderr, d, "S core = [", S->ixlo, S->size, " ]\n");
    
    /* Check: */
    for (ppv_axis_t k = 0; k < d; k++) 
      { demand(S->ixlo[k] == (empty ? 0 : T->ixlo[k] + dx[k]), "{kdtom_translate} bad {T.ixlo}"); 
        demand(S->size[k] == (empty ? 0 : T->size[k]), "{kdtom_translate} bad {T.size}"); 
      }
    
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
    kdtom_test_get_sample(S, getsmp);
    
    free(S);

    return;
    
    ppv_sample_t getsmp(ppv_index_t ix[])
      { ppv_index_t jx[d];
        for (ppv_axis_t k = 0; k < d; k++) { jx[k] = ix[k] - dx[k]; }
        assert(kdtom_index_is_in_box(T->d, jx, T->ixlo, T->size));
        return kdtom_get_sample(T, jx);
      }
  }

ppv_sample_t *kdtom_test_pick_max_samples(int32_t nms)
  { ppv_sample_t *ms = notnull(malloc(nms*sizeof(ppv_sample_t)), "no mem");
    for (int32_t k = 0; k < nms; k++)
      { ppv_nbits_t bps = (ppv_nbits_t)floor(32*((double)k)/(nms-1) + 0.5);
        if (bps == 0)
          { ms[k] = 0; }
        else
          { ppv_sample_t maxmaxsmp = ppv_max_sample(bps);
            ms[k] = (ppv_sample_t)(uint64_abrandom(maxmaxsmp/2+1, maxmaxsmp));
            assert(ppv_min_bps(ms[k]) == bps);
          }
        fprintf(stderr, " bps = %2d  ms[%2d] = " ppv_sample_t_FMT "\n", bps, k, ms[k]);
      }
    return ms;
  }

ppv_array_t *kdtom_test_array_make(ppv_dim_t d, ppv_size_t size[], ppv_sample_t maxsmp)
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");

    ppv_array_t *A = ppv_array_new(d, size, maxsmp);
    
    if (maxsmp == 0) { /* Array must be all zeros anyway: */ return A; }
    
    /* Fill the array with random samples: */
    srandom(4615);
    ppv_throw_noise(A); 

    if (d == 0) { /* Only one sample: */  return A; }

    /* Compute min size {szmin} and cords {ctr} of domain center: */
    ppv_size_t szmin = ppv_MAX_SIZE;
    ppv_size_t szmax = 0;
    double ctr[d]; /* Center of bullseye. */
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_size_t szk = A->size[k];
        if (szk < szmin) { szmin = szk; }
        if (szk > szmax) { szmax = szk; }
        ctr[k] = 0.5*(double)szk;
      }
    if (szmax == 0) { /* Empty aray: */ return A; }
    if (szmin >= 3) 
      { /* Paint the test pattern: */
        double R = 0.50*(double)szmin; /* Radius of max inscribed ball. */
        kdtom_test_paint_bullseye(A, ctr, R);
      }
    return A;
  }
    
void kdtom_test_paint_bullseye(ppv_array_t *A, double ctr[], double R)
  { 
    ppv_dim_t d = A->d;
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid array dimension {d}");
    ppv_sample_t maxsmp = A->maxsmp;
    
    /* Bullseye raddi (may be negative): */
    double R3 = R;           /* Outer radius of outer noise ring. */
    double R2 = R3-sqrt(d);  /* Inner radius of outer noise ring. */
    double R1 = 0.5*R2;      /* Outer radius of inner noise ring. */
    double R0 = R1-sqrt(d);  /* Inner radius of inner noise ring. */

    assert(maxsmp > 0);
    ppv_sample_t smpA = 0;
    ppv_sample_t smpB = (ppv_sample_t)(maxsmp+1)/2;
    assert(smpA != smpB);

    auto bool_t bullpaint(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
      /* Sets the area within {R0} and outside {R3} to {smpA}.
        Sets the area between {R1} and {R2} to {smpB}.
        Leaves other voxels unchanged. */
    ppv_enum(bullpaint, FALSE, A, NULL, NULL);

    /* Invert the center sample: */
    ppv_index_t ctrix[d];
    for (ppv_axis_t k = 0; k < d; k++) 
      { ctrix[k] = A->size[k]/2;
        assert(ctrix[k] < A->size[k]);
      }
    ppv_sample_t ctrval = ppv_get_sample(A, ctrix);
    assert(ctrval <= maxsmp);
    ppv_sample_t newval = (maxsmp - ctrval < ctrval ? 0 : maxsmp);
    ppv_set_sample(A, ctrix, newval);
    return;

    /* Internal proc implementations: */

    bool_t bullpaint(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      {
        /* Compute squared distance from pixel center to domain center: */
        double dist2 = 0;
        for (ppv_axis_t k = 0; k < d; k++)
          { double dk = ((double)ix[k]) + 0.5 - ctr[k];
            dist2 += dk*dk;
          }
        double dist = sqrt(dist2);
        if ((dist <= R0) || (dist >= R3))
          { ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, smpA); }
        else if ((dist >= R1) && (dist <= R2))
          { ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, smpB); }
        else
          { /* Leave sample unchaged. */ }
        return FALSE;
      }
  }

void kdtom_test_choose_array_size(ppv_dim_t d, ppv_size_t size[])
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");
    ppv_sample_count_t npos = 1;
    if (d > 0)
      { /* Choose relative sizes: */
        double rsz_prod = 1;
        double rsz[d];
        for (ppv_axis_t k = 0; k <d; k++)
          { rsz[k] = pow(3.0, ((double)k)/d);
            rsz_prod *= rsz[k];
          }
        double rsz_avg = pow(rsz_prod, 1.0/d);
        /* Compute actual sizes: */
        double xnpos = 1000000.0;         /* Expected number of samples. */
        double xsize = pow(xnpos, 1.0/d); /* Expected size per axis. */
        fprintf(stderr, "array size:");
        for (ppv_axis_t k = 0; k <d; k++)
          { size[k] = (ppv_size_t)floor(xsize*rsz[k]/rsz_avg);
            fprintf(stderr, " %lu", size[k]);
            npos *= size[k];
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "array has %lu samples\n", npos);
    return;
  }

void kdtom_test_plot(kdtom_t *T)
  {
    demand(T->d == 2, "Can't plot this");
    fprintf(stderr, "Pretend that I am plotting...\n");    
    fprintf(stderr, "*plot*, *plot*, *plit*, *plut*...\n");
    return;
  }

void kdtom_test_check_tree(kdtom_t *T)
  {
    ppv_dim_t d = T->d;
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");

    size_t rec_bytes = kdtom_bytesize(T, FALSE);
    size_t tot_bytes = kdtom_bytesize(T, TRUE);
    fprintf(stderr, "byte sizes: root record = %lu  whole tree = %lu\n", rec_bytes, tot_bytes);
    
    ppv_sample_t maxsmp = T->maxsmp;
    demand(maxsmp <= ppv_MAX_SAMPLE_VAL, "invalid {maxsmp}");
    ppv_sample_t fill = T->fill;
    demand(fill <= maxsmp, "invalid {fill}");
    
    /* Either all axes have zero size, none has: */
    for (ppv_axis_t k = 0; k <d; k++)
      { demand((T->size[k] == 0) == (T->size[0] == 0), "inconsisten empty {size}"); }
    
    fprintf(stderr, "!! {kdtom_test_check_tree} incomplete");
    return;
  }

void kdtom_test_show_box(FILE *wr, ppv_dim_t d, char *pref, ppv_index_t ixlo[], ppv_size_t size[], char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    for (ppv_axis_t k = 0; k <d; k++)
      { if (k > 0) { fputs(" ", wr); }
        fprintf(wr, ppv_index_t_FMT "(" ppv_size_t_FMT ")", ixlo[k], size[k]);
      }
    if (suff != NULL) {fputs(suff, wr); }
    return;
  }
