/* See {kdtom_test.h}. */
/* Last edited on 2024-12-05 10:33:00 by stolfi */

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
#include <ix.h>
#include <ixbox.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_array.h>
#include <kdtom_const.h>

#include <kdtom_test.h>

#define ixFMT ppv_index_t_FMT
#define smpFMT ppv_sample_t_FMT
#define szFMT ppv_size_t_FMT

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
            val_T = kdtom_get_sample((kdtom_t *)T, ix);
            val_X = getsmp(ix);
            assert(val_T == val_X);
            ix[k]++;

            ix[k] -= T->size[k];
          }
      }
    return;
  }

int32_t ktdom_test_remove_dup_indices(int32_t n, ppv_index_t ix[])
  { int32_t nu = 0;
    for (uint32_t k = 0;  k < n; k++)
      { ppv_index_t ixk = ix[k];
        /* Search for the proper place of {ixk}: */
        int32_t i = nu; while ((i > 0) && (ix[i-1] > ixk)) { i--; }
        if ((i <= 0) || (ix[i-1] < ixk))
          { /* Not repeated.  Insert {ixk} in proper place: */
            i = nu; while ((i > 0) && (ix[i-1] > ixk)) { ix[i] = ix[i-1]; i--; }
            ix[i] = ixk; nu++;
          }
      }
    return nu;
  }

void kdtom_test_enum_ranges_single(ppv_index_t ixlo, ppv_size_t size, kdtom_test_range_proc_t process)
  {
    /* Interesting index values for low and high indices of clipping box: */
    #define NH 5
    ppv_index_t ixhi = ixlo + size - 1;
    ppv_index_t ixlo_hot[NH] = { ixlo - 2, ixlo, ixlo + 1, ixhi, ixhi + 1};
    ppv_index_t ixhi_hot[NH] = { ixlo - 1, ixlo, ixhi - 1, ixhi, ixhi + 2};
    /* Remove duplicate interesting values: */
    int32_t nlo = ktdom_test_remove_dup_indices(NH, ixlo_hot);
    int32_t nhi = ktdom_test_remove_dup_indices(NH, ixhi_hot);
    kdtom_test_enum_ranges(nlo, ixlo_hot, nhi, ixhi_hot, process);
    return;
    #undef NH
  }
    
void kdtom_test_enum_ranges_split
  ( ppv_index_t ixlo,
    ppv_size_t size, 
    ppv_size_t size0,
    kdtom_test_range_proc_t process
  )
  {
    /* Interesting index values for low and high indices of clipping box: */
    #define NH 7
    ppv_index_t ixhi = ixlo + size - 1;
    ppv_index_t ixlo_hot[NH] = 
      { ixlo - 2, ixlo, ixlo + 1, 
        ixlo + size0 - 2, ixlo + size0,           
        ixhi, ixhi + 1
      };
    ppv_index_t ixhi_hot[NH] = 
      { ixlo - 1, ixlo,           
        ixlo + size0 - 1, ixlo + size0 + 1, 
        ixhi - 1, ixhi, ixhi + 2
      };
    /* Remove duplicate interesting values: */
    int32_t nlo = ktdom_test_remove_dup_indices(NH, ixlo_hot);
    int32_t nhi = ktdom_test_remove_dup_indices(NH, ixhi_hot);
    kdtom_test_enum_ranges(nlo, ixlo_hot, nhi, ixhi_hot, process);
    return;
    #undef NH
  }

void kdtom_test_enum_ranges
  ( int32_t nlo, 
    ppv_index_t ixlo_hot[], 
    int32_t nhi, 
    ppv_index_t ixhi_hot[], 
    kdtom_test_range_proc_t process
  )
  {
    fprintf(stderr, "ixlo_hot =");
    for (uint32_t lohot = 0;  lohot < nlo; lohot++) { fprintf(stderr, " " ixFMT, ixlo_hot[lohot]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "ixhi_hot =");
    for (uint32_t hihot = 0;  hihot < nhi; hihot++) { fprintf(stderr, " " ixFMT, ixhi_hot[hihot]); }
    fprintf(stderr, "\n");
    
    bool_t tried_empty = FALSE; /* Have we tried an empty clipbox? */
    
    for (uint32_t lohot = 0;  lohot < nlo; lohot++)
      { for (uint32_t hihot = 0;  hihot < nhi; hihot++)
          { ppv_index_t ixlo_r = ixlo_hot[lohot];
            ppv_index_t ixhi_r = ixhi_hot[hihot];
            if ((! tried_empty) || (ixlo_r <= ixhi_r))
              { ppv_size_t size_r = (ixlo_r > ixhi_r ? 0 : ixhi_r - ixlo_r + 1);
                process(ixlo_r,  size_r);
                tried_empty = tried_empty | (size_r == 0);
              }
          }
      }
  }

void kdtom_test_clip_core(kdtom_t *T, ppv_axis_t ax, ppv_index_t ixlo_ax, ppv_size_t size_ax)
  {
    bool_t debug = TRUE;
    #define dbtag "{kdtom_test_clip_core}: "

    fprintf(stderr, "--- testing {kdtom_clip_core} ax = %d ixlo = " ixFMT, ax, ixlo_ax);
    fprintf(stderr, " size = " szFMT " ---\n", size_ax);
   
    ppv_dim_t d = T->d;
    
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
      /* Requires {ix} to be inside {T.DK}, and returns {T.V[ix]}. */

    if (debug) { ixbox_print(stderr, d, dbtag "original box = [ ", T->ixlo, T->size, " ]\n"); }

    ppv_index_t ixlo_box[d]; /* {ixlo} value of clip box. */
    ppv_size_t size_box[d];  /* Size vector of clip box. */
    ppv_index_t ixlo_new[d]; /* Expected {T.h.ixlo} after clipping. */
    ppv_size_t size_new[d];  /* Expected {T.h.size} after clipping. */
    
    bool_t empty_box = (size_ax == 0);  /* Clipping box is empty. */
    if (empty_box)
      { ixlo_box[ax] = 0;
        size_box[ax] = 0;
      }
    else
      { /* Set up the clipping box on axis {ax}: */ 
        ixlo_box[ax] = ixlo_ax; 
        size_box[ax] = size_ax;
      }
    assert(T->d == d);  /* !!! DEBUG */
                
    /* Set clipping box index range along other axes as non-clipping: */
    for (ppv_axis_t k = 0; k < d;  k++)
      { if (k != ax)
          { ixlo_box[k] = (empty_box ? 0 : T->ixlo[k] - 1);
            size_box[k] = (empty_box ? 0 : T->size[k] + 2);
          }
      }

    if (debug) { ixbox_print(stderr, d, dbtag "clipping box = [ ", ixlo_box, size_box, " ]\n"); }

    /* Compute the expected clipped core box {ixlo_new,size_new}: */
    ixbox_intersect(d, T->ixlo, T->size, ixlo_box, size_box, ixlo_new, size_new); 
    if (debug) { ixbox_print(stderr, d, dbtag "expected box = [ ", ixlo_new, size_new, " ]\n"); }

    kdtom_t *S = kdtom_clip_core((kdtom_t *)T, ixlo_box, size_box);
    assert(T->d == d);  /* !!! DEBUG */
    if (debug) { ixbox_print(stderr, d, dbtag "clipped box =  [ ", S->ixlo, S->size, " ]\n"); }

    demand(S->d == T->d, "{kdtom_clip_core} bug: {S.d}");
    demand(S->maxsmp == T->maxsmp, "{kdtom_clip_core} bug: {S.maxsmp}");
    demand(S->fill == T->fill, "{kdtom_clip_core} bug: {S.fill}");
    demand(ixbox_is_contained(d, S->ixlo, S->size, ixlo_new, size_new), "{kdtom_clip_core} bug: domain");

    kdtom_test_get_sample(S, getsmp);
    if (debug) { fprintf(stderr, "\n"); }
    if (S != T) { free(S); }
    assert(T->d == d);  /* !!! DEBUG */

    return;
    
    auto ppv_sample_t getsmp(ppv_index_t ix[])
      { 
        assert(ixbox_has(d, ix, S->ixlo, S->size));
        return kdtom_get_sample(T, ix);
      }
    #undef dbtag
  }

void kdtom_test_translate(kdtom_t *T)
  {
    bool_t debug = TRUE;
    #define dbtag "{kdtom_test_translate}: "
    
    fprintf(stderr, "--- testing {kdtom_translate} ---\n");
    ppv_dim_t d = T->d;
    if (debug) { ixbox_print(stderr, d, dbtag "T core in = [ ", T->ixlo, T->size, " ]\n"); }
    
    bool_t empty = kdtom_has_empty_core((kdtom_t *)T);
    ppv_index_t dx[d];       /* Displacement to apply. */
    if (debug) { fprintf(stderr, dbtag "displ = [");  }
    ppv_index_t ixlo_old[d];
    for (ppv_axis_t k = 0; k < d;  k++) 
      { dx[k] = 400 -10*k;
        ixlo_old[k] = T->ixlo[k];
        if (debug) { fprintf(stderr, " " ixFMT, dx[k]); }
      }
    if (debug) { fprintf(stderr, " ]\n"); }

    kdtom_t *S = kdtom_clone(T);
    if (debug) { ixbox_print(stderr, d, dbtag "S core bf = [", S->ixlo, S->size, " ]\n"); }
    
    kdtom_translate((kdtom_t *)S, dx);
    if (debug) { ixbox_print(stderr, d, dbtag "S core af = [", S->ixlo, S->size, " ]\n"); }
    
    /* Check: */
    for (ppv_axis_t k = 0; k < d; k++) 
      { demand(S->ixlo[k] == (empty ? 0 : ixlo_old[k] + dx[k]), "{kdtom_translate} bad {T.ixlo}"); 
        demand(S->size[k] == (empty ? 0 : T->size[k]), "{kdtom_translate} bad {T.size}"); 
      }
    
    auto ppv_sample_t getsmp(ppv_index_t ix[]);
    kdtom_test_get_sample(S, getsmp);
    
    free(S);

    return;
    
    ppv_sample_t getsmp(ppv_index_t ix[])
      { ppv_index_t jx[d];
        for (ppv_axis_t k = 0; k < d; k++) { jx[k] = ix[k] - dx[k]; }
        assert(ixbox_has(T->d, jx, T->ixlo, T->size));
        return kdtom_get_sample(T, jx);
      }
    #undef dbtag
  }

ppv_sample_t *kdtom_test_pick_max_samples(int32_t nms)
  { bool_t debug = FALSE;
    #define dbtag "{kdtom_test_pick_max_samples}: "
    
    ppv_sample_t *ms = notnull(malloc(nms*sizeof(ppv_sample_t)), "no mem");
    for (uint32_t k = 0;  k < nms; k++)
      { ppv_nbits_t bps = (ppv_nbits_t)floor(32*((double)k)/(nms-1) + 0.5);
        if (bps == 0)
          { ms[k] = 0; }
        else
          { ppv_sample_t maxmaxsmp = ppv_max_sample(bps);
            ms[k] = (ppv_sample_t)(uint64_abrandom(maxmaxsmp/2+1, maxmaxsmp));
            assert(ppv_min_bps(ms[k]) == bps);
          }
        if (debug) { fprintf(stderr, dbtag " bps = %2d  ms[%2d] = " smpFMT "\n", bps, k, ms[k]); }
      }
    return ms;
    #undef dbtag
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
    bool_t debug = FALSE;
    #define dbtag "{kdtom_test_choose_array_size}: "
 
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
        if (debug) { fprintf(stderr, dbtag "array size = ["); }
        for (ppv_axis_t k = 0; k <d; k++)
          { size[k] = (ppv_size_t)floor(xsize*rsz[k]/rsz_avg);
            if (debug) { fprintf(stderr, " %lu", size[k]); }
            npos *= size[k];
          }
        if (debug) { fprintf(stderr, " ]\n"); }
      }
    if (debug) { fprintf(stderr, dbtag "array has %lu samples\n", npos); }
    return;
    #undef dbtag
  }

void kdtom_test_plot(kdtom_t *T)
  {
    demand(T->d == 2, "Can't plot this");
    fprintf(stderr, "Pretend that I am plotting...\n");    
    fprintf(stderr, "*plot*, *plat*, *plit*, *plet*, *plut*...\n");
    return;
  }

void kdtom_test_check_tree(kdtom_t *T)
  {
    bool_t debug = TRUE;
    #define dbtag "{kdtom_test_check_tree}: "

    ppv_dim_t d = T->d;
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");

    size_t rec_bytes = kdtom_bytesize(T, FALSE);
    size_t tot_bytes = kdtom_bytesize(T, TRUE);
    if (debug) 
      { fprintf(stderr, dbtag "byte sizes: root record = %lu", rec_bytes); 
        fprintf(stderr, " whole tree = %lu\n", tot_bytes);
      }
    
    ppv_sample_t maxsmp = T->maxsmp;
    demand(maxsmp <= ppv_MAX_SAMPLE_VAL, "invalid {maxsmp}");
    ppv_sample_t fill = T->fill;
    demand(fill <= maxsmp, "invalid {fill}");
    
    /* Either all axes have zero size, none has: */
    for (ppv_axis_t k = 0; k <d; k++)
      { demand((T->size[k] == 0) == (T->size[0] == 0), "inconsisten empty {size}"); }
    
    if (debug) { fprintf(stderr, dbtag "!! incomplete"); }
    return;
    #undef dbtag
  }

    
