/* See {kdtom_test.h}. */
/* Last edited on 2021-07-12 22:27:55 by jstolfi */

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

ppv_array_t *kdtom_test_make_array(ppv_dim_t d, ppv_size_t sz[], ppv_sample_t maxsmp)
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");

    ppv_array_t *A = ppv_array_new(d, sz, maxsmp);
    
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

void kdtom_test_choose_array_size(ppv_dim_t d, ppv_size_t sz[])
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
          { sz[k] = (ppv_size_t)floor(xsize*rsz[k]/rsz_avg);
            fprintf(stderr, " %lu", sz[k]);
            npos *= sz[k];
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
