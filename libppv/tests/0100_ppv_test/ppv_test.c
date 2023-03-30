/* Last edited on 2023-03-18 10:57:57 by stolfi */ 
/* Test of the PPV library. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <ix.h>
#include <ix_io.h>

#include <ppv_array.h>

#define posFMT ppv_pos_t_FMT
#define smpFMT ppv_sample_t_FMT
#define ixFMT  ppv_index_t_FMT
#define szFMT  ppv_size_t_FMT
#define stFMT  ppv_step_t_FMT
#define sctFMT ppv_sample_count_t_FMT

#define bug(FMT,...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "** " FMT "\n", __VA_ARGS__); \
    exit(1); \
  } while(0)

#define bug0(MSG) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fputs(MSG, stderr); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

void do_tests(ppv_dim_t d, ppv_sample_t maxsmp);

ppv_sample_t ix_to_sample(ppv_dim_t d, const ppv_index_t ix[], int32_t K, int32_t L, ppv_sample_t maxsmp);
  /* A pseudorandom sample value that depends on the index vector {ix[0..d-1]}
    and on the parameters {K} and {L}.  The sample will be in the range
    {0..maxsmp}. */

ppv_sample_t pos_to_sample(ppv_pos_t pos, int32_t K, int32_t L, ppv_sample_t maxsmp);
  /* A pseudorandom sample value that depends on the element position {pos}
    and on the parameters {K} and {L}.  The sample will be in the range
    {0..maxsmp}. */
    
void enum_by_hand(ix_index_op_t op, ppv_array_t *A);
  /* Calls {op(ix)} on all index vectors of the array {A}, in lex order.
    The {op} argument must be a procedure of type {ix_index_op_t}, that
    receives an index vector and returns a {bbol_t} ({TRUE} to stop the
    enumration). */

void test_sample_distr(ppv_array_t *A);
  /* Computes the histogram of the sample values in {A}. Prints
    the histogram to {stderr} (if not too large).  Prints
    statistics of the histogram.  
    
    Warning: requires memory proportional to {2^A->bps}, and is not
    useful unless the number of samples in {A} is many times that value.
    Will do nothing if {A->bps} is zero or way too large. */

void test_best_bpw(void);
void test_new_array(ppv_array_t *A, ppv_size_t *sz, ppv_sample_t maxsmp);
void test_packing(ppv_array_t *A);
void test_sample_pos(ppv_array_t *A);
void test_enum(ppv_array_t *A);
void test_assign(ppv_array_t *A, ppv_sample_t maxsmpB);
void test_sample_range(ppv_array_t *A);
void test_throw_noise(ppv_array_t *A, ppv_sample_count_t seed);
void test_throw_balls(ppv_array_t *A, ppv_sample_count_t seed);

void test_crop(ppv_array_t *A);
void test_subsample(ppv_array_t *A);
void test_reverse(ppv_array_t *A);
void test_replicate(ppv_array_t *A);
void test_swap_indices(ppv_array_t *A);
void test_flip_indices(ppv_array_t *A);
void test_slice(ppv_array_t *A);
void test_diagonal(ppv_array_t *A);
void test_chop(ppv_array_t *A);

void check_size(ppv_dim_t d, ppv_size_t *sza, ppv_size_t *szb);

int32_t main (int32_t argn, char **argv)
  {
    test_best_bpw();
    for (ppv_nbits_t bps = 0; bps <= ppv_MAX_BPS; bps++)
      { ppv_sample_t maxmaxsmp = ppv_max_sample(bps);
        int64_t m0 = int64_abrandom(0, maxmaxsmp);
        int64_t m1 = int64_abrandom(0, maxmaxsmp);
        ppv_sample_t maxsmp = (ppv_sample_t)imax(m0, m1);
        do_tests(2,  maxsmp);
      }
    do_tests(5,          65536);  
    do_tests(1,      (1<<25)-1);
    do_tests(6, ppv_MAX_SAMPLE_VAL);
    do_tests(0,             17);
    return 0;
  }

void test_best_bpw(void)
  {
    fprintf(stderr, "Testing {ppv_best_bpw}...\n");
    fprintf(stderr, "%5s %5s %3s %3s %6s\n", "bps", "bpw", "wps", "spw", "loss");
    fprintf(stderr, "%5s %5s %3s %3s %6s\n", "-----", "-----", "---", "---", "------");
    for (ppv_nbits_t bps = 0; bps <= ppv_MAX_BPS; bps++)
      { ppv_nbits_t bpw = ppv_best_bpw(bps);
        uint32_t wps; /* Words per sample. */
        uint32_t spw; /* Samples in {wps} words. */
        double waste; /* Fraction of memory wasted. */
        if (bps == 0)
          { wps = 0; spw = 0; waste = 0.0; }
        else
          { wps = (bps + bpw - 1)/bpw;
            spw = (wps*bpw)/bps;    
            waste = ((double)(wps*bpw - spw*bps))/((double)(wps*bpw));
          }
        fprintf(stderr, "%5u %5u %3u %3u %6.3f\n", bps, bpw, wps, spw, waste);
      }
    fprintf(stderr, "\n");
  }

ppv_sample_t ix_to_sample(ppv_dim_t d, const ppv_index_t ix[], int32_t K, int32_t L, ppv_sample_t maxsmp)
  { double s = 0;
    for (ppv_axis_t ax = 0; ax < d; ax++)
      { s += sin(K*ax + L)*((double)ix[ax]); }
    s = s - floor(s);
    ppv_sample_t smp = (ppv_sample_t)floor((maxsmp + 0.99999)*s);
    assert((smp >= 0) && (smp <= maxsmp));
    return smp;
  }

ppv_sample_t pos_to_sample(ppv_pos_t pos, int32_t K, int32_t L, ppv_sample_t maxsmp)
  { ppv_sample_t smp;
    if (pos == 0)
      { smp = maxsmp; }
    else
      { double ang = (double)(((int64_t)K) + ((int64_t)pos)*((int64_t)L));
        double s = sin(ang);
        s = s - floor(s);
        smp = (ppv_sample_t)floor((maxsmp + 0.99999)*s);
        assert((smp >= 0) && (smp <= maxsmp));
      }
    return smp;
  }

void do_tests(ppv_dim_t d, ppv_sample_t maxsmp)
  {
    fprintf(stderr, "--- {do-tests} d = %d maxsmp = " smpFMT " ---\n", d, maxsmp);

    ppv_dim_t dmax = 6;
    assert(d <= dmax);
    ppv_size_t szmax[6] = { 3, 10, 20, 30, 7, 2 };
    
    /* Create the array: */
    ppv_size_t sz[d];
    for (ppv_axis_t ax = 0; ax < d; ax++) { sz[ax] = szmax[ax]; }
    ppv_array_t *A = ppv_array_new(d, sz,  maxsmp);

    test_new_array(A, sz, maxsmp);
    test_packing(A);
    test_sample_pos(A);
    
    test_enum(A);
    test_assign(A, A->maxsmp);
    if (A->maxsmp < ppv_MAX_SAMPLE_VAL) 
      { ppv_sample_t maxsmpB = (ppv_sample_t)uint64_abrandom(A->maxsmp + 1, ppv_MAX_SAMPLE_VAL);
        test_assign(A, maxsmpB);
      }
    test_sample_range(A);
    test_throw_noise(A, 4615);
    test_throw_noise(A, 0);
    test_throw_balls(A, 4615);
    
    test_crop(A);
    test_subsample(A);
    test_reverse(A);
    test_replicate(A);
    test_swap_indices(A);
    test_flip_indices(A);
    test_slice(A);
    test_diagonal(A);
    test_chop(A);
    
    fprintf(stderr, "done.\n");
    return;
  }

void enum_by_hand(ix_index_op_t op, ppv_array_t *A)
  {
    ppv_dim_t dmax = 6;
    ppv_dim_t d = A->d;
    assert(d <=dmax);
    
    ppv_size_t szA[dmax];  /* Vector {A->size} extended with 1s to {dmax} entries. */
    for (ppv_axis_t ax = 0; ax < dmax; ax++) { szA[ax] = (ax < d ? A->size[ax] : 1); }
    ppv_index_t ixA[dmax];
    for (ixA[5] = 0; ixA[5] < szA[5]; ixA[5]++)
      for (ixA[4] = 0; ixA[4] < szA[4]; ixA[4]++)
        for (ixA[3] = 0; ixA[3] < szA[3]; ixA[3]++)
          for (ixA[2] = 0; ixA[2] < szA[2]; ixA[2]++)
            for (ixA[1] = 0; ixA[1] < szA[1]; ixA[1]++)
              for (ixA[0] = 0; ixA[0] < szA[0]; ixA[0]++)
                { bool_t stop = op(ixA);
                  if (stop) { return; }
                }
    return;
  }

void test_new_array(ppv_array_t *A, ppv_size_t *sz, ppv_sample_t maxsmp)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking num of axes...\n");

    fprintf(stderr, "Checking {ppv_array_new,ppv_descriptor_is_valid}...\n");
    bool_t verbose = TRUE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    if (! ppv_descriptor_is_valid(A, TRUE)) { bug0("{ppv_array_new}: returns invalid desc"); }
    if (A->base != 0) { bug("base = " posFMT, A->base); }
    if (A->maxsmp != maxsmp) { bug("maxsmp = " smpFMT, A->maxsmp); }
    /* Check steps for a vanilla array: */
    ppv_sample_count_t stp = 1;
    for (ppv_axis_t ax = 0; ax < d; ax++) 
      { if (A->step[ax] != stp) { bug("step[%d] = "  stFMT, ax, A->step[ax]); }
        stp *= sz[ax];
      }

    fprintf(stderr, "Checking {ppv_index_is_valid}...\n");
    ppv_index_t ixA[d];
    ppv_index_clear(d, ixA);
    if (! ppv_index_is_valid(ixA, A)) { bug0("bogus ppv_index_is_valid"); } 
    for (ppv_axis_t ax = 0; ax < d; ax++)
      { ixA[ax] = sz[ax]-1;
        if (! ppv_index_is_valid(ixA, A)) { bug0("bogus ppv_index_is_valid"); } 
        ixA[ax] = sz[ax];
        if (ppv_index_is_valid(ixA, A)) { bug0("bogus ppv_index_is_valid"); } 
        ixA[ax] = 0;
      }
  }

void test_packing(ppv_array_t *A)
  {
    fprintf(stderr, "Checking {ppv_{get,set}_sample_at_pos} on a few elems...\n");
    bool_t verbose = TRUE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t npos = 5; /* Number of positions to test */
    fprintf(stderr, "\n");
    for (ppv_pos_t pos = 0; pos < npos; pos++)
      { ppv_sample_t smp_set = pos_to_sample(pos, 417, 4615, A->maxsmp);
        if (verbose) { fprintf(stderr, "setting sample " posFMT " to %u\n", pos, smp_set); }
        ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp_set);
        if (verbose) { ppv_dump_storage(stderr, A->el, pos, A->bps, A->bpw, 2); }
        if (verbose) { fprintf(stderr, "\n"); }
        ppv_sample_t smp_get = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
        if (smp_get != smp_set) 
          { bug("set/get sample mismatch set = " smpFMT " get = " smpFMT, smp_set, smp_get); } 
      }
    fprintf(stderr, "\n");
  }

void test_sample_pos(ppv_array_t *A)
  {
    fprintf(stderr, "Checking {ppv_sample_pos,ppv_{get,set}_sample{,_at_pos}} on all elems...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    
    ppv_dim_t d = A->d;

    fprintf(stderr, "  filling array...\n");
    auto bool_t fill1(const ppv_index_t ixA[]);
      /* Checks {ppv_sample_pos}, sets/gets/compares each sample. */ 
    ppv_pos_t xpos = 0;
    enum_by_hand(fill1, A);
    
    fprintf(stderr, "  re-checking the value of all samples...\n");
    auto bool_t check1(const ppv_index_t ixA[]);
      /* Checks if samples set by {fill1} are still there. */ 
    enum_by_hand(check1, A);
    
    return;
    
    bool_t fill1(const ppv_index_t ixA[])
      { 
        if (! ppv_index_is_valid(ixA, A)) { bug0("not ppv_index_is_valid"); } 
        ppv_pos_t pos = ppv_sample_pos(A, ixA);
        if (pos != xpos) 
          { bug("{ppv_sample_pos}: mismatch pos = " posFMT " xpos = " posFMT, pos, xpos); } 

        ppv_sample_t smp1 = ix_to_sample(d, ixA, 3, +17, A->maxsmp);
        ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp1);
        ppv_sample_t xsmp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
        if (xsmp != smp1) 
          { bug("{set_sample_at_pos/get_sample_at_pos}: mismatch smp = " smpFMT " xsmp = " smpFMT, smp1, xsmp); } 
      
        ppv_sample_t smp2 = ix_to_sample(d, ixA, 7, +11, A->maxsmp);
        ppv_set_sample(A,ixA, smp2);
        ppv_sample_t ysmp = ppv_get_sample(A, ixA);
        if (ysmp != smp2) 
          { bug("{set_sample/get_sample}: mismatch smp = " smpFMT " ysmp = " smpFMT, smp2, ysmp); } 

        xpos++;
        return FALSE;
      }
    
    bool_t check1(const ppv_index_t ixA[])
      { 
        ppv_sample_t smp = ix_to_sample(d, ixA, 7, +11, A->maxsmp);
        ppv_pos_t pos = ppv_sample_pos(A, ixA);
        ppv_sample_t xsmp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
        if (xsmp != smp) 
          { bug("{get_sample_at_pos} mismatch on 2nd pass smp = " smpFMT " xsmp = " smpFMT, smp, xsmp); } 
        ppv_sample_t ysmp = ppv_get_sample(A, ixA);
        if (ysmp != smp) 
          { bug("{ge_sample}: mismatch smp = " smpFMT " ysmp = " smpFMT, smp, ysmp); } 
        return FALSE;
      }
  }

void test_enum(ppv_array_t *A)
  {
    fprintf(stderr, "!! NOT checking {ppv_enum} yet!\n");
  }
  
void test_sample_range(ppv_array_t *A)
  {
    fprintf(stderr, "!! NOT checking {ppv_sample_range} yet!\n");
  }

void test_throw_noise(ppv_array_t *A, ppv_sample_count_t seed)
  { 
    fprintf(stderr, "Checking {ppv_throw_noise}");
    fprintf(stderr, " seed = " ppv_sample_count_t_FMT, seed);
    fprintf(stderr, " A.maxsmp = " smpFMT " ...\n", A->maxsmp);

    srandom((uint32_t)seed);
    ppv_throw_noise(A);
    
    test_sample_distr(A);
    return;
        
  }
    
void test_throw_balls(ppv_array_t *A, ppv_sample_count_t seed)
  { 
    fprintf(stderr, "Checking {ppv_throw_balls}");
    fprintf(stderr, " seed = " ppv_sample_count_t_FMT, seed);
    fprintf(stderr, " A.maxsmp = " smpFMT " ...\n", A->maxsmp);

    srandom((uint32_t)seed);
    if (A->maxsmp >= 1) { ppv_throw_balls(A); }
    
    test_sample_distr(A);
    return;
  }

void test_sample_distr(ppv_array_t *A)
  { 
    ppv_sample_t smp_min, smp_max;
    ppv_sample_range(A, &smp_min, &smp_max);
    fprintf(stderr, "sample range = { " smpFMT, smp_min);
    fprintf(stderr, ".. " smpFMT " }\n", smp_max);
    
    ppv_sample_count_t npos = ppv_sample_count(A, TRUE); /* Voxel count, incl. replics. */
    fprintf(stderr, "sample count = " ppv_sample_count_t_FMT "\n", npos);
    ppv_sample_count_t nposmin = 2; /* Min voxel count. */
    if (npos < nposmin) 
      { fprintf(stderr, "!! no samples -- histogram not computed.\n");
        return;
      }

    uint32_t nvals = (uint32_t)(smp_max - smp_min + 1);
    assert(nvals >= 1);

    bool_t prhist = (nvals <= 256);

    int64_t maxnvals = imax(1, imin(1 << 18, npos/2));
    if (nvals > maxnvals)
      { /* No use mapping samples to a smaller range, so: */
        fprintf(stderr, "!! sample value range too big -- histogram not computed\n");
        return;
      }

    uint64_t hist[nvals]; 
    for (int32_t v = 0; v < nvals; v++) { hist[v] = 0; }

    auto bool_t gather_hist(const ppv_index_t ix[]);
    enum_by_hand(gather_hist, A);
    if (prhist) 
      { fprintf(stderr, "%5s %12s\n", "val", "count");
        fprintf(stderr, "%5s %12s\n", "-----", "------------");
      }

    double hp = 1.0/((double)nvals); /* Prob of a sample falling into a bin. */
    uint64_t hmin = ppv_MAX_SAMPLES;
    uint64_t hmax = 0;
    double sum_h = 0.0;
    for (int32_t v = 0; v < nvals; v++)
      { uint64_t hv = hist[v];
        if (prhist) { fprintf(stderr, "%5u %12lu\n", (uint32_t)(v + smp_min), hv); }
        if (hv > hmax) { hmax = hv;  }
        if (hv < hmin) { hmin = hv;  }
        sum_h += (double)(hv);
      }
    if (prhist) { fprintf(stderr, "\n"); }
    double havg = sum_h/(double)nvals;   /* Actual average bin count. */
    double sum_d2 = 0;
    for (int32_t v = 0; v < nvals; v++)
      { uint64_t hv = hist[v];
        double dv = ((double)hv) - havg;
        sum_d2 += dv*dv;
      }
    double hdev = sqrt(sum_d2/((double)nvals-1));       /* Actual dev of bin counts. */
    
    /* ??? Should compute chi-square statistic. ??? */
    fprintf(stderr, "value count statistics:\n");
    fprintf(stderr, "  max = %lu min = %lu\n", hmin, hmax);
    
    double havg_exp = hp*((double)npos); /* Expected average bin count. */
    assert(havg_exp > 0.0);
    double havg_err = (havg - havg_exp)/havg_exp; /* Error in average. */
    fprintf(stderr, "  avg = %.3f (exp = %.3f, err = %.1f%%)\n", havg, havg_exp, 100*havg_err);
    
    if (nvals >= 2)
      { double hdev_exp = sqrt(hp*(1-hp)*((double)npos)); /* Expected dev of bin counts. */
        assert(hdev_exp > 0.0);
        double hdev_err = (hdev - hdev_exp)/hdev_exp; /* Error in deviation. */
        fprintf(stderr, "  dev = %.3f (exp = %.3f, err = %.1f%%)\n", hdev, hdev_exp, 100*hdev_err);
      }
    return;

    /* Internal procs */

    bool_t gather_hist(const ppv_index_t ix[])
      { ppv_sample_t smp = ppv_get_sample(A, ix);
        assert ((smp_min <= smp) && (smp <= smp_max));
        hist[(int32_t)(smp - smp_min)]++;
        return FALSE;
      }
  }

void test_assign(ppv_array_t *A, ppv_sample_t maxsmpB)
  {
    fprintf(stderr, "Checking {ppv_assign} A.maxsmp = " smpFMT, A->maxsmp);
    fprintf(stderr, " B.maxsmp = " smpFMT " ...\n", maxsmpB);
    
    ppv_dim_t d = A->d;

    /* Check if {A} has replicated elements: */
    bool_t repl = FALSE; 
    for (ppv_axis_t ax = 0; ax < d; ax++)
      { if ((A->size[ax] > 1) && (A->step[ax] == 0))
          { repl = TRUE; }
      }
    if (repl)
      { fprintf(stderr, "test skipped because {A} has replicated elements\n");
        return;
      }
    
    /* Create array {B} with same indices as {A}: */
    ppv_array_t *B = ppv_array_new(d, A->size, maxsmpB);
    fprintf(stderr, "A.el = %016lx B.el = %016lx\n", (uint64_t)A->el, (uint64_t)B->el);
    
    /* Fill {A} and {B} with different garbage: */
    auto bool_t fillAB(const ppv_index_t ixA[]);
    enum_by_hand(fillAB, A);
    
    /* Copy {A} into {B}: */
    ppv_array_assign(B, A);
    
    /* Check whether the assignment succeeded: */
    auto bool_t checkAB(const ppv_index_t ixA[]);
    enum_by_hand(checkAB, A);
    
    return;
    
    bool_t fillAB(const ppv_index_t ixA[])
      { ppv_sample_t smpA = ix_to_sample(d, ixA, 3, +17, A->maxsmp);
        ppv_sample_t smpB = ix_to_sample(d, ixA, 5, -23, B->maxsmp);
        ppv_set_sample(A, ixA, smpA);
        ppv_set_sample(B, ixA, smpB);
        return FALSE;
      }
    
    bool_t checkAB(const ppv_index_t ixA[])
      { ppv_sample_t smpA_get = ppv_get_sample(A, ixA);
        ppv_sample_t smpB_get = ppv_get_sample(B, ixA);
        if (smpB_get != smpA_get)
          { ix_print_indices(stderr, "ix = (", d, ixA, 0, " ", ")\n");
            bug("{ppv_assign}: {A} and {B} differ after assignment: B now = " smpFMT " assigned = " smpFMT, smpB_get, smpA_get);
          }
        ppv_sample_t smpA_set = ix_to_sample(d, ixA, 3, +17, A->maxsmp);
        if (smpA_get != smpA_set)
          { ix_print_indices(stderr, "ix = (", d, ixA, 0, " ", ")\n");
            bug("{ppv_assign}: {A} changed set = " smpFMT " get = " smpFMT "\n", smpA_set, smpA_get);
          }
        return FALSE;
      }
  }

void test_crop(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_crop}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_step_t sh[d]; /* Shift of cropped array along each axis. */
    ppv_size_t sz[d]; /* Size of cropped array along each axis. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {sh} and {sz} for full array: */
        for (ppv_axis_t ax = 0; ax < d; ax++) { sh[ax] = 0; } 
        memcpy(sz, A->size, d*sizeof(ppv_size_t)); 
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Perform first cropping: */
        ppv_axis_t ax1 = (ppv_axis_t)(trial % d);
        ppv_size_t skip1 = ((trial + trial/2) & 3)*sz[ax1]/8;
        ppv_size_t keep1 = ((trial + trial/4) & 3)*(sz[ax1] - skip1)/8;
        sh[ax1] += skip1; sz[ax1] = keep1;
        ppv_crop(B, ax1, skip1, keep1);
        if (verbose) 
          { fprintf(stderr, "ppv_crop(%d, " szFMT  szFMT ")", ax1, skip1, keep1);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_crop}: returns an invalid desc"); }

        /* Perform second cropping (possibly on same axis): */ 
        ppv_axis_t axa2 = (ppv_axis_t)((trial / d) % d);
        ppv_size_t skip2 = ((trial + trial/2) & 3)*sz[axa2]/8;
        ppv_size_t keep2 = ((trial + trial/5) & 3)*(sz[axa2] - skip2)/8;
        sh[axa2] += skip2; sz[axa2] = keep2;
        ppv_crop(B, axa2, skip2, keep2);
        if (verbose) 
          { fprintf(stderr, "ppv_crop(%d, " szFMT ", " szFMT ")", axa2, skip2, keep2);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_crop}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if ((keep1 == 0) || (keep2 == 0))
          { if (B->el != NULL) { bug0("{ppv_crop}: empty array with non-null storage"); } }
        else
          { if (B->el != A->el) { bug0("{ppv_crop}: garbled storage area"); } }
        /* Check whether the size of {B} is correct: */
        check_size(d, B->size, sz);
        /* Now check coincidence of the two arrays: */
        
        auto bool_t check1(const ppv_index_t ixB[]);
        enum_by_hand(check1, B);

        continue;

        /* Internal procedures: */

        bool_t check1(const ppv_index_t ixB[])
          { 
            ppv_index_t ixA[d];
            ppv_index_assign(d, ixA, ixB); ppv_index_shift(d,ixA, sh);
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_crop}: position error posB = " posFMT " posA = " posFMT, posB, posA); }
            return FALSE;
          }

      }
        
    return;
  }

void test_subsample(ppv_array_t *A) 
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_subsample}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t st[d]; /* Subsampling step along each axis. */
    ppv_size_t sz[d]; /* Size of subsampled array along each axis. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {st} and {sz} for full array: */
        for (ppv_axis_t ax = 0; ax < d; ax++) { st[ax] = 1; } 
        memcpy(sz, A->size, d*sizeof(ppv_size_t)); 
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Perform first subsampling: */
        ppv_axis_t ax1 = (ppv_axis_t)(trial % d);
        ppv_size_t step1 = 1 + (trial & 3)*(sz[ax1] - 1)/8;
        ppv_size_t keep1 = (sz[ax1] + step1 - 1)/step1;
        st[ax1] *= step1; sz[ax1] = keep1;
        ppv_subsample(B, ax1, step1);
        if (verbose) 
          { fprintf(stderr, "ppv_subsample(%d, " szFMT ")", ax1, step1);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_subsample}: returns an invalid desc"); }

        /* Perform second subsampling (possibly on same axis): */ 
        ppv_axis_t axa2 = (ppv_axis_t)((trial / d) % d);
        ppv_size_t step2 = 1 + (trial/2 & 3)*(sz[axa2] - 1)/8;
        ppv_size_t keep2 = (sz[axa2] + step2 - 1)/step2;
        st[axa2] *= step2; sz[axa2] = keep2;
        ppv_subsample(B, axa2, step2);
        if (verbose) 
          { fprintf(stderr, "ppv_subsample(%d, " szFMT ")", axa2, step2);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_subsample}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if (B->el != A->el) { bug0("{ppv_subsample}: garbled storage area"); }
        /* Check whether the size of {B} is correct: */
        check_size(d, B->size, sz); 
        /* Now check coincidence of the two arrays: */
        
        auto bool_t check1(const ppv_index_t ixB[]);
        enum_by_hand(check1, B);

        return;
        
        bool_t check1(const ppv_index_t ixB[])       
          { ppv_index_t ixA[d];
            for (ppv_axis_t ax = 0; ax < d; ax++) { ixA[ax] = ixB[ax]*st[ax]; }
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_subsample}: position error posB = " posFMT " posA = " posFMT, posB, posA); } 
            return FALSE;
          }
      }
  }

void test_swap_indices(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    if (d < 2)
      { fprintf(stderr, "NOT checking {ppv_swap_indices}: {d} too small\n");
        return;
      }
    fprintf(stderr, "Checking {ppv_swap_indices}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t tr[d]; /* Axis {ax} of {B} is axis {tr[ax]} of {A}. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {tr} with identity permutation: */
        for (ppv_axis_t ax = 0; ax < d; ax++) { tr[ax] = ax; }
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        for (int32_t pass = 0; pass < 2; pass++)
          { /* Perform a transposition: */
            ppv_axis_t axa1 = (ppv_axis_t)((trial/(1+pass)) % d);
            ppv_axis_t axb1 = (ppv_axis_t)((trial/(2+pass)) % d);
            int32_t d1 = (axa1 < axb1 ? axb1 - axa1 : axa1 - axb1);
            int32_t m1 = d - (axa1 > axb1 ? axa1 : axb1);
            ppv_dim_t n1 = (ppv_dim_t)((trial/(3+pass)) % (d1 == 0 ? m1 : (d1 < m1 ? d1 : m1)));
            assert(axa1 + n1 <= d);
            assert(axb1 + n1 <= d);
            assert((axa1 == axb1) || (axa1 + n1 <= axb1) || (axb1 + n1 <= axa1));
            for (int32_t k = 0; k < n1; k++)
              { ppv_size_t tmp = tr[axa1+k]; tr[axa1+k] = tr[axb1+k]; tr[axb1+k] = tmp; } 
            ppv_swap_indices(B, axa1, axb1, n1);
            if (verbose) 
              { fprintf(stderr, "ppv_swap_indices(%d, %d, %d)", axa1, axb1, n1);
                ppv_print_descriptor(stderr, " = { ", B, " }\n");
              } 
            /* Check validity of descriptor: */
            if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_swap_indices}: returns an invalid desc"); }
          }

        /* Check consistency of storage area: */
        if (B->el != A->el) { bug0("{ppv_swap_indices}: garbled storage area"); }
        for (ppv_axis_t ax = 0; ax < d; ax++)
          { /* Check whether the size of {B} is correct: */
            if (B->size[ax] != A->size[tr[ax]]) 
              { bug("size[%d] = " szFMT, ax, B->size[ax]); }
            /* Check whether the increments of {B} are correct: */
            if (B->step[ax] != A->step[tr[ax]]) 
              { bug("step[%d] = " stFMT, ax, B->step[ax]); }
          }
        /* Now check coincidence of the two arrays: */
        auto bool_t check1(const ppv_index_t ixB[]);
        enum_by_hand(check1, A);
        
        return;
        
        bool_t check1(const ppv_index_t ixA[])       
          { ppv_index_t ixB[d];
            for (ppv_axis_t ax = 0; ax < d; ax++) { ixB[ax] = ixA[tr[ax]]; }
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_swap_indices}: pos error posB = " posFMT " posA = " posFMT, posB, posA); } 
            return FALSE;
          }
      }
  }

void test_flip_indices(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_flip_indices}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t tr[d]; /* Axis {ax} of {B} is axis {tr[ax]} of {A}. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {tr} with identity permutation: */
        for (ppv_axis_t ax = 0; ax < d; ax++) { tr[ax] = ax; }
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Perform first transposition: */
        ppv_axis_t axa1 = (ppv_axis_t)(trial % d);
        ppv_axis_t axb1 = (ppv_axis_t)(axa1 + (trial/2 % (d-axa1)));
        { ppv_axis_t axi=axa1, axj=axb1;
          while (axi < axj) 
            { ppv_size_t tmp = tr[axi]; tr[axi] = tr[axj]; tr[axj] = tmp; axi++; axj--; }
        }
        ppv_flip_indices(B, axa1, axb1);
        if (verbose) 
          { fprintf(stderr, "ppv_flip_indices(%d, %d)", axa1, axb1);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_flip_indices}: returns an invalid desc"); }

        /* Perform second transposition (possibly on same axis): */ 
        ppv_axis_t axa2 = (ppv_axis_t)((trial/d) % d);
        ppv_axis_t axb2 = (ppv_axis_t)(axa2 + ((trial/d/2) % (d-axa2)));
        { ppv_axis_t axi=axa2, axj=axb2;
          while (axi < axj) 
            { ppv_size_t tmp = tr[axi]; tr[axi] = tr[axj]; tr[axj] = tmp; axi++; axj--; }
        }
        ppv_flip_indices(B, axa2, axb2);
        if (verbose) 
          { fprintf(stderr, "ppv_flip_indices(%d, %d)", axa2, axb2);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_flip_indices}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if (B->el != A->el) { bug0("{ppv_flip_indices}: garbled storage area"); }
        for (ppv_axis_t ax = 0; ax < d; ax++)
          { /* Check whether the size of {B} is correct: */
            if (B->size[ax] != A->size[tr[ax]]) 
              { bug("size[%d] = " szFMT, ax, B->size[ax]); }
            /* Check whether the increments of {B} are correct: */
            if (B->step[ax] != A->step[tr[ax]]) 
              { bug("step[%d] = " stFMT, ax, B->step[ax]); }
          }
        /* Now check coincidence of the two arrays: */
        
        auto bool_t check1(const ppv_index_t ixA[]);
        enum_by_hand(check1, A);

        return;

        bool_t check1(const ppv_index_t ixA[])       
          { 
            ppv_index_t ixB[d];
            for (ppv_axis_t ax = 0; ax < d; ax++) { ixB[ax] = ixA[tr[ax]]; }
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_flip_indices}: pos error posB = " posFMT " posA = " posFMT, posB, posA); } 
            return FALSE;
          }
      }
  }

void test_reverse(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_reverse}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    bool_t fp[d]; /* Tells whether each axis was flipped or not. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {fp} for unflipped array: */
        for (ppv_axis_t ax = 0; ax < d; ax++) { { fp[ax] = FALSE; } }
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Perform first flip: */
        ppv_axis_t ax1 = (ppv_axis_t)(trial % d);
        fp[ax1] = ! fp[ax1];
        ppv_reverse(B, ax1);
        if (verbose) 
          { fprintf(stderr, "ppv_reverse(%d)", ax1);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_reverse}: returns an invalid desc"); }

        /* Perform second flip (possibly on same axis): */ 
        ppv_axis_t ax2 = (ppv_axis_t)((trial / d) % d);
        fp[ax2] = ! fp[ax2];
        ppv_reverse(B, ax2);
        if (verbose) 
          { fprintf(stderr, "ppv_reverse(%d)", ax2);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_reverse}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if (B->el != A->el) { bug0("{ppv_reverse}: garbled storage area"); }
        /* Check whether the size of {B} is correct: */
        check_size(d, B->size, A->size);
        /* Check whether the increments of {B} are correct: */
        { for (ppv_axis_t ax = 0; ax < d; ax++) 
            if (B->step[ax] != A->step[ax]*(fp[ax]?-1:+1)) 
              { bug("step[%d] =" stFMT, ax, B->step[ax]); }
        }
        /* Now check coincidence of the two arrays: */
        auto bool_t check1(const ppv_index_t ixB[]);
        enum_by_hand(check1, A);
        
        return;
        
        bool_t check1(const ppv_index_t ixB[])       
          { ppv_index_t ixA[d];
            for (ppv_axis_t ax = 0; ax < d; ax++) 
              { ixA[ax] = (fp[ax] ? A->size[ax] - 1 - ixB[ax] : ixB[ax]); }
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_reverse}: position error posB = " posFMT " posA = " posFMT, posB, posA); }
            return FALSE;
          }
                    
      }
  }

void test_slice(ppv_array_t *A)
  {
    fprintf(stderr, "!! NOT checking {ppv_slice}...\n");
  }

void test_diagonal(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    if (d < 2)
      { fprintf(stderr, "NOT checking {ppv_diagonal}: {d} too small\n");
        return;
      }
    fprintf(stderr, "Checking {ppv_diagonal}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t sz[d]; /* Size of diagonalized array along each axis. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        memcpy(sz, A->size, d*sizeof(ppv_size_t)); 
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Choose any two distinct axes {axa1,axb1}: */
        ppv_axis_t axa1 = (ppv_axis_t)(trial % d);
        ppv_axis_t axb1 = (ppv_axis_t)(trial/d % (d-1));
        if (axb1 >= axa1) { axb1++; }
        /* Make sure {sz[axa1] <= sz[axb1]}: */
        if (sz[axa1] > sz[axb1]) { ppv_axis_t t = axa1; axa1 = axb1; axb1 = t; }
        if (sz[axa1] > 0)
          { 
            /* Take diagonal slice of array: */
            ppv_diagonal(B, axa1, axb1);
            if (verbose) 
              { fprintf(stderr, "ppv_diagonal(%d, %d)", axa1, axb1);
                ppv_print_descriptor(stderr, " = { ", B, " }\n");
              } 
            /* Compute expected sizes of {1}: */
            sz[axb1] = sz[axb1] - (sz[axa1] - 1);
            /* Check validity of descriptor: */
            if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_diagonal}: returns an invalid desc"); }

            /* Check consistency of storage area: */
            if (B->el != A->el) { bug0("{ppv_diagonal}: garbled storage area"); }
            /* Check whether the size of {1} is correct: */
            check_size(d, B->size, sz);
            /* Now check coincidence of the two arrays: */
            auto bool_t check1(const ppv_index_t ixB[]);
            enum_by_hand(check1, B);

            return;

            bool_t check1(const ppv_index_t ixB[])       
              { ppv_index_t ixA[d];
                ppv_index_assign(d, ixA, ixB);
                ixA[axb1] += ixA[axa1];
                ppv_pos_t posB = ppv_sample_pos(B, ixB);
                ppv_pos_t posA = ppv_sample_pos(A, ixA);
                if (posB != posA)
                  { bug("{ppv_diagonal}: position error posB = " posFMT " posA = " posFMT, posB, posA); }
                return FALSE;
              }
          }
      }
  }

void test_chop(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_chop}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t sz[d]; /* Size of chopped array along each axis. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        memcpy(sz, A->size, d*sizeof(ppv_size_t)); 
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Choose an axis to chop: */
        ppv_axis_t axa1 = (ppv_axis_t)(trial % d);
        /* Choose a nonzero chunk size: */
        ppv_size_t chunksz1 = 1 + (trial/3 & 3)*(B->size[axa1] - 1)/8;
        /* Ensure that {B->size[axa1]} is a multiple of {chunksz1}: */
        if ((B->size[axa1] % chunksz1) != 0)
          { if (chunksz1 > B->size[axa1]) 
              { chunksz1 = (B->size[axa1] + 1)/2; }
            ppv_crop(B, axa1, 0, (B->size[axa1]/chunksz1)*chunksz1);
          }
        /* Chop array: */
        ppv_array_t *C = ppv_chop(B, axa1, chunksz1);
        assert(C->d == B->d + 1);
        ppv_axis_t axb1 = (ppv_axis_t)(C->d - 1); /* Chunk index axis. */
        /* Compute expected sizes of {C}: */
        sz[axb1] = sz[axa1]/chunksz1; sz[axa1] = chunksz1;
        if (verbose) 
          { fprintf(stderr, "ppv_chop(%d, " szFMT ", %d)", axa1, chunksz1, axb1);
            ppv_print_descriptor(stderr, " = { ", C, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(C, TRUE)) { bug0("{ppv_chop}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if (C->el != A->el) { bug0("{ppv_chop}: garbled storage area"); }
        /* Check whether the size of {1} is correct: */
        check_size(d, C->size, sz);
        /* Now check coincidence of the two arrays: */
        auto bool_t check1(const ppv_index_t ixA[]);
        enum_by_hand(check1, A);

        return;

        bool_t check1(const ppv_index_t ixA[])       
          { ppv_index_t ixC[d];
            ppv_index_assign(d, ixC, ixA);
            ixC[axa1] = ixA[axa1] % chunksz1;
            ixC[axb1] = ixA[axa1] / chunksz1;
            ppv_pos_t posC = ppv_sample_pos(C, ixC);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posC != posA)
              { bug("{ppv_chop}: pos error posC = " posFMT " posA = " posFMT, posC, posA); } 
            return FALSE;
          }
      }
  }

void test_replicate(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    fprintf(stderr, "Checking {ppv_replicate}...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    int32_t ntrials = 2*d*d;
    ppv_size_t sz[d]; /* Size of replicated array along each axis. */
    for (int32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        memcpy(sz, A->size, d*sizeof(ppv_size_t)); 
        /* Start with the standard array: */
        ppv_array_t *B = ppv_array_clone(A);

        /* Choose axis for first replication: */
        ppv_axis_t ax1 = (ppv_axis_t)(trial % d);
        /* Choose a positive replication factor: */
        ppv_size_t rep1 = 1 + (trial/3 & 3)*5;
        /* Chop and replicate along axis {ax1}: */
        ppv_size_t skip1 = B->size[ax1]/2;
        ppv_crop(B, ax1, skip1, 1);
        ppv_replicate(B, ax1, rep1);
        sz[ax1] = rep1;
        if (verbose) 
          { fprintf(stderr, "ppv_crop(%d," szFMT ",1)+ppv_replicate(%d," szFMT ")", ax1, skip1, ax1, rep1);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_replicate}: returns an invalid desc"); }

        /* Choose axis for second replication (may be {ax1}): */
        ppv_axis_t ax2 = (ppv_axis_t)((trial/d) % d);
        /* Choose a replication factor: */
        ppv_size_t rep2 = 1 + (trial/8 & 3)*5;
        /* Chop and replicate along axis {ax2}: */
        ppv_size_t skip2 = B->size[ax2]/2;
        ppv_crop(B, ax2, skip2, 1);
        ppv_replicate(B, ax2, rep2);
        sz[ax2] = rep2;
        if (verbose) 
          { fprintf(stderr, "ppv_crop(%d," szFMT ",1)+ppv_replicate(%d," szFMT ")", ax2, skip2, ax2, rep2);
            ppv_print_descriptor(stderr, " = { ", B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ppv_descriptor_is_valid(B, TRUE)) { bug0("{ppv_replicate}: returns an invalid desc"); }

        /* Check consistency of storage area: */
        if (B->el != A->el) { bug0("{ppv_replicate}: garbled storage area"); }
        /* Check whether the size of {1} is correct: */
        check_size(d, B->size, sz);
        /* Now check coincidence of the two arrays: */
        auto bool_t check1(const ppv_index_t ixB[]);
        enum_by_hand(check1, B);

        return;

        bool_t check1(const ppv_index_t ixB[])       
          { ppv_index_t ixA[d];
            ppv_index_assign(d, ixA, ixB);
            ixA[ax1] = skip1;
            /* The second crop changes {base} only if {ax1 != ax2}: */
            if (ax2 != ax1) { ixA[ax2] = skip2; }
            ppv_pos_t posB = ppv_sample_pos(B, ixB);
            ppv_pos_t posA = ppv_sample_pos(A, ixA);
            if (posB != posA)
              { bug("{ppv_replicate}: position error posB = " posFMT " posA = " posFMT, posB, posA); }
            return FALSE;
          }
      }
  }

void check_size(ppv_dim_t d, ppv_size_t *sza, ppv_size_t *szb)
  {
    for (ppv_axis_t ax = 0; ax < d; ax++)
      { if (sza[ax] != szb[ax]) { bug("size[%d] = " szFMT, ax, sza[0]); } }
  }
