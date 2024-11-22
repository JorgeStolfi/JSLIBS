/* Last edited on 2023-03-18 10:57:50 by stolfi */ 
/* Test of the I/O functions from the PPV library. */

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <bool.h>

#include <ppv_array.h>
#include <ppv_array_write.h>
#include <ppv_array_read.h>

#define posFMT ppv_pos_t_FMT
#define sizeFMT ppv_size_t_FMT
#define smpFMT ppv_sample_t_FMT

#define dMAX (7)

#define bug(nerrP,FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    (*nerrP)++; \
  } while(0)

ppv_array_t *make_test_array(ppv_dim_t d, ppv_nbits_t bps);
void test_fill_array(ppv_array_t *A);

void check_array_equality(ppv_array_t *A, ppv_array_t *B, int64_t *nerrP);
  /* Compares arrays {A} and {B}, calls {bug} if discrepant.
    Increments {*nerrP} on errors. */

void do_test_write(ppv_array_t *A, char *fname);
  /* Writes the array {A} to file {fname}. */
  
void do_test_read(ppv_array_t *A, char *fname, int64_t *nerrP);
  /* Reads an array {B} from file {fname}, which is could be either a
    2005 format or a 2021 format file. If {A} is not {NULL}, compares it
    to {A} with {check_array_equality(A,B,nerrP)}. */

int32_t main (int32_t argn, char **argv)
  {
    int64_t nerr = 0; /* Number of errors. */
    
    for (ppv_dim_t d = 0; d <= dMAX; d = (ppv_dim_t)(3*d + 2)/2)
      { ppv_nbits_t bps_lo = (d == 0 ? 3 : 0);
        ppv_nbits_t bps_hi = (d == 0 ? 3 : 32);
        
        for (ppv_nbits_t bps = bps_lo; bps <= bps_hi; bps++)
          { ppv_array_t *A = make_test_array(d, bps);
            char *fname = jsprintf("out/test_d%d_bps%02d.ppv", d, bps);
            
            /* Write a file in the current format: */
            do_test_write(A, fname);

            /* Read it back: */
            do_test_read(A, fname, &nerr);

            free(fname);
          }
      }

    /* Read a version in the old format: */
    do_test_read(NULL, "in/test_2005.ppv", &nerr);

    demand(nerr == 0, "there were errors");
    fprintf(stderr, "done.\n");
    return 0;
  }

void do_test_write(ppv_array_t *A, char *fname)
  {
    fprintf(stderr, "Checking {ppv_array_write_file} fname = %s...\n", fname);
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    FILE *wr = open_write(fname, TRUE);
    ppv_array_write_file(wr, A, TRUE);
    ppv_array_write_file(wr, A, FALSE);
    fclose(wr);
  }
  
void do_test_read(ppv_array_t *A, char *fname, int64_t *nerrP)
  {
    fprintf(stderr, "Checking {ppv_array_read_file} (same array)...\n");
    ppv_array_t *B;
    FILE *rd = open_read(fname, TRUE);
    for (int32_t i = 0; i < 2; i++)
      { B = ppv_array_read_file(rd);
        if (! ppv_descriptor_is_valid(B, TRUE)) 
          { bug(nerrP,"{ppv_array_read_file} returns an invalid desc"); }
        if (A != NULL) { check_array_equality(A, B, nerrP); }
        free(B->el);
      }
    fclose(rd);
  }
void check_array_equality(ppv_array_t *A, ppv_array_t *B, int64_t *nerrP)
  {
    ppv_dim_t d = A->d;
    demand(B->d == d, "axis counts don't match");
    
    /* Check array size: */
    for (ppv_axis_t i = 0; i < d; i++)
      { if (B->size[i] != A->size[i]) 
          { bug(nerrP, "size[%d] = " sizeFMT "", i, B->size[i]); }
      }
    
    /* Assume the 2021 format, with {B.maxsmp}: */
    if (B->maxsmp != A->maxsmp) { bug(nerrP, "wrong {.maxsmp}"); }
    if (B->bps != A->bps) { bug(nerrP, "wrong {.bps}"); }

    /* Check contents: */
    auto bool_t comp(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(comp, FALSE, A, B, NULL);
    return;
    
    bool_t comp(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t smpA = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        ppv_sample_t smpB = ppv_get_sample_at_pos(B->el, B->bps, B->bpw, pB);
        if (smpA != smpB)
          { bug(nerrP, "contents differ at positions posA = " posFMT " posB = " posFMT "", pA, pB); } 
        return FALSE;
      }
  }

ppv_array_t *make_test_array(ppv_dim_t d, ppv_nbits_t bps)
  {
    bool_t verbose = TRUE;

    fprintf(stderr, "creating test array d = %d bps =%d ...\n", d, bps);
    
    ppv_size_t sz[d];
    ppv_sample_count_t npos = 5000;
    ppv_choose_test_size(d, npos, sz);

    ppv_sample_t maxsmp;
    if (bps == 0)
      { maxsmp = 0; }
    else
      { ppv_sample_t maxmaxsmp = ppv_max_sample(bps);
        uint64_t m0 = uint64_abrandom(1, maxmaxsmp);
        uint64_t m1 = uint64_abrandom(1, maxmaxsmp);
        maxsmp = (ppv_sample_t)imax(m0,m1);
      }
    fprintf(stderr, "maxsmp = " smpFMT "\n", maxsmp);

    ppv_array_t *A = ppv_array_new(d, sz,  maxsmp);

    if (verbose) 
      { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 

    assert(ppv_descriptor_is_valid(A, TRUE));
      
    test_fill_array(A);
    return A;
  }

void test_fill_array(ppv_array_t *A)
  {
    /* Fill elements: */
    ppv_sample_t val = 4615;
    auto bool_t init(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(init, FALSE, A, NULL, NULL);
    return;
    
    bool_t init(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t smp = (ppv_sample_t)(val % (A->maxsmp + 1));
        ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, smp);
        val++;
        return FALSE;
      }
  }

