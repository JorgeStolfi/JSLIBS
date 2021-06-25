/* Last edited on 2021-06-22 13:44:54 by jstolfi */ 
/* Test of the I/O functions from the PPV library. */

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <bool.h>

#include <ppv_array.h>
#include <ppv_array_write.h>
#include <ppv_array_read.h>

#define posFMT ppv_pos_t_FMT
#define sizeFMT ppv_size_t_FMT

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

ppv_array_t *make_test_array(void);
void test_fill_array(ppv_array_t *A);

void check_array_equality(ppv_array_t *A, ppv_array_t *B);
void dump_storage(FILE *wr, void *el, int nw, ppv_nbits_t bpw);
void check_size(ppv_dim_t d, ppv_size_t *sza, ppv_size_t *szb);

int main (int argn, char **argv)
  {
    char *name = "test";
    
    ppv_array_t *A = make_test_array();

    fprintf(stderr, "Checking ppv_array_write_file...\n");
    bool_t verbose = FALSE;
    if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
    char *fname = NULL;
    asprintf(&fname, "out/%s.ppv", name);
    FILE *wr = open_write(fname, TRUE);
    ppv_array_write_file(wr, A, TRUE);
    ppv_array_write_file(wr, A, FALSE);
    fclose(wr);
    /* Read it back: */
    fprintf(stderr, "Checking ppv_array_read_file...\n");
    ppv_nbits_t new_bpw = A->bpw;  /* For now. */
    ppv_array_t *B;
    FILE *rd = open_read(fname, TRUE);
    for (int32_t i = 0; i < 2; i++)
      { B = ppv_array_read_file(rd, new_bpw);
        if (! ppv_descriptor_is_valid(B, TRUE)) 
          { bug("ppv_array_read_file returns an invalid desc"); }
        if (B->bpw != new_bpw) 
          { bug("ppv_array_read_file returns wrong bpw"); }
        check_array_equality(A, B);
        free(B->el);
      }
    fclose(rd);
    free(fname);
    return 0;
  }

void check_array_equality(ppv_array_t *A, ppv_array_t *B)
  {
    ppv_dim_t d = A->d;
    assert(d == 6);
    demand(B->d == d, "dims don't match");
    /* Check array size: */
    check_size(d, B->size, A->size);
    ppv_size_t *sz = A->size;
    /* Check sample size: */
    if (B->bps != A->bps) { bug("wrong bps"); }
    /* Check contents: */
    ppv_index_t ix[d];
    for (ix[5] = 0; ix[5] < sz[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < sz[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < sz[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < sz[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < sz[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < sz[0]; ix[0]++)
                { 
                  ppv_pos_t posA = ppv_sample_pos(A, ix);
                  ppv_sample_t smpA = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, posA);
                  
                  ppv_pos_t posB = ppv_sample_pos(B, ix);
                  ppv_sample_t smpB = ppv_get_sample_at_pos(B->el, B->bps, B->bpw, posB);
                  
                  if (smpA != smpB)
                    { bug("contents differ at positions posA = " posFMT " posB = " posFMT "", posA, posB); } 
                }
  }

ppv_array_t *make_test_array(void)
  {
    ppv_dim_t d = 6;

    ppv_size_t sz[d];
    sz[0] = 3; sz[1] = 40; sz[2] = 30; sz[3] = 3; sz[4] = 2; sz[5] = 3;
    ppv_nbits_t bps = 7, bpw = 16;

    fprintf(stderr, "bps = %d bpw = %d\n", bps, bpw);

    ppv_array_t *A = ppv_array_new(d, sz,  bps, bpw);

    fprintf(stderr, "creating test array...\n");
    bool_t verbose = TRUE;

    if (verbose) 
      { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 

    if (! ppv_descriptor_is_valid(A, TRUE))
      { bug("{ppv_array_new} returns invalid desc"); }
      
    test_fill_array(A);
    
    return A;
    
  }

void test_fill_array(ppv_array_t *A)
  {
    ppv_dim_t d = A->d;
    ppv_word_t val = 4615;
    ppv_word_t smask = ((ppv_word_t)1 << A->bps) - 1;
    ppv_index_t ix[d];
    for (ix[5] = 0; ix[5] < A->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < A->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < A->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < A->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < A->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < A->size[0]; ix[0]++)
                { 
                  ppv_pos_t pos = ppv_sample_pos(A, ix);
                  ppv_sample_t xsmp = smask & val;
                  ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, xsmp);
                  /* val = (31415*val + 1703) >> 3; */
                  val++;
                }
  }

void dump_storage(FILE *wr, void *el, int nw, ppv_nbits_t bpw)
  {
    for (int32_t iw = 0; iw < nw; iw++)
      { ppv_word_t w;
        if (bpw == 8) 
            { w = *(((ppv_word_08_t *)el) + iw); }
          else if (bpw == 16) 
            { w = *(((ppv_word_16_t *)el) + iw); }
          else 
            { w = *(((ppv_word_32_t *)el) + iw); }
        if (iw != 0) { fprintf(wr, " "); }
        for (int32_t ib = bpw-1; ib >= 0; ib--)
          { fprintf(wr, "%d", (w >> ib) & 1); }
      }
  }

void check_size(ppv_dim_t d, ppv_size_t *sza, ppv_size_t *szb)
  {
    ppv_axis_t i;
    for (i = 0; i < d; i++)
      { if (sza[i] != szb[i]) { bug("size[%d] = " sizeFMT "", i, sza[0]); } }
  }
