/* Last edited on 2011-09-19 22:44:26 by stolfilocal */ 
/* Test of the PPV library. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <ppv_array.h>

#define N ppv_array_NAXES

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

void test_new_array(ppv_array_t *A, ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw);
void test_packing(ppv_array_t *A);
void test_sample_pos(ppv_array_t *A);
void test_crop(ppv_array_t *A);
void test_subsample(ppv_array_t *A);
void test_flip(ppv_array_t *A);
void test_replicate(ppv_array_t *A);
void test_swap_indices(ppv_array_t *A);
void test_flip_indices(ppv_array_t *A);
void test_slice(ppv_array_t *A);
void test_diagonal(ppv_array_t *A);
void test_chop(ppv_array_t *A);

void dump_storage(FILE *wr, void *el, int nw, ppv_nbits_t bpw);
void check_size(ppv_size_t *sza, ppv_size_t *szb);

int main (int argn, char **argv)
{
  ppv_size_t sz[N];
  sz[0] = 3; sz[1] = 10; sz[2] = 20; sz[3] = 30; sz[4] = 7; sz[5] = 2;
  ppv_nbits_t bps = 7, bpw = 16;
  
  fprintf(stderr, "bps = %d bpw = %d\n", bps, bpw);

  ppv_array_t A = ppv_new_array(sz,  bps, bpw);

  test_new_array(&A, sz, bps, bpw);
  test_packing(&A);
  test_sample_pos(&A);
  test_crop(&A);
  test_subsample(&A);
  test_flip(&A);
  test_replicate(&A);
  test_swap_indices(&A);
  test_flip_indices(&A);
  test_slice(&A);
  test_diagonal(&A);
  test_chop(&A);

  return 0;
}

void test_new_array(ppv_array_t *A, ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw)
{
  fprintf(stderr, "Checking num of axes...\n");
  if (N != 6) { bug("unexpected num of axes - fix test"); }
  
  fprintf(stderr, "Checking ppv_new_array, ppv_descriptor_is_valid...\n");
  bool_t verbose = TRUE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  if (! ppv_descriptor_is_valid(A, TRUE)) { bug("ppv_new_array returns invalid desc"); }
  if (A->base != 0) { bug("base = %llu", A->base); }
  if (A->bps != bps) { bug("bpw = %u", A->bps); }
  if (A->bpw != bpw) { bug("bpw = %u", A->bpw); }
  /* Check steps for a vanilla array: */
  ppv_sample_count_t stp = 1;
  { ppv_axis_t i;
    for (i = 0; i < N; i++) 
      { if (A->step[i] != stp) { bug("step[%d] = %lld", i, A->step[i]); }
        stp *= sz[i];
      }
  }
  
  fprintf(stderr, "Checking ppv_index_is_valid...\n");
  ppv_index_t ixA[N];
  ppv_index_clear (ixA);
  if (! ppv_index_is_valid(ixA, A)) { bug("bogus ppv_index_is_valid"); } 

  { ppv_axis_t i;
    for (i = 0; i < N; i++)
      { ixA[i] = sz[i]-1;
        if (! ppv_index_is_valid(ixA, A)) { bug("bogus ppv_index_is_valid"); } 
        ixA[i] = sz[i];
        if (ppv_index_is_valid(ixA, A)) { bug("bogus ppv_index_is_valid"); } 
        ixA[i] = 0;
      }
  }
}

void test_packing(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_set_sample, ppv_get_sample on a few elems...\n");
  bool_t verbose = TRUE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  ppv_word_t smask = ((ppv_word_t)1 << A->bps) - 1;
  int npos = 5; /* Number of positions to test */
  int ndump = (npos*A->bps + A->bpw - 1)/A->bpw + 1; /* Number of words to dump. */
  fprintf(stderr, "\n");
  ppv_pos_t pos;
  for (pos = 0; pos < npos; pos++)
    { ppv_sample_t xsmp = (smask & (2*pos)) | 1 | (1LU <<(A->bps-1));
      if (verbose) { fprintf(stderr, "setting sample %llu to %u\n", pos, xsmp); }
      ppv_set_sample(A->el, A->bps, A->bpw, pos, xsmp);
      if (verbose) { dump_storage(stderr, A->el, ndump, A->bpw); }
      if (verbose) { fprintf(stderr, "\n"); }
      ppv_sample_t smp = ppv_get_sample(A->el, A->bps, A->bpw, pos);
      if (smp != xsmp) 
        { bug("set/get sample mismatch smp = %u xsmp = %u", smp, xsmp); } 
    }
  fprintf(stderr, "\n");
}

void test_sample_pos(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_sample_pos, setting and checking all elems...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  ppv_word_t smask = ((ppv_word_t)1 << A->bps) - 1;
  ppv_pos_t xpos = 0;
  ppv_index_t ixA[N];
  for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
    for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
      for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
        for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
          for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
            for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
              { 
                if (! ppv_index_is_valid(ixA, A)) { bug("not ppv_index_is_valid"); } 
                ppv_pos_t pos = ppv_sample_pos(A, ixA);
                if (pos != xpos) 
                  { bug("ppv_sample_pos mismatch pos = %llu xpos = %llu", pos, xpos); } 
                ppv_sample_t xsmp = smask & xpos;
                ppv_set_sample(A->el, A->bps, A->bpw, pos, xsmp);
                ppv_sample_t smp = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                if (smp != xsmp) 
                  { bug("set/get sample mismatch smp = %u xsmp = %u", smp, xsmp); } 

                xpos++;
              }

  fprintf(stderr, "Re-checking the value of all samples...\n");
  for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
    for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
      for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
        for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
          for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
            for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
              { 
                ppv_pos_t pos = ppv_sample_pos(A, ixA);
                ppv_sample_t xsmp = smask & pos;
                ppv_sample_t smp = ppv_get_sample(A->el, A->bps, A->bpw, pos);
                if (smp != xsmp) 
                  { bug("set/get sample mismatch smp = %u xsmp = %u", smp, xsmp); } 
              }
}

void test_crop(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_crop...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_step_t sh[N]; /* Shift of cropped array along each axis. */
  ppv_size_t sz[N]; /* Size of cropped array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sh} and {sz} for full array: */
      { ppv_axis_t i; for (i = 0; i < N; i++) { sh[i] = 0; } } 
      memcpy(sz, A->size, N*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first cropping: */
      ppv_axis_t i1 = (trial % N);
      ppv_size_t skip1 = ((trial + trial/2) & 3)*sz[i1]/8;
      ppv_size_t keep1 = ((trial + trial/4) & 3)*(sz[i1] - skip1)/8;
      sh[i1] += skip1; sz[i1] = keep1;
      ppv_crop(&B, i1, skip1, keep1);
      if (verbose) 
        { fprintf(stderr, "ppv_crop(%d, %llu, %llu)", i1, skip1, keep1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_crop returns an invalid desc"); }
      
      /* Perform second cropping (possibly on same axis): */ 
      ppv_axis_t i2 = ((trial / N) % N);
      ppv_size_t skip2 = ((trial + trial/2) & 3)*sz[i2]/8;
      ppv_size_t keep2 = ((trial + trial/5) & 3)*(sz[i2] - skip2)/8;
      sh[i2] += skip2; sz[i2] = keep2;
      ppv_crop(&B, i2, skip2, keep2);
      if (verbose) 
        { fprintf(stderr, "ppv_crop(%d, %llu, %llu)", i2, skip2, keep2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_crop returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if ((keep1 == 0) || (keep2 == 0))
        { if (B.el != NULL) { bug("ppv_crop empty array with non-null storage"); } }
      else
        { if (B.el != A->el) { bug("ppv_crop garbled storage area"); } }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB); ppv_index_shift(ixA, sh);
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_crop position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}

void test_subsample(ppv_array_t *A) 
{
  fprintf(stderr, "Checking ppv_subsample...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t st[N]; /* Subsampling step along each axis. */
  ppv_size_t sz[N]; /* Size of subsampled array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {st} and {sz} for full array: */
      { ppv_axis_t i; for (i = 0; i < N; i++) { st[i] = 1; } } 
      memcpy(sz, A->size, N*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);
      
      /* Perform first subsampling: */
      ppv_axis_t i1 = (trial % N);
      ppv_size_t step1 = 1 + (trial & 3)*(sz[i1] - 1)/8;
      ppv_size_t keep1 = (sz[i1] + step1 - 1)/step1;
      st[i1] *= step1; sz[i1] = keep1;
      ppv_subsample(&B, i1, step1);
      if (verbose) 
        { fprintf(stderr, "ppv_subsample(%d, %llu)", i1, step1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_subsample returns an invalid desc"); }
      
      /* Perform second subsampling (possibly on same axis): */ 
      ppv_axis_t i2 = ((trial / N) % N);
      ppv_size_t step2 = 1 + (trial/2 & 3)*(sz[i2] - 1)/8;
      ppv_size_t keep2 = (sz[i2] + step2 - 1)/step2;
      st[i2] *= step2; sz[i2] = keep2;
      ppv_subsample(&B, i2, step2);
      if (verbose) 
        { fprintf(stderr, "ppv_subsample(%d, %llu)", i2, step2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_subsample returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_subsample garbled storage area"); }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, sz); 
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                  { 
                    ppv_axis_t i;
                    for (i = 0; i < N; i++) { ixA[i] = ixB[i]*st[i]; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_subsample position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}

void test_swap_indices(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_swap_indices...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t tr[N]; /* Axis {i} of {B} is axis {tr[i]} of {A}. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {tr} with identity permutation: */
      { ppv_axis_t i; for (i = 0; i < N; i++) tr[i] = i; }
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      int pass;
      for (pass = 0; pass < 2; pass++)
        { /* Perform a transposition: */
          ppv_axis_t i1 = (trial/(1+pass)) % N;
          ppv_axis_t j1 = (trial/(2+pass)) % N;
          int d1 = (i1 < j1 ? j1 - i1 : i1 - j1);
          int m1 = N - (i1 > j1 ? i1 : j1);
          ppv_dim_t n1 = (trial/(3+pass)) % (d1 == 0 ? m1 : (d1 < m1 ? d1 : m1));
          assert(i1 + n1 <= N);
          assert(j1 + n1 <= N);
          assert((i1 == j1) || (i1 + n1 <= j1) || (j1 + n1 <= i1));
          { int k;
            for (k = 0; k < n1; k++)
              { ppv_axis_t tmp = tr[i1+k]; tr[i1+k] = tr[j1+k]; tr[j1+k] = tmp; } 
          }
          ppv_swap_indices(&B, i1, j1, n1);
          if (verbose) 
            { fprintf(stderr, "ppv_swap_indices(%d, %d, %d)", i1, j1, n1);
              ppv_print_descriptor(stderr, " = { ", &B, " }\n");
            } 
          /* Check validity of descriptor: */
          if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_swap_indices returns an invalid desc"); }
        }
      
      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_swap_indices garbled storage area"); }
      { ppv_axis_t i;
        for (i = 0; i < N; i++)
          { /* Check whether the size of {B} is correct: */
            if (B.size[i] != A->size[tr[i]]) 
              { bug("size[%d] = %llu", i, B.size[i]); }
            /* Check whether the increments of {B} are correct: */
            if (B.step[i] != A->step[tr[i]]) 
              { bug("step[%d] = %lld", i, B.step[i]); }
          }
      }
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
        for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
          for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
            for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
              for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
                for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
                  { 
                    ppv_axis_t i;
                    for (i = 0; i < N; i++) { ixB[i] = ixA[tr[i]]; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_swap_indices position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}

void test_flip_indices(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_flip_indices...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t tr[N]; /* Axis {i} of {B} is axis {tr[i]} of {A}. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {tr} with identity permutation: */
      { ppv_axis_t i; for (i = 0; i < N; i++) tr[i] = i; }
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first transposition: */
      ppv_axis_t i1 = (trial % N);
      ppv_axis_t j1 = i1 + (trial/2 % (N-i1));
      { ppv_axis_t ia=i1, ja=j1;
        while (ia < ja) 
          { ppv_axis_t tmp = tr[ia]; tr[ia] = tr[ja]; tr[ja] = tmp; ia++; ja--; }
      }
      ppv_flip_indices(&B, i1, j1);
      if (verbose) 
        { fprintf(stderr, "ppv_flip_indices(%d, %d)", i1, j1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_flip_indices returns an invalid desc"); }

      /* Perform second transposition (possibly on same axis): */ 
      ppv_axis_t i2 = ((trial/N) % N);
      ppv_axis_t j2 = i2 + ((trial/N/2) % (N-i2));
      { ppv_axis_t ia=i2, ja=j2;
        while (ia < ja) 
          { ppv_axis_t tmp = tr[ia]; tr[ia] = tr[ja]; tr[ja] = tmp; ia++; ja--; }
      }
      ppv_flip_indices(&B, i2, j2);
      if (verbose) 
        { fprintf(stderr, "ppv_flip_indices(%d, %d)", i2, j2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_flip_indices returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_flip_indices garbled storage area"); }
      { ppv_axis_t i;
        for (i = 0; i < N; i++)
          { /* Check whether the size of {B} is correct: */
            if (B.size[i] != A->size[tr[i]]) 
              { bug("size[%d] = %llu", i, B.size[i]); }
            /* Check whether the increments of {B} are correct: */
            if (B.step[i] != A->step[tr[i]]) 
              { bug("step[%d] = %lld", i, B.step[i]); }
          }
      }
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
        for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
          for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
            for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
              for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
                for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
                  { 
                    ppv_axis_t i;
                    for (i = 0; i < N; i++) { ixB[i] = ixA[tr[i]]; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_flip_indices position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}

void test_flip(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_flip...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  bool_t fp[N]; /* Tells whether each axis was flipped or not. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {fp} for unflipped array: */
      { ppv_axis_t i; for (i = 0; i < N; i++) { fp[i] = FALSE; } }
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first flip: */
      ppv_axis_t i1 = (trial % N);
      fp[i1] = ! fp[i1];
      ppv_flip(&B, i1);
      if (verbose) 
        { fprintf(stderr, "ppv_flip(%d)", i1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_flip returns an invalid desc"); }

      /* Perform second flip (possibly on same axis): */ 
      ppv_axis_t i2 = ((trial / N) % N);
      fp[i2] = ! fp[i2];
      ppv_flip(&B, i2);
      if (verbose) 
        { fprintf(stderr, "ppv_flip(%d)", i2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_flip returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_flip garbled storage area"); }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, A->size);
      /* Check whether the increments of {B} are correct: */
      { ppv_axis_t i; 
        for (i = 0; i < N; i++) 
          if (B.step[i] != A->step[i]*(fp[i]?-1:+1)) 
            { bug("step[%d] = %lld", i, B.step[i]); }
      }
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_axis_t i;
                    for (i = 0; i < N; i++) 
                      { ixA[i] = (fp[i] ? A->size[i] - 1 - ixB[i] : ixB[i]); }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_flip position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}

void test_slice(ppv_array_t *A)
{
  fprintf(stderr, "NOT checking ppv_slice...\n");
}

void test_diagonal(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_diagonal...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t sz[N]; /* Size of diagonalized array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, N*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Choose any two distinct axes {i1,j1}: */
      ppv_axis_t i1 = (trial % N);
      ppv_axis_t j1 = (trial/N % (N-1));
      if (j1 >= i1) { j1++; }
      /* Make sure {sz[i1] <= sz[j1]}: */
      if (sz[i1] > sz[j1]) { ppv_axis_t t = i1; i1 = j1; j1 = t; }
      if (sz[i1] > 0)
        { 
          /* Take diagonal slice of array: */
          ppv_diagonal(&B, i1, j1);
          if (verbose) 
            { fprintf(stderr, "ppv_diagonal(%d, %d)", i1, j1);
              ppv_print_descriptor(stderr, " = { ", &B, " }\n");
            } 
          /* Compute expected sizes of {1}: */
          sz[j1] = sz[j1] - (sz[i1] - 1);
          /* Check validity of descriptor: */
          if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_diagonal returns an invalid desc"); }

          /* Check consistency of storage area: */
          if (B.el != A->el) { bug("ppv_diagonal garbled storage area"); }
          /* Check whether the size of {1} is correct: */
          check_size(B.size, sz);
          /* Now check coincidence of the two arrays: */
          ppv_index_t ixA[N], ixB[N];
          for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
            for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
              for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
                for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
                  for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                    for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                      { 
                        ppv_index_assign(ixA, ixB);
                        ixA[j1] += ixA[i1];
                        ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                        ppv_pos_t posA = ppv_sample_pos(A, ixA);
                        if (posB != posA)
                          { bug("ppv_diagonal position error posB = %llu posA = %llu", posB, posA); } 
                      }
        }
    }
}

void test_chop(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_chop...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t sz[N]; /* Size of chopped array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, N*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Choose two DISTINCT axes: */
      ppv_axis_t i1 = (trial % N);
      ppv_axis_t j1 = (trial/N % (N-1));
      if (j1 >= i1) { j1++; }
      /* Choose a nonzero chunk size: */
      ppv_size_t chunksz1 = 1 + (trial/3 & 3)*(B.size[i1] - 1)/8;
      /* Ensure that {B.size[i1]} is a multiple of {chunksz1}: */
      if ((B.size[i1] % chunksz1) != 0)
        { if (chunksz1 > B.size[i1]) 
            { chunksz1 = (B.size[i1] + 1)/2; }
          ppv_crop(&B, i1, 0, (B.size[i1]/chunksz1)*chunksz1);
        }
      /* Ensure that {B} is trivial along axis {j1}: */
      ppv_size_t skip1 = B.size[j1]/2;
      ppv_crop(&B, j1, skip1, 1);
      /* Chop array: */
      ppv_chop(&B, i1, chunksz1, j1);
      /* Compute expected sizes of {B}: */
      sz[j1] = sz[i1]/chunksz1; sz[i1] = chunksz1;
      if (verbose) 
        { fprintf(stderr, "ppv_chop(%d, %llu, %d)", i1, chunksz1, j1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_chop returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_chop garbled storage area"); }
      /* Check whether the size of {1} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB);
                    ixA[i1] = ixA[i1] + chunksz1 * ixB[j1];
                    ixA[j1] = skip1;
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_chop position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}
  
void test_replicate(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_replicate...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*N*N;
  ppv_size_t sz[N]; /* Size of replicated array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, N*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);
      
      /* Choose axis for first replication: */
      ppv_axis_t i1 = (trial % N);
      /* Choose a positive replication factor: */
      ppv_size_t rep1 = 1 + (trial/3 & 3)*5;
      /* Chop and replicate along axis {i1}: */
      ppv_size_t skip1 = B.size[i1]/2;
      ppv_crop(&B, i1, skip1, 1);
      ppv_replicate(&B, i1, rep1);
      sz[i1] = rep1;
      if (verbose) 
        { fprintf(stderr, "ppv_crop(%d,%llu,1)+ppv_replicate(%d,%llu)", i1, skip1, i1, rep1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_replicate returns an invalid desc"); }
      
      /* Choose axis for second replication (may be {i1}): */
      ppv_axis_t i2 = (trial/N) % N;
      /* Choose a replication factor: */
      ppv_size_t rep2 = 1 + (trial/8 & 3)*5;
      /* Chop and replicate along axis {i2}: */
      ppv_size_t skip2 = B.size[i2]/2;
      ppv_crop(&B, i2, skip2, 1);
      ppv_replicate(&B, i2, rep2);
      sz[i2] = rep2;
      if (verbose) 
        { fprintf(stderr, "ppv_crop(%d,%llu,1)+ppv_replicate(%d,%llu)", i2, skip2, i2, rep2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_descriptor_is_valid(&B, TRUE)) { bug("ppv_replicate returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("ppv_replicate garbled storage area"); }
      /* Check whether the size of {1} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[N], ixB[N];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB);
                    ixA[i1] = skip1;
                    /* The second crop changes {base} only if {i1 != i2}: */
                    if (i2 != i1) { ixA[i2] = skip2; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("ppv_replicate position error posB = %llu posA = %llu", posB, posA); } 
                  }
    }
}
  
void dump_storage(FILE *wr, void *el, int nw, ppv_nbits_t bpw)
{
  int iw, ib;
  for (iw = 0; iw < nw; iw++)
    { ppv_word_t w;
      if (bpw == 8) 
          { w = *(((ppv_word_08_t *)el) + iw); }
        else if (bpw == 16) 
          { w = *(((ppv_word_16_t *)el) + iw); }
        else 
          { w = *(((ppv_word_32_t *)el) + iw); }
      if (iw != 0) { fprintf(wr, " "); }
      for (ib = bpw-1; ib >= 0; ib--)
        { fprintf(wr, "%d", (w >> ib) & 1); }
    }
}

void check_size(ppv_size_t *sza, ppv_size_t *szb)
{
  ppv_axis_t i;
  for (i = 0; i < N; i++)
    { if (sza[i] != szb[i]) { bug("size[%d] = %llu", i, sza[0]); } }
}
