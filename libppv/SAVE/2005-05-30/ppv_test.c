/* Last edited on 2005-05-28 15:50:26 by stolfi */ 
/* Test of the PPV library. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ppv_array.h>

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
void test_transpose(ppv_array_t *A);
void test_flip(ppv_array_t *A);
void test_diagonal(ppv_array_t *A);
void test_chop(ppv_array_t *A);
void test_replicate(ppv_array_t *A);

void dump_storage(FILE *wr, void *el, int nw, ppv_nbits_t bpw);
void check_size(ppv_size_t *sza, ppv_size_t *szb);

int main (int argn, char **argv)
{
  ppv_size_t sz[ppv_NAX];
  sz[0] = 3; sz[1] = 10; sz[2] = 20; sz[3] = 30; sz[4] = 7; sz[5] = 2;
  ppv_nbits_t bps = 7, bpw = 16;
  
  fprintf(stderr, "bps = %d bpw = %d\n", bps, bpw);

  ppv_array_t A = ppv_new_array(sz,  bps, bpw);

  test_new_array(&A, sz, bps, bpw);
  test_packing(&A);
  test_sample_pos(&A);
  test_crop(&A);
  test_subsample(&A);
  test_transpose(&A);
  test_flip(&A);
  test_diagonal(&A);
  test_chop(&A);
  test_replicate(&A);

  return 0;
}

void test_new_array(ppv_array_t *A, ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw)
{
  fprintf(stderr, "Checking num of axes...\n");
  if (ppv_NAX != 6) { bug("unexpected num of axes - fix test"); }
  
  fprintf(stderr, "Checking ppv_new_array, ppv_valid...\n");
  bool_t verbose = TRUE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  if (! ppv_valid(A)) { bug("ppv_new_array returns invalid desc"); }
  if (A->base != 0) { bug("base = %u", A->base); }
  if (A->bps != bps) { bug("bpw = %u", A->bps); }
  if (A->bpw != bpw) { bug("bpw = %u", A->bpw); }
  /* Check steps for a vanilla array: */
  ppv_tot_size_t stp = 1;
  { ppv_axis_t ax;
    for (ax = 0; ax < ppv_NAX; ax++) 
      { if (A->step[ax] != stp) { bug("step[%d] = %u", ax, A->step[ax]); }
        stp *= sz[ax];
      }
  }
  
  fprintf(stderr, "Checking ppv_inside...\n");
  ppv_index_t ixA[ppv_NAX];
  ppv_index_clear (ixA);
  if (! ppv_inside(A, ixA)) { bug("bogus ppv_inside"); } 

  { ppv_axis_t ax;
    for (ax = 0; ax < ppv_NAX; ax++)
      { ixA[ax] = sz[ax]-1;
        if (! ppv_inside(A, ixA)) { bug("bogus ppv_inside"); } 
        ixA[ax] = sz[ax];
        if (ppv_inside(A, ixA)) { bug("bogus ppv_inside"); } 
        ixA[ax] = 0;
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
    { ppv_sample_t xsmp = (smask & (2*pos)) | 1 | (1<<(A->bps-1));
      if (verbose) { fprintf(stderr, "setting sample %d to %d\n", pos, xsmp); }
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
  ppv_index_t ixA[ppv_NAX];
  for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
    for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
      for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
        for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
          for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
            for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
              { 
                if (! ppv_inside(A, ixA)) { bug("not ppv_inside"); } 
                ppv_pos_t pos = ppv_sample_pos(A, ixA);
                if (pos != xpos) 
                  { bug("ppv_position mismatch pos = %u xpos = %u", pos, xpos); } 
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
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t sh[ppv_NAX]; /* Shift of cropped array along each axis. */
  ppv_size_t sz[ppv_NAX]; /* Size of cropped array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sh} and {sz} for full array: */
      { ppv_axis_t ax; for (ax = 0; ax < ppv_NAX; ax++) { sh[ax] = 0; } } 
      memcpy(sz, A->size, ppv_NAX*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first cropping: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      ppv_size_t skip1 = ((trial + trial/2) & 3)*sz[ax1]/8;
      ppv_size_t keep1 = ((trial + trial/4) & 3)*(sz[ax1] - skip1)/8;
      sh[ax1] += skip1; sz[ax1] = keep1;
      ppv_crop(&B, ax1, skip1, keep1);
      if (verbose) 
        { fprintf(stderr, "crop(%d, %d, %d)", ax1, skip1, keep1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_crop returns an invalid desc"); }
      
      /* Perform second cropping (possibly on same axis): */ 
      ppv_axis_t ax2 = ((trial / ppv_NAX) % ppv_NAX);
      ppv_size_t skip2 = ((trial + trial/2) & 3)*sz[ax2]/8;
      ppv_size_t keep2 = ((trial + trial/5) & 3)*(sz[ax2] - skip2)/8;
      sh[ax2] += skip2; sz[ax2] = keep2;
      ppv_crop(&B, ax2, skip2, keep2);
      if (verbose) 
        { fprintf(stderr, "crop(%d, %d, %d)", ax2, skip2, keep2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_crop returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if ((keep1 == 0) || (keep2 == 0))
        { if (B.el != NULL) { bug("empty array with non-null storage"); } }
      else
        { if (B.el != A->el) { bug("garbled storage area"); } }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
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
                      { bug("crop position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}

void test_subsample(ppv_array_t *A) 
{
  fprintf(stderr, "Checking ppv_subsample...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t st[ppv_NAX]; /* Subsampling step along each axis. */
  ppv_size_t sz[ppv_NAX]; /* Size of subsampled array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {st} and {sz} for full array: */
      { ppv_axis_t ax; for (ax = 0; ax < ppv_NAX; ax++) { st[ax] = 1; } } 
      memcpy(sz, A->size, ppv_NAX*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);
      
      /* Perform first subsampling: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      ppv_size_t step1 = 1 + (trial & 3)*(sz[ax1] - 1)/8;
      ppv_size_t keep1 = (sz[ax1] + step1 - 1)/step1;
      st[ax1] *= step1; sz[ax1] = keep1;
      ppv_subsample(&B, ax1, step1);
      if (verbose) 
        { fprintf(stderr, "subsample(%d, %d)", ax1, step1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_flip returns an invalid desc"); }
      
      /* Perform second subsampling (possibly on same axis): */ 
      ppv_axis_t ax2 = ((trial / ppv_NAX) % ppv_NAX);
      ppv_size_t step2 = 1 + (trial/2 & 3)*(sz[ax2] - 1)/8;
      ppv_size_t keep2 = (sz[ax2] + step2 - 1)/step2;
      st[ax2] *= step2; sz[ax2] = keep2;
      ppv_subsample(&B, ax2, step2);
      if (verbose) 
        { fprintf(stderr, "subsample(%d, %d)", ax2, step2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_subsample returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, sz); 
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                  { 
                    ppv_axis_t ax;
                    for (ax = 0; ax < ppv_NAX; ax++) { ixA[ax] = ixB[ax]*st[ax]; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("crop position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}

void test_transpose(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_transpose...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t tr[ppv_NAX]; /* Axis {i} of {B} is axis {tr[i]} of {A}. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {tr} with identity permutation: */
      { ppv_axis_t ax; for (ax = 0; ax < ppv_NAX; ax++) tr[ax] = ax; }
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first transposition: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      ppv_axis_t bx1 = (trial/2 % ppv_NAX);
      { ppv_axis_t tmp = tr[ax1]; tr[ax1] = tr[bx1]; tr[bx1] = tmp; }
      ppv_transpose(&B, ax1, bx1);
      if (verbose) 
        { fprintf(stderr, "transpose(%d, %d)", ax1, bx1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_flip returns an invalid desc"); }

      /* Perform second transposition (possibly on same axis): */ 
      ppv_axis_t ax2 = ((trial/ppv_NAX) % ppv_NAX);
      ppv_axis_t bx2 = ((trial/ppv_NAX/2) % ppv_NAX);
      { ppv_axis_t tmp = tr[ax2]; tr[ax2] = tr[bx2]; tr[bx2] = tmp; }
      ppv_transpose(&B, ax2, bx2);
      if (verbose) 
        { fprintf(stderr, "transpose(%d, %d)", ax2, bx2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_transpose returns an invalid desc"); }
      
      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      { ppv_axis_t ax;
        for (ax = 0; ax < ppv_NAX; ax++)
          { /* Check whether the size of {B} is correct: */
            if (B.size[ax] != A->size[tr[ax]]) 
              { bug("size[%d] = %u", ax, B.size[ax]); }
            /* Check whether the increments of {B} are correct: */
            if (B.step[ax] != A->step[tr[ax]]) 
              { bug("step[%d] = %u", ax, B.step[ax]); }
          }
      }
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixA[5] = 0; ixA[5] < A->size[5]; ixA[5]++)
        for (ixA[4] = 0; ixA[4] < A->size[4]; ixA[4]++)
          for (ixA[3] = 0; ixA[3] < A->size[3]; ixA[3]++)
            for (ixA[2] = 0; ixA[2] < A->size[2]; ixA[2]++)
              for (ixA[1] = 0; ixA[1] < A->size[1]; ixA[1]++)
                for (ixA[0] = 0; ixA[0] < A->size[0]; ixA[0]++)
                  { 
                    ppv_axis_t ax;
                    for (ax = 0; ax < ppv_NAX; ax++) { ixB[ax] = ixA[tr[ax]]; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("crop position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}

void test_flip(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_flip...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  bool_t fp[ppv_NAX]; /* Tells whether each axis was flipped or not. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {fp} for unflipped array: */
      { ppv_axis_t ax; for (ax = 0; ax < ppv_NAX; ax++) { fp[ax] = FALSE; } }
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Perform first flip: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      fp[ax1] = ! fp[ax1];
      ppv_flip(&B, ax1);
      if (verbose) 
        { fprintf(stderr, "flip(%d)", ax1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_flip returns an invalid desc"); }

      /* Perform second flip (possibly on same axis): */ 
      ppv_axis_t ax2 = ((trial / ppv_NAX) % ppv_NAX);
      fp[ax2] = ! fp[ax2];
      ppv_flip(&B, ax2);
      if (verbose) 
        { fprintf(stderr, "flip(%d)", ax2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_flip returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      /* Check whether the size of {B} is correct: */
      check_size(B.size, A->size);
      /* Check whether the increments of {B} are correct: */
      { ppv_axis_t ax; 
        for (ax = 0; ax < ppv_NAX; ax++) 
          if (B.step[ax] != A->step[ax]*(fp[ax]?-1:+1)) 
            { bug("step[%d] = %u", ax, B.step[ax]); }
      }
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_axis_t ax;
                    for (ax = 0; ax < ppv_NAX; ax++) 
                      { ixA[ax] = (fp[ax] ? A->size[ax] - 1 - ixB[ax] : ixB[ax]); }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("crop position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}

void test_diagonal(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_diagonal...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t sz[ppv_NAX]; /* Size of diagonalized array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, ppv_NAX*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Choose any two axes: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      ppv_axis_t bx1 = (trial/ppv_NAX % ppv_NAX);
      /* Diagonalize array: */
      ppv_diagonal(&B, ax1, bx1);
      /* Compute expected sizes of {1}: */
      if (sz[bx1] < sz[ax1]) { sz[ax1] = sz[bx1]; }
      if (bx1 != ax1) { sz[bx1] = 1; }
      if (verbose) 
        { fprintf(stderr, "diagonal(%d, %d)", ax1, bx1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_diagonal returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      /* Check whether the size of {1} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB);
                    ixA[bx1] = ixA[ax1];
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("diagonal position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}

void test_chop(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_chop...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t sz[ppv_NAX]; /* Size of chopped array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, ppv_NAX*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);

      /* Choose two DISTINCT axes: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      ppv_axis_t bx1 = (trial/ppv_NAX % (ppv_NAX-1));
      if (bx1 >= ax1) { bx1++; }
      /* Choose a nonzero chunk size: */
      ppv_size_t chunksz1 = 1 + (trial/3 & 3)*(A->size[ax1] - 1)/8;
      /* Ensure that the argument of {chop} is trivial along axis {bx1}: */
      ppv_size_t skip1 = A->size[bx1]/2;
      ppv_crop(&B, bx1, skip1, 1);
      /* Chop array: */
      ppv_chop(&B, ax1, chunksz1, bx1);
      /* Compute expected sizes of {1}: */
      sz[bx1] = sz[ax1]/chunksz1; sz[ax1] = chunksz1;
      if (verbose) 
        { fprintf(stderr, "chop(%d, %d, %d)", ax1, chunksz1, bx1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_chop returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      /* Check whether the size of {1} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB);
                    ixA[ax1] = ixA[ax1] + chunksz1 * ixB[bx1];
                    ixA[bx1] = skip1;
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("crop position error posB = %d posA = %d", posB, posA); } 
                  }
    }
}
  
void test_replicate(ppv_array_t *A)
{
  fprintf(stderr, "Checking ppv_replicate...\n");
  bool_t verbose = FALSE;
  if (verbose) { ppv_print_descriptor(stderr, "A = { ", A, " }\n"); } 
  int trial, ntrials = 2*ppv_NAX*ppv_NAX;
  ppv_size_t sz[ppv_NAX]; /* Size of replicated array along each axis. */
  for (trial = 0; trial < ntrials; trial++)
    { /* Initialze {sz} for full array: */
      memcpy(sz, A->size, ppv_NAX*sizeof(ppv_size_t)); 
      /* Start with the standard array: */
      ppv_array_t B = (*A);
      
      /* Choose axis for first replication: */
      ppv_axis_t ax1 = (trial % ppv_NAX);
      /* Choose a positive replication factor: */
      ppv_size_t rep1 = 1 + (trial/3 & 3)*5;
      /* Chop and replicate along axis {ax1}: */
      ppv_size_t skip1 = B.size[ax1]/2;
      ppv_crop(&B, ax1, skip1, 1);
      ppv_replicate(&B, ax1, rep1);
      sz[ax1] = rep1;
      if (verbose) 
        { fprintf(stderr, "chop(%d,%d,1)+replicate(%d,%d)", ax1, skip1, ax1, rep1);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_replicate returns an invalid desc"); }
      
      /* Choose axis for second replication (may be {ax1}): */
      ppv_axis_t ax2 = (trial/ppv_NAX) % ppv_NAX;
      /* Choose a replication factor: */
      ppv_size_t rep2 = 1 + (trial/8 & 3)*5;
      /* Chop and replicate along axis {ax2}: */
      ppv_size_t skip2 = B.size[ax2]/2;
      ppv_crop(&B, ax2, skip2, 1);
      ppv_replicate(&B, ax2, rep2);
      sz[ax2] = rep2;
      if (verbose) 
        { fprintf(stderr, "chop(%d,%d,1)+replicate(%d,%d)", ax2, skip2, ax2, rep2);
          ppv_print_descriptor(stderr, " = { ", &B, " }\n");
        } 
      /* Check validity of descriptor: */
      if (! ppv_valid(&B)) { bug("ppv_replicate returns an invalid desc"); }

      /* Check consistency of storage area: */
      if (B.el != A->el) { bug("garbled storage area"); }
      /* Check whether the size of {1} is correct: */
      check_size(B.size, sz);
      /* Now check coincidence of the two arrays: */
      ppv_index_t ixA[ppv_NAX], ixB[ppv_NAX];
      for (ixB[5] = 0; ixB[5] < B.size[5]; ixB[5]++)
        for (ixB[4] = 0; ixB[4] < B.size[4]; ixB[4]++)
          for (ixB[3] = 0; ixB[3] < B.size[3]; ixB[3]++)
            for (ixB[2] = 0; ixB[2] < B.size[2]; ixB[2]++)
              for (ixB[1] = 0; ixB[1] < B.size[1]; ixB[1]++)
                for (ixB[0] = 0; ixB[0] < B.size[0]; ixB[0]++)
                  { 
                    ppv_index_assign(ixA, ixB);
                    ixA[ax1] = skip1;
                    /* The second crop changes {base} only if {ax1 != ax2}: */
                    if (ax2 != ax1) { ixA[ax2] = skip2; }
                    ppv_pos_t posB = ppv_sample_pos(&B, ixB);
                    ppv_pos_t posA = ppv_sample_pos(A, ixA);
                    if (posB != posA)
                      { bug("replicate position error posB = %d posA = %d", posB, posA); } 
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
  ppv_axis_t ax;
  for (ax = 0; ax < ppv_NAX; ax++)
    { if (sza[ax] != szb[ax]) { bug("size[%d] = %u", ax, sza[0]); } }
}
