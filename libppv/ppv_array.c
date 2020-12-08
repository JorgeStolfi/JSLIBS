/* See ppv_array.h */
/* Last edited on 2020-12-07 23:24:53 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <ppv_array.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <affirm.h>
    
#define N ppv_array_NAXES
  /* For short. */

bool_t ppv_index_is_valid ( ppv_index_t ix[], ppv_array_t *A )
  { return ix_is_valid(N, ix, A->size); }

void ppv_index_clear ( ppv_index_t ix[] )
  { ix_fill(N, ix, 0); }

void ppv_index_assign ( ppv_index_t ix[], ppv_index_t *val )
  { ix_assign(N, ix, val); }

bool_t ppv_index_first ( ppv_index_t ix[], ppv_array_t *A )
  { return ix_assign_min(N, ix, A->size); } 

bool_t ppv_index_last ( ppv_index_t ix[], ppv_array_t *A )
  { return ix_assign_max(N, ix, A->size); } 

void ppv_index_shift ( ppv_index_t ix[], ppv_index_t *inc )
  { ix_shift(N, ix, inc); } 

sign_t ppv_index_compare ( ppv_index_t ixa[], ppv_index_t ixb[] )
  { return ix_compare(N, ixa, ixb, ix_order_L); }
    
bool_t ppv_index_next ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p )
  { demand(d <= N, "d too big"); 
    return ix_next(d, ix, A->size, ix_order_L, A->step, p, NULL, NULL, NULL, NULL);
  }

bool_t ppv_index_prev ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p )
  { demand(d <= N, "d too big"); 
    return ix_prev(d, ix, A->size, ix_order_L, A->step, p, NULL, NULL, NULL, NULL);
  }

bool_t ppv_descriptor_is_valid ( ppv_array_t *A, bool_t die )
  { assert(ppv_MAX_POS == ix_MAX_POS); 
    assert(ppv_MAX_ABS_STEP == ix_MAX_ABS_STEP);
    assert(ppv_MAX_SIZE == ix_MAX_SIZE);
    if (! ix_parms_are_valid(N, A->size, A->base, A->step, die)) { return FALSE; };
    if (! ix_positions_are_distinct(N, A->size, A->step, die)) { return FALSE; };
    return TRUE;
  }

ppv_array_t ppv_new_array ( ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw )
  { ppv_array_t A;
    
    /* Check {bps} and store it: */
    demand(bps <= ppv_MAX_BPS, "bits-per-sample too big");
    A.bps = bps;
    /* Check {bpw} and store it: */
    demand(bpw <= ppv_MAX_BPW, "bits-per-word too big");
    A.bpw = bpw; 
    /* The {base} is always 0: */
    A.base = 0;
    
    /* Save {A.size}, check if too large, see if array is empty: */
    ppv_axis_t i;
    bool_t empty = FALSE;
    for (i = 0; i < N; i++) 
      { ppv_size_t szi = sz[i]; 
        demand (szi <= ppv_MAX_SIZE, "bad array size");
        A.size[i] = szi;
        if (szi == 0) { empty = TRUE; }
      }
    /* Set {A.step}, {A->el}: */
    if (empty || (bps == 0))
      { /* Memoryless array; reset all steps to zero. */ 
        for (i = 0; i < N; i++) { A.step[i] = 0; }
        /* No storage area: */
        A.el = NULL;
      }
    else
      { /* Array needs storage area. */
        /* Compute number {npos} of distinct positions, set {A.step}: */
        ppv_sample_count_t npos = 1;
        for (i = 0; i < N; i++) 
          { ppv_size_t szi = sz[i];
            if (szi == 1) 
              { A.step[i] = 0; }
            else
              { /* Store position increment along axis {i}, if relevant: */
                A.step[i] = npos;
                /* Check overflow BEFORE multiplication: */
                if (npos > 0) { demand (szi <= ppv_MAX_SAMPLES/npos, "too many samples"); }
                npos *= szi;
              }
          }
        /* Compute number of words {nw}: */
        size_t nw;
        if (bps < bpw)
          { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
            nw = (npos + spw - 1) / spw;
          }
        else if (bps > bpw)
          { ppv_nbits_t wps = (ppv_nbits_t)(bps + bpw - 1)/bpw; /* Words per sample. */
            nw = npos * wps;
          }
        else
          { nw = npos; }
        /* Allocate storage are of correct size: */
        size_t cpw; /* Bytes per word. */
        if (bpw == 8)
          { cpw = sizeof(ppv_word_08_t); }
        else if (bpw == 16)
          { cpw = sizeof(ppv_word_16_t); }
        else if (bpw == 32)
          { cpw = sizeof(ppv_word_32_t); }
        else
          { demand(FALSE, "unsupported bits-per-word"); /* GCC pacifier: */ cpw = 0; }
        demand (nw <= ppv_MAX_BYTES/cpw, "too many bytes");
        A.el = malloc(nw * cpw);
        /* Check for allocation failure: */
        affirm(A.el != NULL, "no mem for new array");
      }
    return A;
  }

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, ppv_index_t ix[] )
  { return ix_position(N, ix, A->base, A->step); }

ppv_sample_t ppv_get_sample 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos 
  )
  { if (bps == 0)
      { return 0; }
  else if (bps < bpw)
      { int spw = bpw/bps; /* Samples per word. */
        /* Fetch the word {qw} containing the sample: */
        ppv_word_t qw;
        if (bpw == 8) 
          { qw = *(((ppv_word_08_t *)el) + pos/spw); }
        else if (bpw == 16) 
          { qw = *(((ppv_word_16_t *)el) + pos/spw); }
        else 
          { qw = *(((ppv_word_32_t *)el) + pos/spw); }
        ppv_nbits_t shift = (ppv_nbits_t)((spw - 1 - pos%spw)*bps); /* Shift to place sample at low end: */
        ppv_word_t mask = ((ppv_word_t)1<<bps) - 1;   /* Mask to extract low-end sample. */
        return (qw >> shift) & mask; 
      }
    else if (bps > bpw)
      { int wps = (bps + bpw - 1) / bpw; /* Words per sample. */
        int k;
        /* Gather the sample words into {qw}: */
        ppv_word_t qw = 0;
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos*wps;
            for (k = 0; k < wps; k++) { qw = (qw << bpw) | (*q8); }
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos*wps;
            for (k = 0; k < wps; k++) { qw = (qw << bpw) | (*q16); }
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos*wps;
            for (k = 0; k < wps; k++) { qw = (qw << bpw) | (*q32); }
          }
        ppv_word_t mask = ((ppv_word_t)1 << bps) - 1; /* Mask to remove padding. */
        return qw & mask;
      }
    else /* bps == bpw */
      { if (bpw == 8) 
          { return *(((ppv_word_08_t *)el) + pos); }
        else if (bpw == 16) 
          { return *(((ppv_word_16_t *)el) + pos); }
        else 
          { return *(((ppv_word_32_t *)el) + pos); }
      }
  }

void ppv_set_sample 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos, 
    ppv_sample_t qv 
  )
  { ppv_word_t mask = ((ppv_word_t)1 << bps) - 1; /* Mask to remove padding. */
    demand(qv <= mask, "bad pixel value");
    if (bps == 0)
      { return; }
    else if (bps < bpw)
      { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
        ppv_nbits_t shift = (ppv_nbits_t)((spw - 1 - pos%spw)*bps); /* Shift to apply to sample: */
        ppv_word_t mask = (((ppv_word_t)1<<bps) - 1) << shift; /* Mask to extract sample. */
        /* Insert the sample into the appropriate word: */
        /* fprintf(stderr, "mask = %u ~mask = %u shift = %d\n", mask, (~mask), shift); */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos/spw; 
            (*q8) = (ppv_word_08_t)(((*q8) & (~mask)) | ((qv << shift) & mask));
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos/spw;
            (*q16) = (ppv_word_16_t)(((*q16) & (~mask)) | ((qv << shift) & mask));
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos/spw;
            (*q32) = (ppv_word_32_t)(((*q32) & (~mask)) | ((qv << shift) & mask));
          }
      }
    else if (bps > bpw)
      { ppv_nbits_t wps = (ppv_nbits_t)((bps + bpw - 1) / bpw); /* Words per sample. */
        ppv_word_t maskw = (1 << bpw) - 1; /* Mask to chop a wordful from the sample. */
        ppv_word_t masks = (1 << (bps - wps*bpw)) - 1; /* Mask to protect padding. */
        int k;
        /* Split the sample into words and store them */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + (pos + wps - 1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q8) = (ppv_word_08_t)(qv & maskw); qv >>= bpw; q8--; }
            (*q8) = (ppv_word_08_t)(((*q8) & (!masks)) | (qv & masks)); 
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + (pos + wps -1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q16) = (ppv_word_16_t)(qv & maskw); qv >>= bpw; q16--; }
            (*q16) =  (ppv_word_16_t)(((*q16) & (!masks)) | (qv & masks)); 
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + (pos + wps -1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q32) =  (ppv_word_32_t)(qv & maskw); qv >>= bpw; q32--; }
            (*q32) =  (ppv_word_32_t)(((*q32) & (!masks)) | (qv & masks)); 
          }
      }
    else /* bps == bpw */
      { if (bpw == 8) 
          { *(((ppv_word_08_t *)el) + pos) = (ppv_word_08_t)qv; }
        else if (bpw == 16) 
          { *(((ppv_word_16_t *)el) + pos) = (ppv_word_16_t)qv; }
        else 
          { *(((ppv_word_32_t *)el) + pos) = (ppv_word_32_t)qv; }
      }
  }

/* DESCRIPTOR MANIPULATION */

void ppv_crop ( ppv_array_t *A, ppv_axis_t i, ppv_size_t skip, ppv_size_t keep )
  { ix_crop(N, A->size, &(A->base), A->step, i, skip, keep);
    /* If the array has become empty, {A.el} is now useless: */
    if (keep == 0) { A->el = NULL; }
  }

void ppv_subsample ( ppv_array_t *A, ppv_axis_t i, ppv_size_t stride )
  { ix_subsample(N, A->size, &(A->base), A->step, i, stride); }

void ppv_flip ( ppv_array_t *A, ppv_axis_t i )
  { ix_flip(N, A->size, &(A->base), A->step, i); }

void ppv_replicate ( ppv_array_t *A, ppv_axis_t i, ppv_size_t sz )
  { ix_replicate(N, A->size, &(A->base), A->step, i, sz); }

void ppv_swap_indices ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j, ppv_dim_t n )
  { ix_swap_indices(N, A->size, &(A->base), A->step, i, j, n); }

void ppv_flip_indices ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j )
  { ix_flip_indices(N, A->size, &(A->base), A->step, i, j); }

void ppv_slice ( ppv_array_t *A, ppv_dim_t n, ppv_axis_t ax[], ppv_index_t ix[] )
  { ix_slice(N, A->size, &(A->base), A->step, n, ax, ix); }

void ppv_chop ( ppv_array_t *A, ppv_axis_t i, ppv_size_t sz, ppv_axis_t j )
  { ix_chop(N, A->size, &(A->base), A->step, i, sz, j); }
  
void ppv_diagonal ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j )
  { ix_diagonal(N, A->size, &(A->base), A->step, i, j); }

/* ELEMENT ENUMERATION */

void ppv_enum 
  ( ppv_op_t *op, 
    bool_t reverse, 
    ppv_array_t *A, 
    ppv_array_t *B, 
    ppv_array_t *C 
  )
  {
    /* Obtain a non-null operand {X}: */
    ppv_array_t *X = A;
    if (X == NULL) { X = B; }
    if (X == NULL) { X = C; } 
    /* If all three operands are null, there is nothing to do: */
    if (X == NULL) { return; }
    /* If two or more operands are non-null, their sizes must agree: */
    if ((B != NULL) && (B != X)) 
      { demand(ix_same_size(N, B->size, X->size, FALSE), "B array size mismatch"); }
    if ((C != NULL) && (C != X)) 
      { demand(ix_same_size(N, C->size, X->size, FALSE), "C array size mismatch"); }
      
    /* Get the bases and steps of non-null arguments: */
    ppv_pos_t bA = (A == NULL ? 0 : A->base);
    ppv_step_t *sA = (A == NULL ? NULL : A->step);
    
    ppv_pos_t bB = (B == NULL ? 0 : B->base);
    ppv_step_t *sB = (B == NULL ? NULL : B->step);
    
    ppv_pos_t bC = (C == NULL ? 0 : C->base);
    ppv_step_t *sC = (C == NULL ? NULL : C->step);
    
    ix_enum(op, N, X->size, ix_order_L, reverse, bA, sA, bB, sB, bC, sC);
  }

/* ERROR MESSAGES */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf )
  {
    fprintf(stderr, "%s", pf); 
    fprintf(stderr, "bps = %u bpw = %u", A->bps, A->bpw);
    ppv_axis_t i;
    for (i = 0; i < N; i++)
      { fprintf(stderr, " ");
        fprintf(stderr, ppv_size_t_FMT, A->size[i]);
        fprintf(stderr, "(×");
        fprintf(stderr, ppv_step_t_FMT, A->step[i]);
        fprintf(stderr, ")");
      }
    fprintf(stderr, " base = ");
    fprintf(stderr, ppv_pos_t_FMT, A->base);
    fprintf(stderr, " el = %lu", ((uint64_t)A->el));
    fprintf(stderr, "%s", sf);
  }
