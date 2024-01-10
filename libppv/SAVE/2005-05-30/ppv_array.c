/* See ppv_array.h */
/* Last edited on 2005-05-28 15:50:09 by stolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <ppv_array.h>
#include <stdlib.h>
#include <stdio.h>
    
bool_t ppv_inside ( ppv_array_t *A, ppv_index_t *ix )
  { return 
      (ix[0] >= 0) && (ix[0] < A->size[0]) &&
      (ix[1] >= 0) && (ix[1] < A->size[1]) &&
      (ix[2] >= 0) && (ix[2] < A->size[2]) &&
      (ix[3] >= 0) && (ix[3] < A->size[3]) &&
      (ix[4] >= 0) && (ix[4] < A->size[4]) &&
      (ix[5] >= 0) && (ix[5] < A->size[5]);
  }

void ppv_index_clear ( ppv_index_t *ix )
  {
    ix[0] = 0;
    ix[1] = 0;
    ix[2] = 0;
    ix[3] = 0;
    ix[4] = 0;
    ix[5] = 0;
  }

void ppv_index_assign ( ppv_index_t *ix, ppv_index_t *val )
  {
    ix[0] = val[0];
    ix[1] = val[1];
    ix[2] = val[2];
    ix[3] = val[3];
    ix[4] = val[4];
    ix[5] = val[5];
  }

void ppv_index_shift ( ppv_index_t *ix, ppv_index_t *inc )
  {
    ix[0] += inc[0];
    ix[1] += inc[1];
    ix[2] += inc[2];
    ix[3] += inc[3];
    ix[4] += inc[4];
    ix[5] += inc[5];
  }

sign_t ppv_index_cmp ( ppv_index_t *ixa, ppv_index_t *ixb )
  {
    if (ixa[0] < ixb[0]) 
      { return NEG; } 
    else if (ixa[0] > ixb[0]) 
      { return POS; } 
    if (ixa[1] < ixb[1]) 
      { return NEG; } 
    else if (ixa[1] > ixb[1]) 
      { return POS; } 
    if (ixa[2] < ixb[2]) 
      { return NEG; } 
    else if (ixa[2] > ixb[2]) 
      { return POS; } 
    if (ixa[3] < ixb[3]) 
      { return NEG; } 
    else if (ixa[3] > ixb[3]) 
      { return POS; } 
    if (ixa[4] < ixb[4]) 
      { return NEG; } 
    else if (ixa[4] > ixb[4]) 
      { return POS; } 
    if (ixa[5] < ixb[5]) 
      { return NEG; } 
    else if (ixa[5] > ixb[5]) 
      { return POS; } 
    return 0;
  }
    
#define ppv_INVALID ppv_error("invalid desc")
/* #define ppv_INVALID return FALSE */

bool_t ppv_valid ( ppv_array_t *A )
  { /* Get and check number of distinct positions: */
    ppv_tot_size_t npos = A->npos;
    if (npos > ppv_MAX_SAMPLES) { ppv_INVALID; }
    /* Get and check the position of element {[0,0,0,0,0]}: */
    ppv_pos_t base = A->base;
    if (base > ppv_MAX_POS) { ppv_INVALID; }
    /* Gather and check the number of bits per sample: */
    ppv_nbits_t bps = A->bps;
    if (bps > ppv_MAX_BPS) { ppv_INVALID; }
    /* Gather and check the number of bits per word: */
    ppv_nbits_t bpw = A->bpw;
    if ((bpw != 8) && (bpw != 16) && (bpw != 32)) { ppv_INVALID; }
    /* Memoryless arrays should have no storage and {base==0}: */
    if ((bps == 0) || (npos == 0))
      { if (A->el != NULL) { ppv_INVALID; }
        if (base != 0)  { ppv_INVALID; }
      }
    /* Gather and check {size,step}: */
    ppv_step_t step[ppv_NAX]; /* Steps in increasing order */
    ppv_size_t size[ppv_NAX]; /* Sizes in the same order */
    int ax, nx=0;
    for (ax = 0; ax < ppv_NAX; ax++) 
      { ppv_step_t st = A->step[ax];
        ppv_size_t sz = A->size[ax];
        if (labs(st) > ppv_MAX_ABS_STEP) { ppv_INVALID; }
        if (sz > ppv_MAX_SIZE) { ppv_INVALID; }
        /* An empty array must have {npos == 0}: */
        if ((sz == 0) && (npos != 0)) { ppv_INVALID; }
        if (st != 0) 
          { /* Memoryless arrays and trivial indices must have {step==0}: */
            if ((npos == 0) || (bps == 0) || (sz == 1)) { ppv_INVALID; }
            /* If the index is not replicated, its size must be a factor of {npos}: */
            if ((npos % sz) != 0) { ppv_INVALID; }
            npos /= sz;
            if (st < 0) 
              { /* Flip index to ensure {st >= 0}: */
                st = -st;
                if (base/(sz - 1) < st) { ppv_INVALID; }
                base -= (sz - 1)*st;
              }
            /* Save {size,step} in increasing order: */
            int bx = nx;
            while ((bx > 0) && (step[bx-1] > st))
              { step[bx] = step[bx-1];
                size[bx] = size[bx-1];
                bx--;
              }
            step[bx] = st;
            size[bx] = sz;
            nx++;
          }
      }
    if (npos > 1) { ppv_INVALID; }
    /* Now check whether all non-replicated positions are distinct and valid: */
    npos = 1;
    for (ax = 0; ax < nx; ax++)
      { if (step[ax] < npos) { ppv_INVALID; }
        npos *= size[ax];
      }
    /* Well, everything seems allright... */
    return TRUE;
  }

ppv_array_t ppv_new_array ( ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw )
  { ppv_array_t A;
    /* Save the size vector: */
    A.size[0] = sz[0];
    A.size[1] = sz[1];
    A.size[2] = sz[2];
    A.size[3] = sz[3];
    A.size[4] = sz[4];
    A.size[5] = sz[5];
    /* Check {bps} and store it: */
    if (bps > ppv_MAX_BPS) { ppv_error("bits-per-sample too big"); }
    A.bps = bps;
    /* Check {bpw} and store it: */
    if (bpw > ppv_MAX_BPW) { ppv_error("bits-per-word too big"); }
    A.bpw = bpw; 
    /* The {base} is always 0: */
    A.base = 0;
    /* Compute the total num of samples {npos} and the natural {step} vector: */
    if ((sz[0] == 0) || (sz[1] == 0) || (sz[2] == 0) || (sz[3] == 0) || (sz[4] == 0) || (sz[5] == 0))
      { A.npos = 0; }
    else if (bps == 0)
      { A.npos = 1; }
    else
      { /* Array needs storage area. */
        /* Compute number of samples {npos} and {step[ax]} for each {ax}: */
        ppv_tot_size_t npos = 1;
        int ax;
        for (ax = 0; ax < ppv_NAX; ax++) 
          { ppv_size_t sz = A.size[ax];
            if (sz == 1) 
              { A.step[ax] = 0; }
            else
              { if (sz > ppv_MAX_SIZE) { ppv_error("bad array size"); }
                /* Store position increment along axis {ax}, if relevant: */
                A.step[ax] = npos;
                /* Check overflow BEFORE multiplication: */
                if ((npos > 0) && (sz > ppv_MAX_SAMPLES/npos)) 
                  { ppv_error("too many samples"); }
                npos *= sz;
              }
          }
        A.npos = npos;
        /* Compute number of words {nw}: */
        size_t nw;
        if (bps < bpw)
          { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
            nw = (npos + spw - 1) / spw;
          }
        else if (bps > bpw)
          { ppv_nbits_t wps = (bps + bpw - 1)/bpw; /* Words per sample. */
            nw = npos * wps;
          }
        else
          { nw = npos; }
        /* Allocate storage are of correct size: */
        if (nw > ppv_MAX_BYTES/(bpw/8)) 
          { ppv_error("too many bytes"); }
        if (bpw == 8)
          { A.el = malloc(nw * sizeof(ppv_word_08_t)); }
        else if (bpw == 16)
          { A.el = malloc(nw * sizeof(ppv_word_16_t)); }
        else if (bpw == 32)
          { A.el = malloc(nw * sizeof(ppv_word_32_t)); }
        else
          { ppv_error("unsupported bits-per-word"); }
        /* Check for allocation failure: */
        if (A.el == NULL) { ppv_error("no mem for new ppv_array_t"); }
        return A;
      }
    /* Array is memoryless: */
    A.step[0] = 0;
    A.step[1] = 0;
    A.step[2] = 0;
    A.step[3] = 0;
    A.step[4] = 0;
    A.el = NULL;
    return A;
  }

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, ppv_index_t *ix )
  { return
      A->base + 
      ix[0]*A->step[0] + 
      ix[1]*A->step[1] + 
      ix[2]*A->step[2] + 
      ix[3]*A->step[3] + 
      ix[4]*A->step[4] + 
      ix[5]*A->step[5];
  }

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
        ppv_nbits_t shift = (spw - 1 - pos%spw)*bps; /* Shift to place sample at low end: */
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
    if (qv > mask) { ppv_error("bad pixel value"); }
    if (bps == 0)
      { return; }
    else if (bps < bpw)
      { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
        ppv_nbits_t shift = (spw - 1 - pos%spw)*bps; /* Shift to apply to sample: */
        ppv_word_t mask = (((ppv_word_t)1<<bps) - 1) << shift; /* Mask to extract sample. */
        /* Insert the sample into the appropriate word: */
        /* fprintf(stderr, "mask = %u ~mask = %u shift = %d\n", mask, (~mask), shift); */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos/spw; 
            (*q8) = ((*q8) & (~mask)) | ((qv << shift) & mask);
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos/spw;
            (*q16) = ((*q16) & (~mask)) | ((qv << shift) & mask);
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos/spw;
            (*q32) = ((*q32) & (~mask)) | ((qv << shift) & mask);
          }
      }
    else if (bps > bpw)
      { ppv_nbits_t wps = (bps + bpw - 1) / bpw; /* Words per sample. */
        ppv_word_t maskw = (1 << bpw) - 1; /* Mask to chop a wordful from the sample. */
        ppv_word_t masks = (1 << (bps - wps*bpw)) - 1; /* Mask to protect padding. */
        int k;
        /* Split the sample into words and store them */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + (pos + wps - 1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q8) = (qv & maskw); qv >>= bpw; q8--; }
            (*q8) = ((*q8) & (!masks)) | (qv & masks); 
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + (pos + wps -1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q16) = (qv & maskw); qv >>= bpw; q16--; }
            (*q16) =  ((*q16) & (!masks)) | (qv & masks); 
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + (pos + wps -1)*wps;
            for (k = 1; k < wps; k++) 
              { (*q32) =  (qv & maskw); qv >>= bpw; q32--; }
            (*q32) =  ((*q32) & (!masks)) | (qv & masks); 
          }
      }
    else /* bps == bpw */
      { if (bpw == 8) 
          { *(((ppv_word_08_t *)el) + pos) = qv; }
        else if (bpw == 16) 
          { *(((ppv_word_16_t *)el) + pos) = qv; }
        else 
          { *(((ppv_word_32_t *)el) + pos) = qv; }
      }
  }

/* DESCRIPTOR MANIPULATION */

void ppv_crop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t skip, ppv_size_t sz )
  { /* Check axis: */
    if (ax >= ppv_NAX) { ppv_error("bad axis"); }
    /* Save the original size along axis {ax}: */
    ppv_size_t osz = A->size[ax];
    /* Check and store the new size: */
    if (sz > osz) { ppv_error("bad new size"); }
    A->size[ax] = sz;
    /* Check the starting index: */
    if (skip > osz - sz) { ppv_error("bad start index"); }
    if (sz == 0)
      { /* Array becomes empty: */
        if (A->npos > 0) 
          { /* Array was non-empty, make it empty: */
            A->base = 0;
            if (A->npos > 1) 
              { /* Array was non-trivial, normalize all steps: */
                A->step[0] = 0;
                A->step[1] = 0;
                A->step[2] = 0;
                A->step[3] = 0;
                A->step[4] = 0;
                A->step[5] = 0;
              }
            A->npos = 0;
            A->el = NULL;
          }
      }
    else if (A->step[ax] > 0)
      { /* Update {A->npos}, {A->step[ax]}: */
        /* This can't overflow because {skip*A->step[ax]} is a valid position: */
        A->base += skip * A->step[ax];
        /* Eliminate the {osz} factor from {npos}: */
        A->npos /= osz;
        if (sz > 1) 
          { /* Multiply the {sz} factor into {npos}: */
            A->npos *= sz;
          }
        else
          { /* Normalize step to 0: */
            A->step[ax] = 0;
          }
      }
  }

void ppv_subsample ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t stride )
  { /* Check axis: */
    if (ax >= ppv_NAX) { ppv_error("bad axis"); }
    /* Check stride: */
    if (stride == 0) { ppv_error("bad stride"); }
    /* Check for trivial subsampling: */
    if (stride == 1) { return; }
    /* Save the original size along axis {ax}: */
    ppv_size_t osz = A->size[ax];
    /* If array is empty or trivial along that index, there is nothing to do: */
    if (osz <= 1) { return; }
    /* Compute and store the new size: */
    ppv_size_t sz = (osz - 1)/stride + 1;
    A->size[ax] = sz;
    /* If array is empty or trivial, we are done: */
    if (A->npos <= 1) { return; }
    /* Eliminate the {osz} factor from {npos}: */
    if (A->step[ax] != 0) { A->npos /= osz; }
    if (sz == 1) 
      { /* Normalize step to 0: */
        A->step[ax] = 0;
      }
    else
      { /* Update the step (this should never overflow): */
        A->step[ax] *= stride;
        /* Multiply the {sz} factor into {npos}: */
        if (A->step[ax] != 0) { A->npos *= sz; }
      }
  }

void ppv_replicate ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz )
  { /* Check axis: */
    if (ax >= ppv_NAX) { ppv_error("bad axis"); }
    /* Check the original size along axis {ax}: */
    if (A->size[ax] != 1) { ppv_error("bad orig size"); }
    /* Check and set the new size: */
    if ((sz == 0) || (sz > ppv_MAX_SIZE)) { ppv_error("bad new size"); }
    A->size[ax] = sz;
    /* The increment must already be zero, but just in case... */
    A->step[ax] = 0;
    /* This operation does not change {npos}. */
  }

void ppv_transpose ( ppv_array_t *A, ppv_axis_t ax, ppv_axis_t bx )
  { /* Check axes: */
    if (ax >= ppv_NAX) { ppv_error("bad axis a"); }
    if (bx >= ppv_NAX) { ppv_error("bad axis b"); }
    /* Trivial transpose: */
    if (ax == bx) { return; }
    /* Exchange sizes and steps: */
    { ppv_size_t tmp = A->size[ax]; A->size[ax] = A->size[bx]; A->size[bx] = tmp; }
    { ppv_step_t tmp = A->step[ax]; A->step[ax] = A->step[bx]; A->step[bx] = tmp; }
    /* This operation does not change {npos}. */
  }

void ppv_flip ( ppv_array_t *A, ppv_axis_t ax )
  { /* Check axis: */
    if (ax >= ppv_NAX) { ppv_error("bad axis"); }
    /* If there is no motion along that axis, we have nothing to d: */
    ppv_step_t st = A->step[ax];
    if (st == 0) { return; }
    /* Now {size[ax]} must be nonzero. */
    /* Adjust {base} to the last valid index value: */
    A->base += (A->size[ax]-1)*st; 
    /* Negate {step[ax]}: */
    A->step[ax] = -st;
    /* This operation does not change {npos}. */
  }

void ppv_chop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz, ppv_axis_t bx )
  { /* Check axes: */
    if (ax >= ppv_NAX) { ppv_error("bad axis a"); }
    if (bx >= ppv_NAX) { ppv_error("bad axis b"); }
    if (ax == bx) { ppv_error("axes must be different"); }
    /* The chunk size must be positive: */
    if (sz < 1) { ppv_error("bad chunk size"); }
    /* The array must be trivial in direction {bx}: */
    if (A->size[bx] != 1) { ppv_error("axis b not trivial"); }
    /* Eliminate factor {size[ax]} from {npos}: */
    if (A->step[ax] != 0) { A->npos /= A->size[ax]; }
    /* Compute new sizes along {bx} and {ax}: */
    A->size[bx] = A->size[ax]/sz;  A->size[ax] = sz;
    /* Compute new stesp along {bx} and {ax}: */
    if (A->size[bx] == 0) 
      { /* Resuting array is empty: */
        if (A->npos > 0) 
          { /* Array was non-empty, make it empty: */
            A->base = 0;
            if (A->npos > 1) 
              { /* Array was non-trivial, normalize all steps: */
                A->step[0] = 0;
                A->step[1] = 0;
                A->step[2] = 0;
                A->step[3] = 0;
                A->step[4] = 0;
                A->step[5] = 0;
              }
            A->npos = 0;
            A->el = NULL;
          }
        else
          { /* Array was already empty, assume all steps are zero. */ }
      } 
    else 
      { if (A->size[bx] == 1)
          { /* Array is still trivial along {bx}, assume {step[bx] = 0} */ }
        else
          { /* This computation cannot overflow: */
            A->step[bx] = sz * A->step[ax];
          }
        /* If {step[ax]} was zero, so is {step[bx]}: */
        if (A->step[ax] != 0) { A->npos *= sz*A->size[bx]; }
      }
    if (sz == 1) { A->step[ax] = 0; }
  }
  
void ppv_diagonal ( ppv_array_t *A, ppv_axis_t ax, ppv_axis_t bx )
  { /* Check axes: */
    if (ax >= ppv_NAX) { ppv_error("bad axis a"); }
    if (bx >= ppv_NAX) { ppv_error("bad axis b"); }
    /* Trivial case: */
    if (ax == bx) { return; }
    /* Eliminate factors {size[ax]} and {size[bx]} from {npos}: */
    if (A->step[ax] != 0) { A->npos /= A->size[ax]; }
    if (A->step[bx] != 0) { A->npos /= A->size[bx]; }
    /* The new size in {ax} is the lesser of the two sizes: */
    if (A->size[bx] < A->size[ax]) 
      { A->size[ax] = A->size[bx]; }
    /* The step along {ax} is the sum of both steps, unless the array is trivial: */
    if (A->size[ax] <= 1)
      { A->step[ax] = 0; } 
    else
      { A->step[ax] += A->step[bx];
        A->npos *= A->size[ax];
      } 
    /* The array becomes trivial in direction {bx}: */
    A->size[bx] = 1; A->step[bx] = 0;
    /* Note that the final array is memoryless if and only if the
      original array was memoryless; so we don't need to check for
      this condition, or change {A.el}. If {A} was empty in the
      direction {bx}, it will become non-empty (size 1) in that
      direction, but will also become empty in the direction {ax}; so,
      that case is not a problem either. */
  }

/* ELEMENT ENUMERATION */

void ppv_enum ( ppv_array_t *A, ppv_array_t *B, ppv_array_t *C, ppv_op_t op )
  {
    /* Obtain a non-null operand {X}: */
    ppv_array_t *X = A;
    if (X == NULL) { X = B; }
    if (X == NULL) { X = C; } 
    /* If all three operands are null, there is nothing to do: */
    if (X == NULL) { return; }
    /* If two or more operands are non-null, their sizes must agree: */
    if ((B != NULL) && (B != X)) 
      { ppv_axis_t ax;
        for (ax = 0; ax < ppv_NAX; ax++)
          if (X->size[0] != B->size[0]) 
            { ppv_error("array size mismatch"); }
      }
    if ((C != NULL) && (C != X)) 
      { ppv_axis_t ax;
        for (ax = 0; ax < ppv_NAX; ax++)
          if (X->size[0] != C->size[0]) 
            { ppv_error("array size mismatch"); }
      }
      
    /* Now loop on samples, first index faster: */
    ppv_index_t ix[ppv_NAX];
    for (ix[5] = 0; ix[5] < X->size[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < X->size[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < X->size[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < X->size[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < X->size[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < X->size[0]; ix[0]++)
                { 
                  /* Get sample positions: */
                  ppv_pos_t posA = (A == NULL ? 0 : ppv_sample_pos(A, ix));
                  ppv_pos_t posB = (B == NULL ? 0 : ppv_sample_pos(B, ix));
                  ppv_pos_t posC = (C == NULL ? 0 : ppv_sample_pos(C, ix));
                  op(ix, posA, posB, posC);
                }
  }

/* ERROR MESSAGES */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf )
  {
    fprintf(stderr, "%s", pf); 
    fprintf(stderr, "bps = %u bpw = %u", A->bps, A->bpw);
    ppv_axis_t ax;
    for (ax = 0; ax < ppv_NAX; ax++)
      { fprintf(stderr, " %u(×%d)", A->size[ax], A->step[ax]); }
    fprintf(stderr, " base = %u npos = %u el = %u", 
      A->base, A->npos, ((unsigned int)A->el)
    );
    fprintf(stderr, "%s", sf);
  }
  
void ppv_progerror 
  ( const char *msg, 
    const char *file, 
    const unsigned int line, 
    const char* proc 
  )
  { fprintf (stderr, "%s:%u: *** (%s) %s\n", file, line, proc, msg);
    exit(1);
  }
