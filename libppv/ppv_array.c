/* See ppv_array.h */
/* Last edited on 2021-07-03 15:04:58 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <indexing.h>

#include <ppv_array.h>

/* INTERNAL PROOTYPES: */
    
ppv_array_t *ppv_array_new_desc ( ppv_dim_t d );
  /* Allocates a new descriptor {A} for an array with dimension {d}. The
    {A.step} and {A.size} vectors are allocated in the same {malloc}
    record and set to all zeros, as {A.base}. The element area pointer
    {A.el} is set to {NULL}. */

/* IMPLEMENTATIONS: */

bool_t ppv_index_is_valid ( const ppv_index_t ix[], ppv_array_t *A )
  { return ix_is_valid(A->d, ix, A->size); }

void ppv_index_clear ( ppv_dim_t d, ppv_index_t ix[] )
  { ix_fill(d, ix, 0); }

void ppv_index_assign ( ppv_dim_t d, ppv_index_t ix[], const ppv_index_t val[] )
  { ix_assign(d, ix, val); }

void ppv_index_shift ( ppv_dim_t d, ppv_index_t ix[], ppv_index_t *inc )
  { ix_shift(d, ix, inc); } 

sign_t ppv_index_compare ( ppv_dim_t d, const ppv_index_t ixa[], const ppv_index_t ixb[] )
  { return ix_compare(d, ixa, ixb, ix_order_L); }
    
bool_t ppv_index_first ( ppv_index_t ix[], ppv_array_t *A )
  { return ix_assign_min(A->d, ix, A->size); } 

bool_t ppv_index_last ( ppv_index_t ix[], ppv_array_t *A )
  { return ix_assign_max(A->d, ix, A->size); } 

bool_t ppv_index_next ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p )
  { demand(d <= A->d, "d too big"); 
    return ix_next(d, ix, A->size, ix_order_L, A->step, p, NULL, NULL, NULL, NULL);
  }

bool_t ppv_index_prev ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p )
  { demand(d <= A->d, "d too big"); 
    return ix_prev(d, ix, A->size, ix_order_L, A->step, p, NULL, NULL, NULL, NULL);
  }

bool_t ppv_descriptor_is_valid ( ppv_array_t *A, bool_t die )
  { assert(ppv_MAX_DIM == ix_MAX_DIM); 
    assert(ppv_MAX_POS == ix_MAX_POS); 
    assert(ppv_MAX_ABS_STEP == ix_MAX_ABS_STEP);
    assert(ppv_MAX_SIZE == ix_MAX_SIZE);
    if (! ix_parms_are_valid(A->d, A->size, A->base, A->step, die)) { return FALSE; };
    if (! ix_positions_are_distinct(A->d, A->size, A->step, die)) { return FALSE; };
    return TRUE;
  }

ppv_array_t *ppv_array_new ( ppv_dim_t d, ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw )
  { 
    ppv_array_t *A = ppv_array_new_desc(d);
    
    /* Check {bps} and store it: */
    demand(bps <= ppv_MAX_BPS, "bits-per-sample too big");
    A->bps = bps;
    /* Check {bpw} and store it: */
    demand(bpw <= ppv_MAX_BPW, "bits-per-word too big");
    A->bpw = bpw; 
    /* The {base} is always 0: */
    A->base = 0;
    
    /* Save {A->size}, see if array is empty: */
    for (ppv_axis_t ax = 0; ax < A->d; ax++) { A->size[ax] = sz[ax]; }
    /* Compute number {npos} of distinct positions, set {A->step}: */
    ppv_sample_count_t npos = ppv_compute_npos_steps(d, A->size, A->step);

    size_t nbytes = ppv_tot_sample_bytes(npos, bps, bpw);
    /* Set {A->step,A->el}: */
    if (nbytes == 0)
      { /* Memoryless array: */ 
        A->el = NULL;
      }
    else
      { demand (nbytes <= ppv_MAX_BYTES, "too many bytes");
        A->el = notnull(malloc(nbytes), "no mem for new array");
      }
    return A;
  }

ppv_sample_count_t ppv_compute_npos_steps ( ppv_dim_t d, ppv_size_t size[], ppv_step_t step[] )
  { 
    ppv_sample_count_t npos = 1;
    for (ppv_axis_t ax = 0; ax < d; ax++) 
      { ppv_size_t szi = size[ax];
        demand (szi <= ppv_MAX_SIZE, "bad array size");
        if (szi == 0) 
          { npos = 0; break; }
        else if (szi == 1) 
          { if (step != NULL) { step[ax] = 0; } }
        else
          { /* Store position increment along axis {ax}, if relevant: */
            if (step != NULL) { step[ax] = npos; }
            /* Check overflow BEFORE multiplication: */
            if (npos > 0) { demand (szi <= ppv_MAX_SAMPLES/npos, "too many samples"); }
            npos *= szi;
          }
      }
    
    /* If the array is empty, all steps should be zero: */
    if ((npos == 0) && (step != NULL))
      { for (ppv_axis_t ax = 0; ax < d; ax++) { step[ax] = 0; } }
      
    return npos;
  }
 
ppv_array_t *ppv_array_clone ( ppv_array_t *A )
  { 
    ppv_dim_t d = A->d;
    ppv_array_t *C = ppv_array_new_desc(d);
    for (ppv_axis_t ax = 0; ax < d; ax++) { C->step[ax] = A->step[ax]; C->size[ax] = A->size[ax]; }
    C->base = A->base;
    C->el = A->el;
    C->bps = A->bps;
    C->bpw = A->bpw;
    return C;
  }

ppv_sample_count_t ppv_sample_count( ppv_array_t *A, bool_t reptoo )
  {
    ppv_sample_count_t nv = 1;
    for (ppv_axis_t ax = 0; ax < A->d; ax++)
      { ppv_step_t stepk = A->step[ax];
        ppv_size_t sizek = A->size[ax];
        if ((sizek >= 2) && (stepk == 0) && (! reptoo)) { sizek = 1; }
        if (sizek == 0) 
          { nv = 0; break; }
        else
          { assert(nv <= ppv_MAX_SIZE/sizek);
            nv = nv*sizek;
          }
      }
    assert((A->el == NULL) == ((nv == 0) || (A->bps == 0)));
    return nv;
  }

bool_t ppv_is_empty( ppv_array_t *A )
  { for (ppv_axis_t ax = 0; ax < A->d; ax++)
      { if (A->size[ax] == 0) { return TRUE; } }
    return FALSE;
  }

ppv_array_t *ppv_array_new_desc ( ppv_dim_t d )
  { 
    size_t head_bytes = sizeof(ppv_array_t);   /* Mem size wihout {A.size,A.step} vectors. */
    size_t size_bytes = d * sizeof(ppv_size_t);     /* Mem size of {A.size} vector. */
    size_t step_bytes = d * sizeof(ppv_step_t);     /* Mem size of {A.step} vector. */
    size_t desc_bytes = head_bytes + size_bytes + step_bytes;
    
    ppv_array_t *A = notnull(malloc(desc_bytes), "no mem");
    A->size = (ppv_size_t *)(((char*)A) + head_bytes);
    A->step = (ppv_step_t *)(((char*)A->size) + size_bytes);
    
    /* The following check assumes {A.el} is the last field of the descriptor.
      The placement of {A.el} in the final position should ensure that
      the addresses of {A->size} and {A->step} are synchronized to 64 bits. */
    assert(((char*)A->size) == ((char*)&(A->el)) + sizeof(void*));
    
    /* Clear the fields, just in case: */
    A->d = d;
    for (ppv_axis_t ax = 0; ax < d; ax++) { A->size[ax] = 0; A->step[ax] = 0; }
    A->base = 0;
    A->bps = 0;
    A->bpw = 0;
    A->el = NULL;
    
    return A;
  }

ppv_nbits_t ppv_best_bpw( ppv_nbits_t bps )
  { 
    ppv_nbits_t bpw = ppv_MAX_BPW; /* Default. */
    /* Cases when {bpw} is best equal to {bps}: */
    if (bps == 0)
      { /* Irrelevant: */ bpw = 8; }
    else if ((bps == 8) || (bps == 16) || (bps == 32))
      { bpw = bps; }
    else if ((bps % 32) == 0)
      { bpw = 32; }
    else if ((bps % 16) == 0)
      { bpw = 16; }
    else if ((bps % 8) == 0)
      { bpw = 8; }
    else
      { double waste_min = +INF; /* Best fraction of space wasted. */
        int32_t wps_min;         /* Words per sample at that {bpw}. */
        for (ppv_nbits_t b = 8; b <= 32; b = (ppv_nbits_t)(2*b))
           { int32_t wps = (bps + b - 1)/b; /* Count of {b}-bit words per sample. */
             int32_t spw = (bps >= b ? 1 : (wps*b)/bps); /* Samples that fit in those words. */
             double waste = ((double)(wps*b - spw*bps))/((double)(wps*b)); /* Wastage.*/
             if ((waste < waste_min) || ((waste == waste_min) && (wps < wps_min)))
               { bpw = b; waste_min = waste; wps_min = wps; }
           }
      }
    return bpw;
  }

ppv_sample_t ppv_get_sample ( ppv_array_t *A, const ppv_index_t ix[] )
  { ppv_pos_t pos = ppv_sample_pos(A, ix);
    ppv_sample_t qv = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
    return qv;
  }

void ppv_set_sample ( ppv_array_t *A, const ppv_index_t ix[], ppv_sample_t qv )
  { ppv_pos_t pos = ppv_sample_pos(A, ix);
    ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, qv);
    return;
  }

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, const ppv_index_t ix[] )
  { return ix_position(A->d, ix, A->base, A->step); }

ppv_sample_t ppv_get_sample_at_pos 
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
        ppv_sample_t qv = (qw >> shift) & mask;
        /* fprintf(stderr, "get: pos = " ppv_pos_t_FMT,  pos); */
        /* fprintf(stderr, " el = %016lx  qw = %08x", (uint64_t)el, qw); */
        /* fprintf(stderr, " shift = %d mask = %08x qv = %08x\n", shift, mask, qv); */
        return qv; 
      }
    else if (bps > bpw)
      { int wps = (bps + bpw - 1) / bpw; /* Words per sample. */
        ppv_word_t mask = ((ppv_word_t)1 << bps) - 1; /* Mask to remove padding. */
        /* Gather the sample words into {qw}: */
        ppv_word_t qw = 0;
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos*wps;
            /* fprintf(stderr, "get: pos = " ppv_pos_t_FMT " el = %016lx", pos, (uint64_t)el); */
            /* fprintf(stderr, " q8 ="); */
            for (int32_t k = 0; k < wps; k++) 
              { /* fprintf(stderr, " %02x", (*q8)); */
                qw = (qw << bpw) | (*q8);
                /* fprintf(stderr, " -> %08x", qw); */
                q8++;
              }
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos*wps;
            for (int32_t k = 0; k < wps; k++) { qw = (qw << bpw) | (*q16); q16++; }
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos*wps;
            for (int32_t k = 0; k < wps; k++) { qw = (qw << bpw) | (*q32); q32++; }
          }
        ppv_sample_t qv = qw & mask;
        /* fprintf(stderr, " qv = %08x", qv); */
        return qv;
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

void ppv_set_sample_at_pos 
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
        
        /* fprintf(stderr, "set: mask = %u ~mask = %u shift = %d\n", mask, (~mask), shift); */
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
            /* fprintf(stderr, "set: pos = " ppv_pos_t_FMT, pos); */
            /* fprintf(stderr, " el = %016lx qw = %08x", (uint64_t)el, (*q32)); */
            /* fprintf(stderr, " qv = %08x mask = %08x shift = %d", qv, mask, shift); */
            (*q32) = (ppv_word_32_t)(((*q32) & (~mask)) | ((qv << shift) & mask));
            /* fprintf(stderr, " qw = %08x\n", (*q32)); */
          }
      }
    else if (bps > bpw)
      { ppv_nbits_t wps = (ppv_nbits_t)((bps + bpw - 1) / bpw); /* Words per sample. */
        ppv_word_t maskw = (1 << bpw) - 1; /* Mask to chop a wordful from the sample. */
        ppv_word_t masks = (1 << (bps - wps*bpw)) - 1; /* Mask to protect padding. */
        /* Split the sample into words and store them */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos*wps + wps - 1;
            /* fprintf(stderr, "set: pos = " ppv_pos_t_FMT " el = %016lx", pos, (uint64_t)el); */
            /* fprintf(stderr, " qv = %08x", qv); */
            for (int32_t k = 1; k < wps; k++) 
              { /* fprintf(stderr, " q8 = %02x", (*q8)); */
                (*q8) = (ppv_word_08_t)(qv & maskw);
                /* fprintf(stderr, " -> %02x", (*q8)); */
                qv >>= bpw;
                q8--; }
            /* fprintf(stderr, " q8 = %02x", (*q8)); */
            (*q8) = (ppv_word_08_t)(((*q8) & (~masks)) | (qv & masks)); 
            /* fprintf(stderr, " -> %02x\n", (*q8)); */
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos*wps + wps -1;
            for (int32_t k = 1; k < wps; k++) 
              { (*q16) = (ppv_word_16_t)(qv & maskw); qv >>= bpw; q16--; }
            (*q16) =  (ppv_word_16_t)(((*q16) & (~masks)) | (qv & masks)); 
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos*wps + wps -1;
            for (int32_t k = 1; k < wps; k++) 
              { (*q32) =  (ppv_word_32_t)(qv & maskw); qv >>= bpw; q32--; }
            (*q32) =  (ppv_word_32_t)(((*q32) & (~masks)) | (qv & masks)); 
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

void ppv_sample_range(ppv_array_t *A, ppv_sample_t *vminP, ppv_sample_t *vmaxP)
  { 
    ppv_sample_t vmin = (ppv_sample_t)((1 << A->bps) - 1);
    ppv_sample_t vmax = 0;
    auto bool_t update(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC);
      /* Updates {vmin,vmax} with the sample at {pA}. */
      
    if (A->bps > 0) { ppv_enum(update, FALSE, A, NULL, NULL); }
    (*vminP) = vmin;
    (*vmaxP) = vmax;
    return;
    
    bool_t update(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC)
      { ppv_sample_t val = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        if (val < vmin) { vmin = val; }
        if (val > vmax) { vmax = val; }
        return FALSE;
      }
  }



/* DESCRIPTOR MANIPULATION */

void ppv_crop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t skip, ppv_size_t keep )
  { ix_crop(A->d, A->size, &(A->base), A->step, ax, skip, keep);
    /* If the array has become empty, {A.el} is now useless: */
    if (keep == 0) { A->el = NULL; }
  }

void ppv_subsample ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t stride )
  { ix_subsample(A->d, A->size, &(A->base), A->step, ax, stride); }

void ppv_reverse ( ppv_array_t *A, ppv_axis_t ax )
  { ix_flip(A->d, A->size, &(A->base), A->step, ax); }

void ppv_replicate ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz )
  { ix_replicate(A->d, A->size, &(A->base), A->step, ax, sz); }

void ppv_swap_indices ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1, ppv_dim_t n )
  { ix_swap_indices(A->d, A->size, &(A->base), A->step, ax0, ax1, n); }

void ppv_flip_indices ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1 )
  { ix_flip_indices(A->d, A->size, &(A->base), A->step, ax0, ax1); }

void ppv_diagonal ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1 )
  { ix_diagonal(A->d, A->size, &(A->base), A->step, ax0, ax1); }

ppv_array_t *ppv_slice ( ppv_array_t *A, ppv_axis_t ax, ppv_index_t ix )
  { ppv_dim_t d = A->d;
    demand(ax < d, "invalid axis");
    
    /* Make a local copy of {A}'s fields: */
    ppv_pos_t baseA = A->base;
    ppv_size_t sizeA[d];
    ppv_step_t stepA[d];
    for(ppv_axis_t ax = 0; ax < d; ax++) { sizeA[ax] = A->size[ax]; stepA[ax] = A->step[ax]; } 
    
    /* Apply the slicing operation to them: */
    ix_axis_t axv[1];  axv[0] = ax;
    ix_index_t ixv[1]; ixv[0] = ix;
    ix_slice(d, sizeA, &baseA, stepA, 1, axv, ixv);
    assert(sizeA[d-1] == 1);
    assert(stepA[d-1] == 0);
    
    /* Copy the parameters in a new descriptor, with one dimension less: */
    ppv_dim_t dS = (ppv_dim_t)(d-1);
    ppv_array_t *S = ppv_array_new_desc(dS);
    S->base = baseA;
    for(ppv_axis_t ax = 0; ax < dS; ax++) { S->size[ax] = sizeA[ax]; S->step[ax] = stepA[ax]; }
    S->bps = A->bps;
    S->bpw = A->bpw;
    S->el = A->el;
    
    return S;
  }

ppv_array_t *ppv_chop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz )
  { 
    ppv_dim_t d = A->d;
    demand(d < ppv_MAX_DIM, "too many axes already");
    demand(ax < d, "invalid axis");
    demand((sz > 0) && (sz < ppv_MAX_SIZE), "invalid segment size");
    
    /* Make a local copy of {A}'s fields, with one dimension more: */
    ppv_dim_t dS =(ppv_dim_t)(d+1);
    ppv_pos_t baseA = A->base;
    ppv_size_t sizeA[dS];
    ppv_step_t stepA[dS];
    for(ppv_axis_t ax = 0; ax < d; ax++) { sizeA[ax] = A->size[ax]; stepA[ax] = A->step[ax]; } 
    sizeA[d] = 1;
    stepA[d] = 0;
    
    /* Apply the chopping operation to them: */
    ix_chop(dS, sizeA, &baseA, stepA, ax, sz, d);
    assert(sizeA[ax] == sz);
    assert(sizeA[d] == A->size[ax]/sz);
    
    /* Copy the parameters in a new descriptor, with one dimension more: */
    ppv_array_t *S = ppv_array_new_desc(dS);
    S->base = baseA;
    for(ppv_axis_t ax = 0; ax < dS; ax++) { S->size[ax] = sizeA[ax]; S->step[ax] = stepA[ax]; }
    S->bps = A->bps;
    S->bpw = A->bpw;
    S->el = A->el;
    
    return S;
  }
  
/* ELEMENT ENUMERATION */

bool_t ppv_enum 
  ( ppv_index_pos3_op_t *op, 
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
    if (X == NULL) { return FALSE; }
    /* If two or more operands are non-null, their sizes must agree: */
    if ((B != NULL) && (B != X)) 
      { demand(ix_same_size(A->d, B->size, X->size, FALSE), "B array size mismatch"); }
    if ((C != NULL) && (C != X)) 
      { demand(ix_same_size(A->d, C->size, X->size, FALSE), "C array size mismatch"); }
      
    /* Get the bases and steps of non-null arguments: */
    ppv_pos_t bA = (A == NULL ? 0 : A->base);
    ppv_step_t *sA = (A == NULL ? NULL : A->step);
    
    ppv_pos_t bB = (B == NULL ? 0 : B->base);
    ppv_step_t *sB = (B == NULL ? NULL : B->step);
    
    ppv_pos_t bC = (C == NULL ? 0 : C->base);
    ppv_step_t *sC = (C == NULL ? NULL : C->step);
    
    return ix_enum(op, A->d, X->size, ix_order_L, reverse, bA, sA, bB, sB, bC, sC);
  }

void ppv_array_assign ( ppv_array_t *A, ppv_array_t *B  )
  {
    ppv_sample_t mask = (ppv_sample_t)((((ppv_sample_t)1) << A->bps) - 1); /* {A} sample mask. */
    
    auto bool_t copy_voxel( const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
      /* Sets the sample of {A} at position {pA} to the sample value of {B} at positon {pB}.
        Ignores {ix} and {pC}. */
      
    (void) ppv_enum(copy_voxel, FALSE, A, B, NULL);
    return;
    
    /* Internal procedures: */
    bool_t copy_voxel( const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC )
      { ppv_sample_t vsmp = ppv_get_sample_at_pos(B->el, B->bps, B->bpw, pB);
        ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, vsmp & mask);
        return FALSE;
      }
  }

/* ERROR MESSAGES */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf )
  {
    fprintf(stderr, "%s", pf); 
    fprintf(stderr, "bps = %u bpw = %u", A->bps, A->bpw);
    for (ppv_axis_t ax = 0; ax < A->d; ax++)
      { fprintf(stderr, " ");
        fprintf(stderr, ppv_size_t_FMT, A->size[ax]);
        fprintf(stderr, "(×");
        fprintf(stderr, ppv_step_t_FMT, A->step[ax]);
        fprintf(stderr, ")");
      }
    fprintf(stderr, " base = ");
    fprintf(stderr, ppv_pos_t_FMT, A->base);
    fprintf(stderr, " el = %lu", ((uint64_t)A->el));
    fprintf(stderr, "%s", sf);
  }

size_t ppv_tot_sample_bytes(ppv_sample_count_t npos, ppv_nbits_t bps, ppv_nbits_t bpw)
  {
    /* Compiler/machine sanity check: */
    assert(sizeof(ppv_word_08_t) == 1);
    assert(sizeof(ppv_word_16_t) == 2);
    assert(sizeof(ppv_word_32_t) == 4);
    demand((bpw == 8) || (bpw == 16) || (bpw == 32), "unsupported bits-per-word");
    
    int32_t Bpw = bpw/8; /* Bytes per word. */

    /* Compute number of words {nwords}: */
    size_t nwords; 
    if (bps == 0)
      { /* Zero-length samples - no storage required: */
        nwords = 0;
      }
    else if (bps < bpw)
      { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
        nwords = (npos + spw - 1) / spw;
        demand( nwords < ppv_MAX_BYTES/Bpw, "too many bytes in array");
      }
    else if (bps > bpw)
      { ppv_nbits_t wps = (ppv_nbits_t)(bps + bpw - 1)/bpw; /* Words per sample. */
        demand(npos < ppv_MAX_BYTES/(wps*Bpw), "too many bytes in array");
        nwords = npos * wps;
      }
    else
      { demand(npos < ppv_MAX_BYTES/Bpw, "too many bytes in array");
        nwords = npos;
      }
    /* Allocate storage are of correct size: */
    size_t nbytes = nwords * Bpw;
    return nbytes;
  }

void ppv_choose_test_size(ppv_dim_t d, ppv_sample_count_t npos, ppv_size_t sz[])
  { demand(npos <= ppv_MAX_SAMPLES, "too many samples");
    if (d == 0)
      { /* Nothing to do. */}
    else if (d == 1) 
      { sz[0] = (ppv_size_t)imax(1, imin(npos, ppv_MAX_SIZE)); }
    else
      { double xsize_avg = (npos == 0 ? 0 : pow((double)npos, 1.0/d)); /* Geom avg size. */
        double xfmax = sqrt(3);
        double xfmin = 1.0/xfmax;
        for (ppv_axis_t ax = 0; ax < d; ax++)
          { double xf = xfmin * pow(xfmax/xfmin, ((double)ax)/(d-1)); 
            sz[ax] = (ppv_size_t)imax(1, imin((int64_t)ceil(xf*xsize_avg), ppv_MAX_SIZE)); 
          }
      }
  }

void ppv_throw_noise(ppv_array_t *A)
  {
    bool_t debug = TRUE;
    ppv_nbits_t bps = A->bps;
    ppv_nbits_t bpw = A->bpw;
    ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1);

    auto bool_t throw_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_sample_count_t nset = 0;  /* Count of samples assigned. */
    if (debug) { fprintf(stderr, "samples = "); }
    ppv_enum(throw_voxel, FALSE, A, NULL, NULL);
    
    if (debug) { fprintf(stderr, "\n"); }
    if (debug) 
      { fprintf(stderr, ppv_sample_count_t_FMT " samples set\n", nset);
        ppv_sample_count_t npos = ppv_sample_count(A, TRUE); 
        assert(nset == npos);
      }
    return;
    
    bool_t throw_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t val = (ppv_sample_t)uint64_mrandom(maxval);
        ppv_set_sample_at_pos(A->el, bps, bpw, pA, val);
        if (debug) 
          { if (nset < 256) 
              { fprintf(stderr, " " ppv_sample_t_FMT, val); }
            else if (nset == 256) 
              { fprintf(stderr, " ..."); }
          }
        nset++;
        return FALSE;
      }
  }
        
void ppv_throw_balls(ppv_array_t *A)
  {
    bool_t debug = TRUE;
    ppv_dim_t d = A->d;
    ppv_nbits_t bps = A->bps;
    ppv_nbits_t bpw = A->bpw;
    if ((bps == 0) || (d == 0)) { /* Single value or single voxel: */ return; }
    ppv_sample_t maxval = (ppv_sample_t)((1 << bps) - 1);

    /* Choose centers and radii: */
    int32_t nb = 10; /* Number of balls. */
    ppv_sample_count_t npos = ppv_sample_count(A, FALSE);
    double rad_max = 0.2*pow((double)npos, 1.0/d);
    double **ctr = notnull(malloc(nb*sizeof(double *)), "no mem");
    double rad[nb];
    ppv_sample_t val[nb];
    if (debug) { fprintf(stderr, "throwing %d balls...\n", nb); }
    for (int32_t ib = 0; ib < nb; ib++)
      { if (debug) { fprintf(stderr, "  ball %3d center = (", ib); }
        ctr[ib] = notnull(malloc(d*sizeof(double)), "no mem");
        for (ppv_axis_t ax = 0; ax < A->d; ax++) 
          { ctr[ib][ax] = dabrandom(0, (double)A->size[ax]); 
            if (debug) { fprintf(stderr, " %8.2f", ctr[ib][ax]); }
          }
        rad[ib] = dabrandom(0.3*rad_max, rad_max);
        if (debug) { fprintf(stderr, " ) radius = %8.2f", rad[ib]); }
        val[ib] = (ppv_sample_t)uint64_abrandom(1, maxval);
        if (debug) { fprintf(stderr, " value = " ppv_sample_t_FMT "\n", val[ib]); }
      }
    auto bool_t set_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(set_voxel, FALSE, A, NULL, NULL);
    return;
    
    bool_t set_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t v = 0; /* Value to store. */
        /* Check if voxel center is inside some ball: */
        for(int32_t ib = 0; ib < nb; ib++)
          { double d2 = 0.0; /* Squared dist from center of ball {b} */
            for (ppv_axis_t ax = 0; ax < A->d; ax++) 
              { double da = ((double)ix[ax] + 0.5) - ctr[ib][ax]; 
                d2 += da*da;
              }
            if (sqrt(d2) < rad[ib]) { v = val[ib]; }
          }
        ppv_set_sample_at_pos(A->el, bps, bpw, pA, v);
        return FALSE;
      }
    
  
  }
