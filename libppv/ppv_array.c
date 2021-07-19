/* See ppv_array.h */
/* Last edited on 2021-07-16 10:14:00 by jstolfi */
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
#include <ix.h>

#include <ppv_array.h>

/* INTERNAL PROOTYPES: */
    
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
    if (A->maxsmp > ppv_max_sample(A->bps))
      { fail_test(die, "invalid max sample value"); }
    if (! ix_parms_are_valid(A->d, A->size, A->base, A->step, die)) { return FALSE; };
    if (! ix_positions_are_distinct(A->d, A->size, A->step, die)) { return FALSE; };
    if ((ix_is_empty(A->d, A->size) || (A->maxsmp == 0)) != (A->el == 0)) 
      { fail_test(die, "inconsistent element storage area"); }
    return TRUE;
  }

ppv_array_t *ppv_array_new ( ppv_dim_t d, ppv_size_t *sz, ppv_sample_t maxsmp )
  { 
    ppv_array_t *A = ppv_array_new_desc(d);
    
    /* Check {maxsmp} and store it: */
    demand(maxsmp <= ppv_MAX_SAMPLE_VAL, "max sample value too big");
    A->maxsmp = maxsmp;

    /* Compute {bps} and store it: */
    ppv_nbits_t bps = ppv_min_bps(maxsmp);
    assert(bps <= ppv_MAX_BPS);
    A->bps = bps;

    /* Check {bpw} and store it: */
    ppv_nbits_t bpw = ppv_best_bpw(bps);
    assert(bpw <= ppv_MAX_BPW);
    A->bpw = bpw; 

    /* The {base} of a new array is 0: */
    A->base = 0;
    
    /* Save {A->size}: */
    for (ppv_axis_t ax = 0; ax < A->d; ax++) { A->size[ax] = sz[ax]; }

    /* Compute number {npos} of distinct positions, set {A->step}: */
    ppv_sample_count_t npos = ppv_compute_npos_steps(d, A->size, A->step);

    /* Set {A->el}: */
    size_t nbytes = ppv_tot_sample_bytes(npos, bps, bpw);
    if (nbytes == 0)
      { /* Memoryless array: */ 
        A->el = NULL;
      }
    else
      { demand(nbytes <= ppv_MAX_BYTES, "too many bytes");
        A->el = notnull(malloc(nbytes), "no mem for new array");
      }
    return A;
  }

ppv_sample_count_t ppv_compute_npos_steps ( ppv_dim_t d, ppv_size_t size[], ppv_step_t step[] )
  { 
    ppv_sample_count_t npos = 1;
    for (ppv_axis_t ax = 0; ax < d; ax++) 
      { ppv_size_t szi = size[ax];
        demand(szi <= ppv_MAX_SIZE, "bad array size");
        if (szi == 0) 
          { npos = 0; break; }
        else if (szi == 1) 
          { if (step != NULL) { step[ax] = 0; } }
        else
          { /* Store position increment along axis {ax}, if relevant: */
            if (step != NULL) { step[ax] = npos; }
            /* Check overflow BEFORE multiplication: */
            if (npos > 0) { demand(szi <= ppv_MAX_SAMPLES/npos, "too many samples"); }
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
    C->maxsmp = A->maxsmp;
    return C;
  }

ppv_sample_count_t ppv_sample_count( ppv_array_t *A, bool_t reptoo )
  {
    ppv_sample_count_t nv;
    if (reptoo)
      { nv = ix_num_tuples(A->d, A->size); }
    else
      { nv = ix_num_positions(A->d, A->size, A->step); }
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
    size_t head_bytes = iroundup(sizeof(ppv_array_t), 8);  /* Mem size wihout {A.size,A.step} vectors. */
    size_t size_bytes = d * sizeof(ppv_size_t);  /* Mem size of {A.size} vector. */
    size_t step_bytes = d * sizeof(ppv_step_t);  /* Mem size of {A.step} vector. */
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
    A->maxsmp = 0;
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

ppv_sample_t ppv_max_sample( ppv_nbits_t bps )
  { 
    ppv_sample_t maxsmp = (ppv_sample_t)((((uint64_t)1) << bps) - 1);
    return maxsmp;
  }

ppv_nbits_t ppv_min_bps( ppv_sample_t maxsmp )
  { ppv_nbits_t bps = 0;
    while (maxsmp > 0) { maxsmp = (maxsmp >> 1); bps++; }
    return bps;
  }

ppv_sample_t ppv_get_sample ( ppv_array_t *A, const ppv_index_t ix[] )
  { demand(ix_is_valid(A->d, ix, A->size), "invalid index");
    ppv_pos_t pos = ppv_sample_pos(A, ix);
    ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
    assert(smp <= A->maxsmp);
    return smp;
  }

void ppv_set_sample ( ppv_array_t *A, const ppv_index_t ix[], ppv_sample_t smp )
  { demand(ix_is_valid(A->d, ix, A->size), "invalid index");
    ppv_pos_t pos = ppv_sample_pos(A, ix);
    demand(smp <= A->maxsmp, "invalid sample value");
    ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
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
        ppv_word_t mask = (ppv_word_t)((((uint64_t)1)<<bps) - 1);   /* Mask to extract low-end sample. */
        ppv_sample_t smp = (qw >> shift) & mask;
        /* fprintf(stderr, "get: pos = " ppv_pos_t_FMT,  pos); */
        /* fprintf(stderr, " el = %016lx  qw = %08x", (uint64_t)el, qw); */
        /* fprintf(stderr, " shift = %d mask = %08x smp = %08x\n", shift, mask, smp); */
        return smp; 
      }
    else if (bps > bpw)
      { int wps = (bps + bpw - 1) / bpw; /* Words per sample. */
        ppv_word_t mask = (ppv_word_t)((((uint64_t)1) << bps) - 1); /* Mask to remove padding. */
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
        ppv_sample_t smp = qw & mask;
        /* fprintf(stderr, " smp = %08x", smp); */
        return smp;
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
    ppv_sample_t smp 
  )
  { ppv_word_t vmask = (ppv_word_t)((((uint64_t)1) << bps) - 1); /* Mask to remove padding. */
    demand(smp <= vmask, "bad pixel value");
    if (bps == 0)
      { return; }
    else if (bps < bpw)
      { ppv_nbits_t spw = bpw/bps; /* Samples per word. */
        ppv_nbits_t shift = (ppv_nbits_t)((spw - 1 - pos%spw)*bps); /* Shift to apply to sample: */
        ppv_word_t shmask = vmask << shift; /* Mask to extract sample. */
        /* Insert the sample into the appropriate word: */
        
        /* fprintf(stderr, "set: shmask = %u ~shmask = %u shift = %d\n", shmask, (~shmask), shift); */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos/spw; 
            (*q8) = (ppv_word_08_t)(((*q8) & (~shmask)) | ((smp << shift) & shmask));
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos/spw;
            (*q16) = (ppv_word_16_t)(((*q16) & (~shmask)) | ((smp << shift) & shmask));
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos/spw;
            /* fprintf(stderr, "set: pos = " ppv_pos_t_FMT, pos); */
            /* fprintf(stderr, " el = %016lx qw = %08x", (uint64_t)el, (*q32)); */
            /* fprintf(stderr, " smp = %08x shmask = %08x shift = %d", smp, shmask, shift); */
            (*q32) = (ppv_word_32_t)(((*q32) & (~shmask)) | ((smp << shift) & shmask));
            /* fprintf(stderr, " qw = %08x\n", (*q32)); */
          }
      }
    else if (bps > bpw)
      { ppv_nbits_t wps = (ppv_nbits_t)((bps + bpw - 1) / bpw); /* Words per sample. */
        ppv_word_t maskw = (ppv_word_t)((((uint64_t)1) << bpw) - 1); /* Mask to chop a wordful from the sample. */
        ppv_nbits_t bhw = (ppv_nbits_t)(((int32_t)bps - 1)%bpw + 1); /* Number of sample bits in the first word. */
        ppv_word_t masks = (ppv_word_t)((((uint64_t)1) << bhw) - 1); /* Mask to chop the bits of in the first word. */
        /* Split the sample into words and store them */
        if (bpw == 8) 
          { ppv_word_08_t *q8 = ((ppv_word_08_t *)el) + pos*wps + wps - 1;
            /* fprintf(stderr, "set: pos = " ppv_pos_t_FMT " el = %016lx", pos, (uint64_t)el); */
            /* fprintf(stderr, " smp = %08x\n", smp); */
            for (int32_t k = 1; k < wps; k++) 
              { /* fprintf(stderr, " q8 = %02x", (*q8)); */
                (*q8) = (ppv_word_08_t)(smp & maskw);
                /* fprintf(stderr, " -> %02x\n", (*q8)); */
                smp >>= bpw;
                q8--; }
            /* fprintf(stderr, " q8 = %02x", (*q8)); */
            (*q8) = (ppv_word_08_t)(((*q8) & (~masks)) | (smp & masks)); 
            /* fprintf(stderr, " -> %02x\n", (*q8)); */
          }
        else if (bpw == 16) 
          { ppv_word_16_t *q16 = ((ppv_word_16_t *)el) + pos*wps + wps -1;
            for (int32_t k = 1; k < wps; k++) 
              { (*q16) = (ppv_word_16_t)(smp & maskw); smp >>= bpw; q16--; }
            (*q16) =  (ppv_word_16_t)(((*q16) & (~masks)) | (smp & masks)); 
          }
        else 
          { ppv_word_32_t *q32 = ((ppv_word_32_t *)el) + pos*wps + wps -1;
            for (int32_t k = 1; k < wps; k++) 
              { (*q32) =  (ppv_word_32_t)(smp & maskw); smp >>= bpw; q32--; }
            (*q32) =  (ppv_word_32_t)(((*q32) & (~masks)) | (smp & masks)); 
          }
      }
    else /* bps == bpw */
      { if (bpw == 8) 
          { *(((ppv_word_08_t *)el) + pos) = (ppv_word_08_t)smp; }
        else if (bpw == 16) 
          { *(((ppv_word_16_t *)el) + pos) = (ppv_word_16_t)smp; }
        else 
          { *(((ppv_word_32_t *)el) + pos) = (ppv_word_32_t)smp; }
      }
  }

void ppv_sample_range(ppv_array_t *A, ppv_sample_t *smp_minP, ppv_sample_t *smp_maxP)
  { 
    ppv_sample_t smp_min = A->maxsmp;
    ppv_sample_t smp_max = 0;
    auto bool_t update(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC);
      /* Updates {smp_min,smp_max} with the sample at {pA}. */
      
    if (A->bps > 0) { ppv_enum(update, FALSE, A, NULL, NULL); }
    (*smp_minP) = smp_min;
    (*smp_maxP) = smp_max;
    return;
    
    bool_t update(const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC)
      { ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pA);
        assert(smp <= A->maxsmp);
        if (smp < smp_min) { smp_min = smp; }
        if (smp > smp_max) { smp_max = smp; }
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
    S->maxsmp = A->maxsmp;
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
    S->maxsmp = A->maxsmp;
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
    auto bool_t copy_voxel( const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
      /* Sets the sample of {A} at position {pA} to the sample value of {B} at positon {pB}.
        Ignores {ix} and {pC}. */
      
    (void) ppv_enum(copy_voxel, FALSE, A, B, NULL);
    return;
    
    /* Internal procedures: */
    bool_t copy_voxel( const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC )
      { ppv_sample_t smp = ppv_get_sample_at_pos(B->el, B->bps, B->bpw, pB);
        assert(smp <= B->maxsmp);
        demand(smp <= A->maxsmp, "invalid sample value");
        ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pA, smp);
        return FALSE;
      }
  }

/* ERROR MESSAGES */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf )
  {
    fprintf(wr, "%s", pf); 
    fprintf(wr, "d = %u", A->d);
    fprintf(wr, " bps = %u bpw = %u", A->bps, A->bpw);
    fprintf(wr, " maxsmp = " ppv_sample_t_FMT, A->maxsmp);
    fprintf(wr, " size(step) = [ ");
    for (ppv_axis_t ax = 0; ax < A->d; ax++)
      { fprintf(wr, " ");
        fprintf(wr, ppv_size_t_FMT, A->size[ax]);
        fprintf(wr, "(");
        fprintf(wr, ppv_step_t_FMT, A->step[ax]);
        fprintf(wr, ")");
      }
    fprintf(wr, " ] base = ");
    fprintf(wr, ppv_pos_t_FMT, A->base);
    fprintf(wr, " el = %lu", ((uint64_t)A->el));
    ppv_sample_count_t nix = ppv_sample_count(A, TRUE);
    fprintf(wr, " nix = " ppv_sample_count_t_FMT, nix);
    ppv_sample_count_t nel = ppv_sample_count(A, FALSE);
    fprintf(wr, " nel = " ppv_sample_count_t_FMT, nel);
    fprintf(wr, "%s", sf);
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

void ppv_dump_storage(FILE *wr, void *el, ppv_pos_t pos, ppv_nbits_t bps, ppv_nbits_t bpw, int32_t nxw)
  {
    if (bps == 0)
      { fprintf(stderr, "  zero-length pixels -- no storage area to dump, sample is zero.\n");
        return;
      }
    if (el == NULL)
      { fprintf(stderr, "  array has no storage area\n");
        return;
      }

    /* Find word count and word and bit indices {iwp,ibt} of element {pos}: */
    int64_t nw;   /* Number of words that contain the sample. */
    int64_t iwp;  /* Index of word that contains the high bit of element. */
    int64_t ibt;  /* Index of high bit of element in that word. */
    if (bps <= bpw)
      { /* One or more whole samples per word: */
        int64_t spw = bpw/bps;             /* Whole samples per word. */
        iwp = ((int64_t)pos)/spw;
        nw = 1;
        int64_t isw = ((int64_t)pos)%spw;  /* Index of sample in word. */
        ibt = ((int64_t)bpw - spw*bps) + isw*bps;
      }
    else
      { /* Two or more words per sample: */
        int64_t wps = (bps + bpw - 1)/bpw;  /* Words per sample. */
        iwp = ((int64_t)pos)*wps;
        nw = wps;
        ibt = wps*bpw - bps;
      }
    assert(iwp >= 0);
    assert(ibt >= 0);
    
    /* Print {nw} words starting with that one, plus {nxw} words on each side for context: */
    int64_t iw_ini = (iwp < nxw ? 0 : iwp - nxw);
    int64_t iw_fin = iwp + nw - 1 + nxw;
    
    int64_t kb_ini = ibt;           /* Index of first sample bit relative to start of the word {iwp}. */
    int64_t kb_fin = ibt + bps - 1; /* Index of last sample bit relative to start of the word {iwp}. */
    
    /* Print word indices: */
    for (int64_t iw = iw_ini; iw <= iw_fin; iw++)
      { if (iw != 0) { fprintf(wr, " "); }
        fprintf(wr, "%*ld", bpw, iw); 
      }
    fprintf(stderr, "\n");
    
    /* Print the bits, compute the sample: */
    int32_t kb; /* Sequential index of bit relative to start of the word {iwp}. */
    uint64_t smp = 0; /* Sample value. */
    for (uint64_t iw = iw_ini; iw <= iw_fin; iw++)
      { if (iw == iwp) { kb = 0; }
        ppv_word_t w;
        if (bpw == 8) 
            { w = *(((ppv_word_08_t *)el) + iw); }
          else if (bpw == 16) 
            { w = *(((ppv_word_16_t *)el) + iw); }
          else 
            { w = *(((ppv_word_32_t *)el) + iw); }
        if (iw != 0) { fprintf(wr, " "); }
        for (int32_t axb = bpw-1; axb >= 0; axb--)
          { int32_t bit = (w >> axb) & 1;
            fprintf(wr, "%d", bit);
            if ((kb >= kb_ini) && (kb <= kb_fin)) 
              { smp = (smp << 1) | bit; }
            kb++;
          }
      }
    fprintf(stderr, "\n");
    
    /* Highlight the bits of the sample: */
    for (uint64_t iw = iw_ini; iw <= iw_fin; iw++)
      { if (iw == iwp) { kb = 0; }
        if (iw != 0) { fprintf(wr, " "); }
        for (int32_t axb = bpw-1; axb >= 0; axb--)
          { if ((kb >= kb_ini) && (kb <= kb_fin)) 
              { fprintf(wr, "^"); }
            else
              { fprintf(wr, " "); }
            kb++;
          }
      }
    fprintf(stderr, "\n");
   
    fprintf(stderr, "sample value = %lu\n", smp);
    fprintf(stderr, "\n");
  }

void ppv_throw_noise(ppv_array_t *A)
  {
    bool_t debug = FALSE;
    ppv_nbits_t bps = A->bps;
    ppv_nbits_t bpw = A->bpw;
    ppv_sample_t maxsmp = A->maxsmp;

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
      { ppv_sample_t smp = (ppv_sample_t)uint64_mrandom(maxsmp);
        ppv_set_sample_at_pos(A->el, bps, bpw, pA, smp);
        if (debug) 
          { if (nset < 256) 
              { fprintf(stderr, " " ppv_sample_t_FMT, smp); }
            else if (nset == 256) 
              { fprintf(stderr, " ..."); }
          }
        nset++;
        return FALSE;
      }
  }
        
void ppv_throw_balls(ppv_array_t *A)
  {
    bool_t debug = FALSE;
    ppv_dim_t d = A->d;
    ppv_nbits_t bps = A->bps;
    ppv_nbits_t bpw = A->bpw;
    ppv_sample_t maxsmp = A->maxsmp;
    if ((maxsmp == 0) || (d == 0)) { /* Single value or single voxel: */ return; }

    /* Choose centers, radii, and values: */
    int32_t nb = 10; /* Number of balls. */
    ppv_sample_count_t npos = ppv_sample_count(A, FALSE);
    double rad_max = 0.2*pow((double)npos, 1.0/d);
    double **ctr = notnull(malloc(nb*sizeof(double *)), "no mem");
    double rad[nb];
    ppv_sample_t bsmp[nb];
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
        bsmp[ib] = (ppv_sample_t)imax(uint64_abrandom(1, maxsmp),uint64_abrandom(1, maxsmp));
        if (debug) { fprintf(stderr, " value = " ppv_sample_t_FMT "\n", bsmp[ib]); }
      }
    auto bool_t set_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
    ppv_enum(set_voxel, FALSE, A, NULL, NULL);
    return;
    
    bool_t set_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { ppv_sample_t smp = 0; /* Value to store. */
        /* Check if voxel center is inside some ball: */
        for(int32_t ib = 0; ib < nb; ib++)
          { double d2 = 0.0; /* Squared dist from center of ball {b} */
            for (ppv_axis_t ax = 0; ax < A->d; ax++) 
              { double da = ((double)ix[ax] + 0.5) - ctr[ib][ax]; 
                d2 += da*da;
              }
            if (sqrt(d2) < rad[ib]) { smp = bsmp[ib]; }
          }
        ppv_set_sample_at_pos(A->el, bps, bpw, pA, smp);
        return FALSE;
      }
  }
