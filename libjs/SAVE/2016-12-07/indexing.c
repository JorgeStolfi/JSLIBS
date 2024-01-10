/* See indexing.h */
/* Last edited on 2013-10-25 01:20:52 by stolfilocal */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <indexing.h>
#include <affirm.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

/* IMPLEMENTATIONS */

ix_pos_t ix_position ( ix_dim_t d, const ix_index_t ix[], ix_pos_t bp, const ix_step_t st[] )
  { ix_pos_t p = bp;
    /* Should check range {0..ix_MAX_POS} */
    int i; 
    for (i = 0; i < d; i++)
      { p += (ix_pos_t)(ix[i]*st[i]);
        /* Should check range {0..ix_MAX_POS} */
      }
    return p;
  }

ix_pos_t ix_position_safe ( ix_dim_t d, const ix_index_t ix[], const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[] )
  { if (ix_is_valid(d, ix, sz))
      { return ix_position(d, ix, bp, st); }
    else
      { return ix_pos_NONE; }
  }

ix_pos_t ix_packed_position ( ix_dim_t d, const ix_index_t ix[], const ix_size_t sz[], ix_order_t ixor )
  { ix_pos_t p = 0;
    if (d == 0) { return p; }
    ix_step_t st = 1;
    int k = 0;
    while (TRUE)
      { /* Get the axis {i}, which is axis number {k} in order {ixor}: */
        int i = (ixor == ix_order_F ? k : d - 1 - k);
        /* Check validity of index {ix[i]}: */
        if ((ix[i] < 0) || (ix[i] >= sz[i])) { return ix_pos_NONE; }
        /* Update position {p} to include index {ix[i]}: */
        p += ix[i]*st;
        /* !!! Should check range {0..ix_MAX_POS} !!! */
        k++;
        if (k >= d) { break; }
        /* !!! Should worry about exceeding {ix_MAX_ABS_STEP} !!! */
        st *= (ix_step_t)sz[i];
      }
    return p;
  }

void ix_packed_steps ( ix_dim_t d, const ix_size_t sz[], ix_order_t ixor, ix_step_t st[] )
  {
    /* First pass to check limits and whether array is empty: */
    bool_t empty = (d == 0);
    ix_axis_t i, k;
    for (i = 0; i < d; i++) 
      { ix_size_t szi = sz[i]; 
        demand (szi <= ix_MAX_SIZE, "bad array size");
        if (szi == 0) { empty = TRUE; }
      }

    /* Set the position steps {st[0..NA-1]}: */
    if (empty)
      { /* Array has no elements; reset all steps to zero. */ 
        for (i = 0; i < d; i++) { st[i] = 0; }
      }
    else
      { ix_pos_count_t npos = 1;
        for (k = 0; k < d; k++) 
          { /* Get the axis {i}, which is axis number {k} in order {ixor}: */
            i = (ix_axis_t)(ixor == ix_order_F ? k : d - 1 - k);
            ix_size_t szi = sz[i]; 
            if (szi == 1) 
              { /* Trivial axis -- increment is irrelevant: */
                st[i] = 0;
              }
            else
              { /* Non-trivial axis -- increment is current total size: */
                st[i] = npos;
                /* Check overflow BEFORE multiplication: */
                if (npos > 0) { demand (szi <= ix_MAX_POSITIONS/npos, "too many samples"); }
                npos *= szi;
              }
          }
      }
  }

void ix_packed_indices ( ix_dim_t d, ix_pos_t p, const ix_size_t sz[], ix_order_t ixor, ix_index_t ix[] )
  { /* Should check range {0..ix_MAX_POS} */
    /* assert(p >= 0); */
    int k; 
    for (k = 0; k < d; k++) 
      { /* Get the axis {i}, which is axis number {k} in the order {ixor}: */
        int i = (ixor == ix_order_F ? k : d - 1 - k);
        /* Should check {sz[i]>0} */
        ix_size_t w = sz[i]; 
        ix[i] = p % w; 
        /* Should check {ix[i] < ix_MAX_INDEX} */
        p /= w;
      }
    assert(p == 0);
  }

void ix_steps_assign ( ix_dim_t d, ix_step_t sta[], const ix_step_t stb[] )
  { int i;
    for (i = 0; i < d; i++) { sta[i] = stb[i]; }
  }

/* INDEX TUPLE MANIPULATION */

void ix_fill ( ix_dim_t d, ix_index_t ix[], ix_index_t val )
  { 
    int i; 
    for (i = 0; i < d; i++) { ix[i] = val; }
  }

void ix_assign ( ix_dim_t d, ix_index_t ix[], const ix_index_t val[] )
  { 
    int i; 
    for (i = 0; i < d; i++) 
      { /* Should check {val[i]>0} */ ix[i] = val[i]; }
  }

bool_t ix_assign_min ( ix_dim_t d, ix_index_t ix[], const ix_size_t sz[] )
  { 
    int i; 
    bool_t valid = TRUE;
    for (i = 0; i < d; i++) 
      { ix[i] = 0; if (sz[i] == 0) { valid = FALSE; } }
    return valid;
  }
  
bool_t ix_assign_max ( ix_dim_t d, ix_index_t ix[], const ix_size_t sz[] )
  { 
    int i; 
    bool_t valid = TRUE;
    for (i = 0; i < d; i++) 
      { ix[i] = sz[i]-1; if (sz[i] == 0) { valid = FALSE; } }
    return valid;
  }

void ix_shift ( ix_dim_t d, ix_index_t ix[], const ix_index_t inc[] )
  { 
    int i; 
    for (i = 0; i < d; i++) { /* Should check {inc[i]>-ix[i]} */  ix[i] += inc[i]; }
  }

bool_t ix_is_valid ( ix_dim_t d, const ix_index_t ix[], const ix_size_t sz[] )
  { 
    int i; 
    for (i = 0; i < d; i++) 
      { if ((ix[i] < 0) || (ix[i] >= sz[i])) { return FALSE; } }
    return TRUE;
  }

/* SIZE TUPLE MANIPULATION */

void ix_sizes_assign ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] )
  { int i;
    for (i = 0; i < d; i++) { sza[i] = szb[i]; }
  }

bool_t ix_sizes_shrink ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] )
  { int i;
    bool_t changed = FALSE;
    for (i = 0; i < d; i++)
      { if (sza[i] > szb[i]) { sza[i] = szb[i]; changed = TRUE; } }
    return changed;
  }

bool_t ix_sizes_expand ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] )
  { int i;
    bool_t changed = FALSE;
    for (i = 0; i < d; i++)
      { if (sza[i] < szb[i]) { sza[i] = szb[i]; changed = TRUE; } }
    return changed;
  }

ix_size_t ix_max_size ( ix_dim_t d, const ix_size_t sz[] )
  { int i;
    ix_size_t msz = 0;
    for (i = 0; i < d; i++) { if (sz[i] > msz) { msz = sz[i]; } }
    return msz;
  }

ix_size_t ix_min_size ( ix_dim_t d, const ix_size_t sz[] )
  { int i;
    ix_size_t msz = ix_MAX_SIZE;
    for (i = 0; i < d; i++) { if (sz[i] < msz) { msz = sz[i]; } }
    return msz;
  }

/* ELEMENT COUNTING */

bool_t ix_is_empty ( ix_dim_t d, const ix_size_t sz[] )
  { int i; 
    for (i = 0; i < d; i++) 
      { if (sz[i] <= 0) { return TRUE; } }
    return FALSE;
  }

ix_pos_count_t ix_num_tuples ( ix_dim_t d, const ix_size_t sz[] )
  { uint64_t n = 1;
    int i;
    for (i = 0; i < d; i++) 
      { ix_size_t szi = sz[i]; 
        if (szi < 1) { return 0; }
        if (szi > 1)
          { demand(szi <= ix_MAX_POSITIONS / n, "too many tuples to count");
            n *= szi;
          }
      }
    return n;
  }
  
ix_pos_count_t ix_num_positions ( ix_dim_t d, const ix_size_t sz[], const ix_step_t st[] )
  { uint64_t n = 1;
    int i;
    for (i = 0; i < d; i++) 
      { ix_size_t szi = sz[i]; 
        if (szi < 1) { return 0; }
        if ((szi > 1) && (st[i] != 0))
          { demand (szi <= ix_MAX_POSITIONS / n, "overflow");
            n *= szi;
          }
      }
    return n;
  }

ix_pos_t ix_min_pos ( ix_dim_t d, const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[] )
  { ix_pos_t p = bp;
    int i;
    for (i = 0; i < d; i++) 
      { ix_step_t sti = st[i];
        if (sti < 0)
          { ix_size_t szi = sz[i]; 
            if (szi <= 0) { return ix_pos_NONE; }
            p += (szi-1)*sti;
          }
      }
    return p;
  }
  
ix_pos_t ix_max_pos ( ix_dim_t d, const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[] )
  { ix_pos_t p = bp;
    int i;
    for (i = 0; i < d; i++) 
      { ix_step_t sti = st[i];
        if (sti > 0)
          { ix_size_t szi = sz[i]; 
            if (szi <= 0) { return ix_pos_NONE; }
            p += (szi-1)*sti;
          }
      }
    return p;
  }

bool_t ix_same_size ( ix_dim_t d, const ix_size_t sza[], const ix_size_t szb[], bool_t die )
  { ix_axis_t i;
    for (i = 0; i < d; i++) 
      { if (sza[i] != szb[i]) { fail_test(die,"bp too big"); } }
      return TRUE;
  }

bool_t ix_contained ( ix_dim_t d, const ix_size_t sza[], const ix_size_t szb[], bool_t die )
  { ix_axis_t i;
    for (i = 0; i < d; i++) 
      { if (sza[i] > szb[i]) { fail_test(die,"bp too big"); } }
      return TRUE;
  }

/* STEP MANIPULATION */

void ix_crop 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t skip, 
    ix_size_t keep
  )
  { /* Check axis: */
    demand(i < d, "bad axis"); 
    /* Check and store the new size: */
    demand(skip + keep <= sz[i], "bad range");
    sz[i] = keep;
    if (keep == 0)
      { /* Array becomes empty: */
        (*bp) = 0;
        ix_axis_t k;
        for (k = 0; k < d; k++) { st[k] = 0; }
      }
    else
      { /* Update {A->step[ax]}: */
        if (st[i] != 0)
          { if (skip > 0)
              { /* Update {*bp} (can't overflow because element exists): */
                (*bp) += skip * st[i];
              }
            /* If index is now trivial, clear step: */
            if (keep == 1) { st[i] = 0; }
          } 
      } 
  }

void ix_subsample 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t stride 
  )
  { /* Check axis: */
    demand(i < d, "bad axis");
    /* Check stride: */
    demand(stride > 0, "bad stride");
    /* Check for trivial subsampling: */
    if (stride == 1) { return; }
    /* Save the original size along axis {i}: */
    ix_size_t osz = sz[i];
    /* If array is empty or trivial along that index, there is nothing to do: */
    if (osz <= 1) { return; }
    /* Compute and store the new size: */
    ix_size_t nsz = (osz-1)/stride + 1;
    assert(nsz > 0);
    sz[i] = nsz;
    /* Adjust step if needed: */
    if (st[i] != 0)
      { if (nsz == 1) 
          { /* Index {i} became trivial: */
            st[i] = 0;
          }
        else
          { /* Index {i} now is {stride} times as powerful: */
            st[i] *= stride;
          }
      }
  }

void ix_flip 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i
  )
  { /* Check axis: */
    demand(i < d, "bad axis");
    /* If there is no motion along that axis, we have nothing to d: */
    ix_step_t sti = st[i];
    if (sti == 0) { return; }
    /* Now {sz[i]} must be positive. */
    /* Adjust {*bp} to the last valid index value: */
    (*bp) += (sz[i]-1)*sti; 
    /* Negate the step {st[i]}: */
    st[i] = -sti;
  }

void ix_swap_indices 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_axis_t j ,
    ix_dim_t n
  )
  { /* Check axes: */
    demand(i+n <= d, "bad axis i");
    demand(j+n <= d, "bad axis j");
    /* Trivial transpose: */
    if ((i == j) || (n == 0)) { return; }
    demand((i >= j+n) || (j>= i+n), "overlapping axis blocks");  
    /* Exchange parameters of axes {i,j}: */
    while(n > 0)
      { { ix_size_t tmp = sz[i]; sz[i] = sz[j]; sz[j] = tmp; }
        { ix_pos_t  tmp = st[i]; st[i] = st[j]; st[j] = tmp; }
        i++; j++; n--;
      }
  }

void ix_flip_indices
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_axis_t j 
  )
  { /* Check axes: */
    demand(i < d, "bad axis i");
    demand(j < d, "bad axis j");
    /* Swap symmetric entries: */
    while (i < j)
      { /* Exchange parameters of axes {i,j}: */
        { ix_size_t tmp = sz[i]; sz[i] = sz[j]; sz[j] = tmp; }
        { ix_pos_t tmp = st[i]; st[i] = st[j]; st[j] = tmp; } 
        i++; j--;
      }
  }

void ix_slice
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_dim_t n, 
    const ix_axis_t ax[], 
    const ix_index_t ix[] 
  )
  { if (n == 0) { return; }
    /* Crop the selected indices to the specified values: */
    { int k; 
      for (k = 0; k < n; k++) 
        { ix_axis_t i = (ix_axis_t)k;
          if (ax != NULL)
            { i = ax[k];
              demand((i+0 >= 0) && (i < d), "invalid axis");
              if (k > 0) { demand(i > ax[k-1], "axes out of order"); }
            }
          ix_crop(d, sz, bp, st, i, ix[k], 1);
        }
    }
    { /* Pull the remaining indices to the front: */ 
      int k = 0;
      ix_axis_t j = 0;
      ix_axis_t i;
      for (i = 0; i < d; i++) 
        { if ((k < n) && ((ax == NULL) || (i == ax[k]))) 
            { k++; }
          else
            { sz[j] = sz[i]; st[j] = st[i]; j++; } 
        }
      /* Make the remaining indices trivial: */
      while(j < d) { sz[j] = 1; st[j] = 0; j++; }
    }
  }

void ix_replicate 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t size 
  )
  { /* Check axis: */
    demand(i < d, "bad axis");
    /* Axis {i} must be trivial: */
    demand(sz[i] == 1, "index is not trivial");
    /* Check and set the new size: */
    demand((size > 0) && (size <= ix_MAX_SIZE), "bad new size");
    sz[i] = size;
    /* The increment must already be zero, but just in case... */
    st[i] = 0;
  }

void ix_diagonal 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i,
    ix_axis_t j 
  )
  { /* Check axes: */
    demand(i < d, "bad axis i");
    demand(j < d, "bad axis j");
    demand(i != j, "axes are equal");
    /* Require the array to be wide enough: */
    demand(sz[j] >= sz[i], "array has wrong shape");
    /* The step along {i} is the sum of both steps: */
    if (sz[i] > 1) { st[i] += st[j]; }
    /* The size along {j} is reduced by {sz[i]-1}: */
    if (sz[i] == 0)
      { /* Increment {sz[j]} for consistency, but beware fo overflow: */
        demand(sz[j] <= ix_MAX_SIZE - 1, "array too big on axis j");
        sz[j]++;
      }
    else
      { sz[j] = sz[j] - (sz[i] - 1); }
    /* May become trivial (but not empty) along {j}: */
    if (sz[j] <= 1) { st[j] = 0; }
  }

void ix_chop 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t stride, 
    ix_axis_t j 
  )
  { /* Check axes: */
    demand(i < d, "bad axis i");
    demand(j < d, "bad axis j");
    demand(i != j, "axes must be different");
    /* The chunk size must be positive: */
    demand(stride > 0, "bad chunk size");
    /* The array must be trivial in direction {j}: */
    demand(sz[j] == 1, "axis j not trivial");
    /* Save the original size along axis {i}: */
    ix_size_t osz = sz[i];
    /* The size must be a multiple of the chunk size: */
    demand(osz % stride == 0, "non-integral number of chunks");
    /* Compute new sizes along {j} and {i}: */
    sz[j] = osz/stride;  sz[i] = stride;
    /* Compute new steps along {j} and {i}: */
    if (sz[j] == 1)
      { /* Array is still trivial along {j}; just in case: */
        st[j] = 0;
      }
    else
      { /* This computation cannot overflow.: */
        st[j] = stride * st[i];
      }
    /* Check whether axis {i} became trivial: */
    if (stride == 1) { st[i] = 0; } 
  }
 
/* ELEMENT ENUMERATION */

bool_t ix_next
  ( ix_dim_t d, 
    ix_index_t ix[],
    const ix_size_t sz[],
    ix_order_t ixor, 
    const ix_step_t stA[],
    ix_pos_t *pA, 
    const ix_step_t stB[],
    ix_pos_t *pB, 
    const ix_step_t stC[],
    ix_pos_t *pC 
  )
  { /* Simplify the situation: */
    if (stA == NULL) { pA = NULL; }
    if (stB == NULL) { pB = NULL; }
    if (stC == NULL) { pC = NULL; }
    /* Look for an index {ix[i]} that can be incremented: */
    int k;
    for (k = 0; k < d; k++)
      { /* Get the axis {i}, which is axis number {k} in the order {ixor}: */
        int i = (ixor == ix_order_F ? k : d - 1 - k);
        /* Grab the index {ix[i]}: */
        ix_index_t *ixi = &(ix[i]);
        if ((*ixi) < sz[i]-1) 
          { /* Increment index {ix[i]} */
            (*ixi)++; 
            if (pA != NULL) { (*pA) += stA[i]; }
            if (pB != NULL) { (*pB) += stB[i]; }
            if (pC != NULL) { (*pC) += stC[i]; }
            return FALSE;
          }
        /* Index {ix[i]} is at its limit, reset to 0: */
        if (pA != NULL) { (*pA) -= (*ixi)*stA[i]; }
        if (pB != NULL) { (*pB) -= (*ixi)*stB[i]; }
        if (pC != NULL) { (*pC) -= (*ixi)*stC[i]; }
        (*ixi) = 0;
      }
    return TRUE;
  }

bool_t ix_prev
  ( ix_dim_t d, 
    ix_index_t ix[],
    const ix_size_t sz[],
    ix_order_t ixor, 
    const ix_step_t stA[],
    ix_pos_t *pA, 
    const ix_step_t stB[],
    ix_pos_t *pB, 
    const ix_step_t stC[],
    ix_pos_t *pC
  )
  { /* Simplify the situation: */
    if (stA == NULL) { pA = NULL; }
    if (stB == NULL) { pB = NULL; }
    if (stC == NULL) { pC = NULL; }
    /* Look for an index {ix[i]} that can be decremented: */
    int k;
    for (k = 0; k < d; k++)
      { /* Get the axis {i}, which is axis number {k} in the order {ixor}: */
        int i = (ixor == ix_order_F ? k : d - 1 - k);
        /* Grab the index {ix[i]}: */
        ix_index_t *ixi = &(ix[i]);
        if ((*ixi) > 0)
          { /* Decrement index {ix[i]} */
            (*ixi)--;
            if (pA != NULL) { (*pA) -= stA[i]; }
            if (pB != NULL) { (*pB) -= stB[i]; }
            if (pC != NULL) { (*pC) -= stC[i]; }
            return FALSE;
          }
        /* Index {ix[i]} is at its limit, reset to {sz[i]-1}: */
        (*ixi) = sz[i]-1;
        if (pA != NULL) { (*pA) += (*ixi)*stA[i]; }
        if (pB != NULL) { (*pB) += (*ixi)*stB[i]; }
        if (pC != NULL) { (*pC) += (*ixi)*stC[i]; }
      }
    return TRUE; 
  }
  
sign_t ix_compare ( ix_dim_t d, const ix_index_t ixa[], const ix_index_t ixb[], ix_order_t ixor )
  { int k;
    for (k = 0; k < d; k++)
      { /* Get the axis {i}, which is axis number {k} in the REVERSE of order {ixor}: */
        int i = (ixor == ix_order_F ? d - 1 - k : k);
        const ix_index_t *ixai = ixa+i;
        const ix_index_t *ixbi = ixb+i;
        if ((*ixai) != (*ixbi))
          { if ((*ixai) < (*ixbi))
              { return NEG; }
            else
              { return POS; }
          }
      }
    return ZER;
  }
    
int ix_sync_level ( ix_dim_t d, const ix_size_t sz[], const ix_index_t ix[], ix_order_t ixor, bool_t reverse )
  {
    int k = 0;
    while (k < d)
      { int i = (ixor == ix_order_F ? k : d - 1 - k);
        if (ix[i] != (reverse ? sz[i] - 1 : 0)) { break; }
        k++;
      }
    return k;
  }

void ix_enum 
  ( ix_op_t *op,
    ix_dim_t d,
    const ix_size_t sz[], 
    ix_order_t ixor, 
    bool_t reverse,
    ix_pos_t bpA, 
    const ix_step_t stA[], 
    ix_pos_t bpB, 
    const ix_step_t stB[], 
    ix_pos_t bpC, 
    const ix_step_t stC[]
  )
  { ix_index_t ix[d];
    if (! reverse)
      { if (ix_assign_min(d,ix,sz))
          { ix_pos_t pA = bpA;
            ix_pos_t pB = bpB;
            ix_pos_t pC = bpC;
            do 
              { op(ix, pA,pB,pC); } 
            while (! ix_next(d,ix,sz,ixor,stA,&pA,stB,&pB,stC,&pC));
         }
      }
    else
      { if (ix_assign_max(d,ix,sz))
         { ix_pos_t pA = (stA == NULL ? bpA : ix_position(d,ix,bpA,stA));
           ix_pos_t pB = (stB == NULL ? bpB : ix_position(d,ix,bpB,stB));
           ix_pos_t pC = (stC == NULL ? bpC : ix_position(d,ix,bpC,stC));
           do 
             { op(ix, pA,pB,pC); } 
           while (! ix_prev(d,ix,sz,ixor,stA,&pA,stB,&pB,stC,&pC));
         }
      }
  }

/* PARAMETER VALIDATION */

bool_t ix_parms_are_valid ( ix_dim_t d, const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[], bool_t die )
  { /* Check the position of element {[0,0,...,0]}: */
    if (bp > ix_MAX_POS) { fail_test(die,"bp too big"); }
    
    /* Check sizes and number of (supposedly) distinct positions: */
    ix_dim_t dw = 0; /* Number of non-trivial indices. */
    ix_pos_count_t npos = 1; /* Number of index tuples, excl. replications. */
    ix_axis_t i;
    for (i = 0; i < d; i++) 
      { ix_step_t sti = st[i];
        ix_size_t szi = sz[i];
        /* Count number of axes with non-trivial steps: */
        if (sti != 0) { dw++; }
        /* Validate size and step: */
        if (szi == 0) 
          { /* Array is empty: */ npos = 0; }
        else if (szi >= 2)
          { /* Check size range: */
            if (szi > ix_MAX_SIZE) { fail_test(die,"array size too big"); }
            if ((sti > 0) && (npos > 0))
              { /* Axis is not replicated; more positions (worry about overflow): */
                if (szi > ix_MAX_POSITIONS / npos) { fail_test(die,"too many positions"); }
                npos *= szi;
              }
          }
      }
    if (npos == 0)
      { /* Empty array -- all steps and the base must be zero: */
        if (bp != 0)  { fail_test(die,"nonzero base in empty array"); }
        if (dw != 0)  { fail_test(die,"nonzero step in empty array"); }
        return TRUE;
      }
    /* Check number of (supposedly) distinct positions: */
    if (npos > ix_MAX_POSITIONS) { fail_test(die,"too many positions"); }

    /* Check steps and min/max position */
    ix_pos_t lop = bp; /* Minimum position. */
    ix_pos_t phi = bp; /* Maximum position. */
    for (i = 0; i < d; i++) 
      { ix_step_t sti = st[i];
        ix_size_t szi = sz[i];
        /* Validade size and step: */
        assert(szi > 0); /* Else we should have returned already. */
        if (szi == 1)
          { /* Axis is trivial, step must be zero: */
            if (sti != 0) { fail_test(die,"nonzero step in trivial axis"); }
          }
        else if (sti > 0)
          { /* Check step position: */
            if (sti > ix_MAX_ABS_STEP) { fail_test(die,"step too big"); }
            /* Check that the maximum position is within range: */
            if ((ix_MAX_POS - phi)/(szi - 1) < sti) { fail_test(die,"position too big"); }
            phi += (szi - 1)*sti;
          }
        else if (sti < 0)
          { /* Beware that we can't make {sti} positive before checking: */
            if (sti < -ix_MAX_ABS_STEP) { fail_test(die,"step too big"); }
            /* Make step positive: */
            sti = -sti;
            /* Check that minimum position is non-negative: */
            if (lop/(szi - 1) < sti) { fail_test(die,"negative position"); }
            lop -= (szi - 1)*sti;
          }
        else
          { /* Step is zero; axis is replicated, ignore: */ }
      }
    /* Well, everything seems allright... */
    return TRUE; 
  }

bool_t ix_positions_are_distinct ( ix_dim_t d, const ix_size_t sz[], const ix_step_t st[], bool_t die )
  { /* Check {sz,st}: */

    /* The injectivity is hard to check exactly, so this procedure
    checks first whether there is a sequence of transpositions and
    flips that makes the lexicographical order of the index tuples
    (ignoring replicated axes) imply strict ordering of the positions.
    
    More precisely, let {stw[0..dw-1]} be the absolute values of the
    nonzero steps, sorted in increasing order, and {szw[0..dw-1]} the
    corresponding sizes. A sufficient condiiton for injectivity is that
      { stw[i] > SUM { stw[j]*(szw[j]-1) : j \in 0..i-1 } }
    for all {i} in {0..dw-1}. (This means that the displacement due to
    an increment of 1 in the index corresponding to {stw[i]} exceeds
    the maximum displacement that can be obtained with all the
    previous indices combined.) */

    ix_step_t stw[d]; /* Nonzero steps in increasing order */
    ix_size_t szw[d]; /* Corresponding sizes in the same order */
    ix_dim_t dw = 0; /* Number of non-trivial indices. */
    ix_axis_t i;
    for (i = 0; i < d; i++) 
      { ix_step_t sti = st[i];
        ix_size_t szi = sz[i];
        if (sti != 0)
          { /* Assume the sizes are consistent with the steps: */
            assert(szi >= 2); 
            /* Assume {sti} is not too big to negate: */
            if (sti < 0) { sti = -sti; }
            /* Save {szi,sti} in increasing order: */
            int j = dw;
            while ((j > 0) && (stw[j-1] > sti))
              { stw[j] = stw[j-1];
                szw[j] = szw[j-1];
                j--;
              }
            stw[j] = sti;
            szw[j] = szi;
            dw++;
          }
      }
    /* Compute {tsw[i]} = max displ. using indices {0..i}. */
    ix_pos_t tdw[d]; 
    ix_pos_t tdisp = 0; /*  */
    for (i = 0; i < dw; i++)
      { tdisp += stw[i]*(szw[i]-1);
        tdw[i] = tdisp;
      }
    /* Exclude trailing indices which are strictly monotonic: */
    while ((dw > 0) && (stw[dw-1] > tdw[dw-1])) { dw--; }
    /* If no indices are left, the positions are strictly monotonic: */ 
    if (dw == 0) { return TRUE; }
   
    /* Positions are not motonic, rats. Must do the hard test. */
    /* There should be a smart way to do it. We use brute force: */
    /* we count the number {ct[p]} of index tuples per position {p}. */
    ix_pos_count_t npos = tdw[dw-1] + 1;
    uint8_t *ct = (uint8_t *)malloc(npos*sizeof(uint8_t)); 
    affirm(ct != NULL, "no memory for reachability counters");
    ix_pos_t p;
    /* Clear all counters: */
    for (p = 0; p < npos; p++) { ct[p] = 0; }
    /* Enumerate index tuples and bump counters: */
    ix_index_t ixw[dw];
    for (i = 0; i < dw; i++) { ixw[i] = 0; }
    p = 0;
    bool_t collision = FALSE;
    do 
      { /* assert(p >= 0); */
        assert(p < npos);
        if (ct[p] > 0) { collision = TRUE; break; }
        ct[p]++;
      } 
    while(! ix_next(dw,ixw,szw,ix_order_F,stw,&p,NULL,NULL,NULL,NULL));
    /* Free storage before returning anything: */
    free(ct); 
    /* Verdict: */
    if (! collision) { return TRUE; }
    fail_test(die,"positions are not distinct");
  }

ix_index_t ix_reduce ( ix_index_t i, ix_size_t N, ix_reduction_t red )
  {
    /* Beware that {s % u} with {s} signed and {u} unsigned is {((unsiged)s) % u}. */
    ix_index_t NN = N;
    switch(red)
      {
        case ix_reduction_SINGLE:  /* ... *,*,*,0,1,2,3,4,5,*,*,*,*,*,*,*,*,*,* ... */
          if ((i < 0) || (i >= NN)) { i = -1; }
          break;

        case ix_reduction_EXTEND:  /* ... 0,0,0,0,1,2,3,4,5,5,5,5,5,5,5,5,5,5,5 ... */
          if (i < 0) { i = 0; }
          if (i >= NN) { i = NN-1; }
          break;
        
        case ix_reduction_REPEAT:  /* ... 3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3 ... */
          i %= NN;
          if (i < 0) { i += NN; }
          break;

        case ix_reduction_MIRROR:  /* ... 2,1,0,0,1,2,3,4,5,5,4,3,2,1,0,0,1,2,3 ... */
          i %= (2*NN);
          if (i < 0) { i += 2*NN; }
          if (i >= NN) { i = NN - 1 - (i - NN); }
          break;
          
        case ix_reduction_PXMIRR:  /* ... 3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,4,5 ... */
          if (NN == 1) 
            { i = 0; }
          else
            { i %= (2*NN - 2);
              if (i < 0) { i += 2*NN - 2; }
              if (i >= N) { i = N - 2 - (i - NN); }
            }
          break;

        default:
          demand(FALSE, "invalid reduction code");
      }
    return i;
  }

void ix_reduce_range ( ix_index_t i0, ix_size_t m, ix_size_t N, ix_reduction_t red, ix_index_t i[] )
  { 
    int k;
    for (k = 0; k < m; k++) { i[k] = ix_reduce(i0+k, N, red); }
  }
