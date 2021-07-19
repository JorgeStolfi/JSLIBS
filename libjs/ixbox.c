/* See {indexing_box.h}. */
/* Last edited on 2021-07-17 08:27:06 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ix.h>
#include <ix_io.h>
#include <jsmath.h>
#include <affirm.h>

#include <ixbox.h>

#define ixFMT   ix_index_t_FMT
#define smpFMT  ix_sample_t_FMT
#define szFMT   ix_size_t_FMT

bool_t ixbox_is_empty(ix_dim_t d, ix_size_t  size[])
  { if (size == NULL) { return TRUE; }
    for (ix_axis_t k = 0; k < d; k++)
      { if (size[k] == 0) { return TRUE; } }
    return FALSE;
  }

bool_t ixbox_has(ix_dim_t d, ix_index_t ix[], ix_index_t ixlo[], ix_size_t size[])
  {
    if (size == NULL) { return FALSE; }
    for (ix_axis_t k = 0; k < d;  k++)
      { ix_index_t dxk = ix[k] - (ixlo == NULL ? 0 : ixlo[k]);
        if ((dxk < 0) || (dxk >= size[k])) { return FALSE; }
      }
    return TRUE;
  }

void ixbox_intersect
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[],
    ix_index_t ixlo_R[], 
    ix_size_t  size_R[]
  )
  {
    bool_t empty_A = ixbox_is_empty(d, size_A);
    bool_t empty_B = ixbox_is_empty(d, size_B);
    bool_t empty_R;
    if (empty_A || empty_B)
      { empty_R = TRUE; }
    else
      { /* Compute the new core domain {ixlo_R,size_R} and {empty_R}: */
        for (ix_axis_t k = 0; k < d;  k++)
          { /* Get the index range  {ixlok_A..ixhik_A} of box{A} on axis {k}: */
            ix_index_t ixlok_A = (ix_index_t)(ixlo_A == NULL? 0 : ixlo_A[k]);
            ix_index_t ixhik_A = (ix_index_t)(ixlok_A + (size_A == NULL ? 0 : size_A[k]) - 1);

            /* Get the index range {ixlok_B..ixhik_B} of  box {B} on axis {k}: */
            ix_index_t ixlok_B = (ix_index_t)(ixlo_B == NULL? 0 : ixlo_B[k]);
            ix_index_t ixhik_B = (ix_index_t)(ixlok_B + (size_B == NULL ? 0 : size_B[k]) - 1);

            /* Get the expected clipped core box {ixlo_R[k],size_R[k]}: */
            ix_index_t ixlok_R = (ix_index_t)imax(ixlok_A, ixlok_B);
            ix_index_t ixhik_R = (ix_index_t)imin(ixhik_A, ixhik_B);
            if (ixlok_R > ixhik_R) { empty_R = TRUE; break; }
            ixlo_R[k] = ixlok_R; 
            size_R[k] = (ix_size_t)(ixhik_R - ixlok_R + 1);
          }
       }
    if (empty_R) { for (ix_axis_t k = 0; k < d;  k++) { ixlo_R[k] = 0; size_R[k] = 0; } }
    return;
  }
  
bool_t ixbox_is_contained
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[]
  )
  {
    bool_t empty_A = ixbox_is_empty(d, size_A);
    if (empty_A) { return TRUE; }
    
    bool_t empty_B = ixbox_is_empty(d, size_B);
    if (empty_B) { return FALSE; }
    
    for (ix_axis_t k = 0; k < d;  k++)
      { /* Get the index range  {ixlok_A..ixhik_A} of box{A} on axis {k}: */
        ix_index_t ixlok_A = (ix_index_t)(ixlo_A == NULL? 0 : ixlo_A[k]);
        ix_index_t ixhik_A = (ix_index_t)(ixlok_A + (size_A == NULL ? 0 : size_A[k]) - 1);
        assert(ixlok_A <= ixhik_A);

        /* Get the index range {ixlok_B..ixhik_B} of  box {B} on axis {k}: */
        ix_index_t ixlok_B = (ix_index_t)(ixlo_B == NULL? 0 : ixlo_B[k]);
        ix_index_t ixhik_B = (ix_index_t)(ixlok_B + (size_B == NULL ? 0 : size_B[k]) - 1);
        assert(ixlok_A <= ixhik_A);
        
        if ((ixlok_A < ixlok_B) || (ixhik_A > ixhik_B)) { return FALSE; }
      }
    return TRUE;
  }
        
bool_t ixbox_equal
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[]
  )
  {
    bool_t empty_A = ixbox_is_empty(d, size_A);
    bool_t empty_B = ixbox_is_empty(d, size_B);
    if (empty_A && empty_B) { return TRUE; }
    if (empty_A || empty_B) { return FALSE; }
    /* Compare the two ranges aong each axis: */
    assert((size_A != NULL) && (size_B != NULL));
    for (ix_axis_t k = 0; k < d;  k++)
      { if (size_A[k] != size_B[k]) { return FALSE; }
        ix_index_t ixlok_A = (ix_index_t)(ixlo_A == NULL? 0 : ixlo_A[k]);
        ix_index_t ixlok_B = (ix_index_t)(ixlo_B == NULL? 0 : ixlo_B[k]);
        if (ixlok_A != ixlok_B) { return FALSE; }
      }
    return TRUE;
  }
  
void ixbox_print
  ( FILE *wr, 
    ix_dim_t d, 
    char *pref, 
    ix_index_t ixlo[], 
    ix_size_t size[], 
    char *suff
  )
  { 
    if (pref != NULL) { fputs(pref, wr); }
    for (ix_axis_t k = 0; k <d; k++)
      { if (k > 0) { fputs(" ", wr); }
        fprintf(wr, ixFMT "(" szFMT ")", ixlo[k], size[k]);
      }
    if (suff != NULL) {fputs(suff, wr); }
    return;
  }

