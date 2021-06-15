/* See {indexing_descr.h}. */
/* Last edited on 2021-06-13 11:32:52 by jstolfi */

#define indexing_descr_C_COPYRIGHT "Copyright © 2003  Jorge Stolfi, State University of Campinas"

#include <indexing_descr.h>
#include <indexing.h>
#include <indexing_io.h>
#include <bool.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <affirm.h>
    
#define NA ix_descr_NAXES
  /* For short. */

ix_descr_t ix_descr_make ( ix_dim_t na, const ix_size_t sz[], const ix_step_t st[], ix_pos_t bp )
  { ix_descr_t D;
    D.na = na;
    ix_sizes_assign(na, D.sz, sz);
    ix_steps_assign(na, D.st, st);
    D.bp = bp;
    return D;
  }

ix_descr_t ix_descr_from_sizes ( ix_dim_t na, const ix_size_t sz[] )
  { ix_descr_t D;
    D.na = na;
    ix_sizes_assign(na, D.sz, sz);
    ix_packed_steps(na, D.sz, ix_order_L, D.st);
    D.bp = 0;
    return D;
  }

void ix_descr_steps_assign ( ix_descr_t *D, ix_step_t sta[], const ix_step_t stb[] )
  { ix_steps_assign(D->na, sta, stb); }

ix_pos_t ix_descr_position ( ix_descr_t *D, const ix_index_t ix[] )
  { return ix_position(D->na, ix, D->bp, D->st); }

void ix_descr_indices_fill ( ix_descr_t *D, ix_index_t ix[], ix_index_t val )
  { ix_fill(D->na, ix, val); }

void ix_descr_indices_assign ( ix_descr_t *D, ix_index_t ix[], const ix_index_t val[] )
  { ix_assign(D->na, ix, val); }

bool_t ix_descr_indices_first ( ix_descr_t *D, ix_index_t ix[] )
  { return ix_assign_min(D->na, ix, D->sz); }

bool_t ix_descr_indices_last ( ix_descr_t *D, ix_index_t ix[] )
  { return ix_assign_max(D->na, ix, D->sz); }

void ix_descr_indices_shift ( ix_descr_t *D, ix_index_t ix[], const ix_index_t inc[] )
  { ix_shift(D->na, ix, inc); }

bool_t ix_descr_indices_are_valid ( ix_descr_t *D, const ix_index_t ix[] )
  { return ix_is_valid(D->na, ix, D->sz); }

void ix_descr_sizes_assign ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] )
  { ix_sizes_assign(D->na, sza, szb); }

bool_t ix_descr_sizes_shrink ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] )
  { return ix_sizes_shrink(D->na, sza, szb); }

bool_t ix_descr_sizes_expand ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] )
  { return ix_sizes_expand(D->na, sza, szb); }

ix_size_t ix_descr_max_size ( ix_descr_t *D )
  { return ix_max_size(D->na, D->sz); }

ix_size_t ix_descr_min_size ( ix_descr_t *D )
  { return ix_min_size(D->na, D->sz); }

bool_t ix_descr_is_empty ( ix_descr_t *D )
  { return ix_is_empty(D->na, D->sz); }

ix_count_t ix_descr_num_tuples ( ix_descr_t *D )
  { return ix_num_tuples(D->na, D->sz); }

ix_count_t ix_descr_num_positions ( ix_descr_t *D )
  { return ix_num_positions(D->na, D->sz, D->st); }

ix_pos_t ix_descr_min_pos ( ix_descr_t *D )
  { return ix_min_pos(D->na, D->sz, D->bp, D->st); }

ix_pos_t ix_descr_max_pos ( ix_descr_t *D )
  { return ix_max_pos(D->na, D->sz, D->bp, D->st); }

bool_t ix_descr_same_size ( ix_descr_t *A, ix_descr_t *B, bool_t die )
  { 
    if (A->na != B->na) { fail_test(die,"different number of indices"); }
    return ix_same_size(A->na, A->sz, B->sz, die);
  }

bool_t ix_descr_contained ( ix_descr_t *A, ix_descr_t *B, bool_t die )
  { 
    if (A->na != B->na) { fail_test(die,"different number of indices"); }
    return ix_contained(A->na, A->sz, B->sz, die);
  }

bool_t ix_descr_next ( ix_descr_t *D, ix_index_t ix[], ix_pos_t *p )
  { return ix_next(D->na, ix, D->sz, ix_order_L, D->st, p, NULL, NULL, NULL, NULL); }

bool_t ix_descr_prev ( ix_descr_t *D, ix_index_t ix[], ix_pos_t *p )
  { return ix_prev(D->na, ix, D->sz, ix_order_L, D->st, p, NULL, NULL, NULL, NULL); }

sign_t ix_descr_compare ( ix_descr_t *D, const ix_index_t ixa[], const ix_index_t ixb[] )
  { return ix_compare(D->na, ixa, ixb, ix_order_L); }
    
/* DESCRIPTOR MANIPULATION */

void ix_descr_crop ( ix_descr_t *D, ix_axis_t i, ix_size_t skip, ix_size_t keep )
  { ix_crop(D->na, D->sz, &(D->bp), D->st, i, skip, keep); }

void ix_descr_subsample ( ix_descr_t *D, ix_axis_t i, ix_size_t stride )
  { ix_subsample(D->na, D->sz, &(D->bp), D->st, i, stride); }

void ix_descr_flip ( ix_descr_t *D, ix_axis_t i )
  { ix_flip(D->na, D->sz, &(D->bp), D->st, i); }

void ix_descr_replicate ( ix_descr_t *D, ix_axis_t i, ix_size_t sz )
  { ix_replicate(D->na, D->sz, &(D->bp), D->st, i, sz); }

void ix_descr_swap_indices ( ix_descr_t *D, ix_axis_t i, ix_axis_t j, ix_dim_t n )
  { ix_swap_indices(D->na, D->sz, &(D->bp), D->st, i, j, n); }

void ix_descr_flip_indices ( ix_descr_t *D, ix_axis_t i, ix_axis_t j )
  { ix_flip_indices(D->na, D->sz, &(D->bp), D->st, i, j); }

void ix_descr_slice( ix_descr_t *D, ix_dim_t n, const ix_axis_t ax[], const ix_index_t ix[] )
  { ix_slice( D->na, D->sz, &(D->bp), D->st, n, ax, ix );
    D->na = (ix_dim_t)(D->na - n);
  }

void ix_descr_diagonal ( ix_descr_t *D, ix_axis_t i, ix_axis_t j )
  { ix_diagonal(D->na, D->sz, &(D->bp), D->st, i, j); }

void ix_descr_chop ( ix_descr_t *D, ix_axis_t i, ix_size_t sz, ix_axis_t j )
  { ix_chop(D->na, D->sz, &(D->bp), D->st, i, sz, j); }

/* ELEMENT ENUMERATION */

bool_t ix_descr_enum 
  ( ix_descr_index_pos3_op_t *op,
    ix_order_t ixor,
    bool_t reverse,
    ix_descr_t *A,
    ix_descr_t *B,
    ix_descr_t *C
  )
  {
    /* Obtain a non-null operand {X}: */
    ix_descr_t *X = A;
    if (X == NULL) { X = B; }
    if (X == NULL) { X = C; } 

    /* If all three operands are null, there is nothing to do: */
    if (X == NULL) { return; }

    /* Get the effective number of indices: */
    ix_dim_t na = A->na;
    demand(B->na == na, "unequal dimensions");
    demand(C->na == na, "unequal dimensions");
    
    /* Get the intersection {sz} of the domains: */
    ix_size_t sz[na];
    ix_descr_sizes_assign(A, sz, X->sz);
    if ((B != NULL) && (B != X)) { ix_descr_sizes_shrink(B, sz, B->sz); }
    if ((C != NULL) && (C != X)) { ix_descr_sizes_shrink(C, sz, C->sz); }
      
    /* Get the bases and steps of non-null arguments: */
    ix_pos_t bpA = (A == NULL ? 0 : A->bp);
    ix_step_t *stA = (A == NULL ? NULL : A->st);
    
    ix_pos_t bpB = (B == NULL ? 0 : B->bp);
    ix_step_t *stB = (B == NULL ? NULL : B->st);
    
    ix_pos_t bpC = (C == NULL ? 0 : C->bp);
    ix_step_t *stC = (C == NULL ? NULL : C->st);
    
    return ix_enum(op, na, sz, ixor, reverse, bpA, stA, bpB, stB, bpC, stC);
  }

bool_t ix_descr_is_valid ( ix_descr_t *D, bool_t die )
  { if (! ix_parms_are_valid(D->na, D->sz, D->bp, D->st, die)) { return FALSE; };
    return TRUE;
  }

bool_t ix_descr_positions_are_distinct ( ix_descr_t *D, bool_t die )
  { if (! ix_positions_are_distinct(D->na, D->sz, D->st, die)) { return FALSE; };
    return TRUE;
  }

void ix_descr_print_descr ( FILE *wr, char *pre, ix_descr_t *D, int wd, char *suf )
  { /* Get the effective dimension {na}: */
    ix_dim_t na = D->na; 
    /* Print the indexing parameters: */
    ix_print_parms(wr, pre, na, &(D->bp), D->sz, D->st, wd, suf);
  }
