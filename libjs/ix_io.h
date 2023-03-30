#ifndef indexing_io_H
#define indexing_io_H

/* Printout and debugging tools for {indexing_h} */
/* Last edited on 2023-03-18 11:21:09 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <sign.h>
#include <bool.h>

#include <ix.h>

void ix_print_dim ( FILE *wr, ix_dim_t d );
void ix_print_pos ( FILE *wr, ix_pos_t p );
  /* These procedures write the given parameter to {wr}, in decimal,
     bracketed by the strings {lp} and {rp} (which default to "" if NULL). */

void ix_print_indices ( FILE *wr, char *lp, ix_dim_t d, const ix_index_t ix[], int32_t wd, char *sp, char *rp );
void ix_print_sizes   ( FILE *wr, char *lp, ix_dim_t d, const ix_size_t sz[],  int32_t wd, char *sp, char *rp );
void ix_print_steps   ( FILE *wr, char *lp, ix_dim_t d, const ix_step_t st[],  int32_t wd, char *sp, char *rp );
  /* These procedures write the first {d} elements of the given tuple to {wr}, in decimal. 
    The tuple is bracketed by the strings {lp} and {rp} (which default to "" if NULL);
    In the case of {ix_print_steps}, the sign "+" or "-" is always printed.
    Each element is left-padded with ' ' to total width {wd} or more.
    The elements are separated by the string {sp} (which defaults to " " if NULL). */

void ix_print_parms
 ( FILE *wr, 
   char *pre, 
   ix_dim_t d, 
   ix_pos_t *bp, 
   const ix_size_t sz[], 
   const ix_step_t st[], 
   int32_t wd, 
   char *suf
 );
 /* Writes the indexing parameters {d}, {bp}, {sz[0..d-1]} and {st[0..d-1} to {wr},
   in readable ASCII format. The output consists of four lines 
    
      | "axes = {d}" 
      | "base = {bp}"
      | "step = {st[0]} {st[1]} ... {st[d-1]}"
      | "size = {sz[0]} {sz[1]} ... {sz[d-1]}"
    
    If any of the parameters {bp}, {sz}, and {st} is NULL, the corresponding line is omitted.
    Each line is prefixed by the string {pre} (which defaults to "" if NULL) 
    and suffixed with the string {suf} (which defaults to "\n" if NULL).  Each size or
    step is padded on the left to total width {wd} or more. Elements are separaed by 
    blanks. Positive steps will have an explicit "+" sign. */

#endif
