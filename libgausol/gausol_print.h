/* gausol_print.h - printout functions for Gaussian full elimination. */
/* Last edited on 2024-11-28 16:40:48 by stolfi */

#ifndef gausol_print_H
#define gausol_print_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void gausol_print_array
  ( FILE *wr,
    uint32_t indent,
    char *fmt,
    char *head,
    uint32_t m, uint32_t prow[], uint32_t rh,
    uint32_t n, uint32_t pcol[], uint32_t rv,
    char *Mname, double M[], 
    char *foot
  );
  /* Writes to {wr} the array {M}, permuted by {prow,pcol} in a
    human-readable format. Each element is printed with format {fmt}.
    The strings {head} and {foot}, if not NULL, are printed on separate
    lines before and after the array. The {Mname}, if not {NULL}, is
    printed to the left of the middle row. Everything is indented by
    {indent} spaces.
    
    Specifically, assumes that {M} is an {m×n} matrix {M[0..m-1,0..n-1]}
    stored as a vector {M[0..m*n-1]}, with entry entry {M[i,j]} of the
    matrix stored in element {M[n*i+j]} of the vector.
    
    First prints {M} as-is. Then, if {pcol} and/or {prow} is not {NULL},
    prints also the permuted version {N} of {M} such that {N[i,j] =
    M[prow[i],pcol[j]}, for {i} in {0..m-1} and {j} in {0..n-1}. If not
    {NULL}, the table {prow[0..m-1]} must be a permutation of {0..m-1};
    if it is {NULL}, it is assumed to be the identity, that is,
    {prow[i]==i} for all those {i}. Likewise, either {pcol} is {NULL},
    or {pcol[0..n-1]} is a permutation of {0..n-1}.
    
    The parameters {rh} and {rv} specify a partition of the permuted array into
    up to four blocks. If {rh} is in {1..m-1}, prints a horizontal line
    separating the first {rh} rows of {N} from the other {m-rh} rows. If
    {rv} is in {1..n-1}, prints a vertical line separating the first {rv}
    columns of {N} form the other {n-rv} columns. */

void gausol_print_system
  ( FILE *wr, 
    uint32_t indent,
    char *fmt, 
    char *head, 
    uint32_t m, uint32_t prow[], uint32_t rh, 
    uint32_t n, uint32_t pcol[], uint32_t rv, char *Aname, double A[], 
    uint32_t p, char *Bname, double B[], 
    uint32_t q, char *Cname, double C[], 
  
    
    char *foot
  );
  /* Writes to {wr} the matrices {A} ({m×n}), {B} ({m×p}), and {C}
    ({m×q}) and the right-hand-side side by side, in a human-readable
    format. Each array element is printed with format {fmt}. The strings
    {head} and {foot}, if not NULL, are printed on separate lines before
    and after the system. The strings {Aname}, {Bname}, and {Cname} if
    not {NULL}, are assumed to be the names of the matrices,} and
    printed to the left of the respective middle rows.  
    
    If a matrix is {NULL}, or has zero columns, it is omitted. If the
    row count is {m} is zero, ony the {head} and {foot} are printed.
    Everything is indented by {indent} spaces.
    
    The vectors {prow} and {pcol} should be either {NULL} or
    Tpermutations of the indices {0..m-1} and {0..n-1}, respectively. If
    Tthey are noy both {NULL}, permuted version of the three arrays are
    Tprinted after the unpermuted ones. The row permutation {prow}
    Tapplies to all three matrices, and if {rh} is a number in {1..m-1},
    Tprints the horizontal line across all three after the first {rh}
    Trows. he column permutation {pcol} applies to {A} only, and if {rv}
    Tis in {1..n-1}, prints a vertical line between the the first {rv}
    Tcolumns and the last {n-rv}. */


void gausol_print_row_has_name_eq
  ( uint32_t i,
    uint32_t m,
    uint32_t rank,
    bool_t *at_dash_P,
    bool_t *at_data_P
  );
  /* Used when printing a matrix with {m} rows to properly center
    the " {name} = "  string at the left of matrix.
    
    Assumes that {rank} is in {0..m}, and the dash line is printed if and
    only if {rank} is in {1,m-1}, between the first {rank} rows and the last
    {m-rank} rows. The procedure also assumes that one is about to print
    row {i}, but will print first the row of dashes if {rank} is in
    {1..m-1} and {i==rank}. 
    
    The procedure sets {*at_dash_P} to true if and only if a dash line
    is to be printed and the " {name} = " should go at left of it.
    
    It sets {*at_data_P} true if and only if the " {name} = " should go
    at the left of the data line. */
    
#endif
