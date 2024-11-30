/* gauss_elim.h - basic tools for Gaussian triangulation and elimination. */
/* Last edited on 2024-11-26 20:17:08 by stolfi */

#ifndef gauss_elim_H
#define gauss_elim_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* In all the procedures below, two-dimensional matrices are stored
  into one-dimensional vectors, in row-by-row order. */

void gauss_elim_print_print_array
  ( FILE *wr,
    uint32_t indent,
    char *fmt,
    char *head,
    uint32_t m,
    uint32_t n,
    char *Mname,
    double M[],
    char *foot
  );
  /* Writes to {wr} the array {M}, in a human-readable format. Each
    element is printed with format {fmt}. The strings {head} and {foot},
    if not NULL, are printed on separate lines before and after the
    array. The {Mname}, if not {NULL}, is printed to the left of the
    middle row. Everything is indented by {indent} spaces. */

void gauss_elim_print_print_system
  ( FILE *wr, 
    uint32_t indent,
    char *fmt, 
    char *head, 
    uint32_t m, 
    uint32_t n, 
    char *Aname,
    double A[], 
    uint32_t p,
    char *Bname,
    double B[], 
    uint32_t q,
    char *Cname,
    double C[], 
    char *foot
  );
  /* Writes to {wr} the matrices {A} ({m×n}), {B} ({m×p}), and {C}
    ({m×q}) and the right-hand-side side by side, in a human-readable
    format. Each array element is printed with format {fmt}. The strings
    {head} and {foot}, if not NULL, are printed on separate lines before
    and after the system. The strings {Aname}, {Bname}, and {Cname} if
    not {NULL}, are assumed to be the names of the matrices,} and
    printed to the left of the respective middle rows.  
    
    If a matrix is {NULL}, it is omitted. If a matrix has zero columns,
    it is printed as brackets only, with no elements. Everything is
    indented by {indent} spaces. */


#endif
