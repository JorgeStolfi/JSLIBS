#ifndef indexing_box_H
#define indexing_box_H

/* Index boxes -- domains of multi-dimensional arrays. */
/* Last edited on 2024-11-15 19:13:59 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <sign.h>
#include <bool.h>

#include <ix.h>

/* For the procedures in this module, a box of {\RZ^d} is described by
  its low index vector {ixlo[0..d-1]} and it size vector {size[0..d-1}
  along each axis.

  A {NULL} {ixlo} or {size} parameter is interpreted as a vector of {d}
  zeros.

  If an input box has {size[k]==0} for any axis {k}, the box is assumed
  to be empty. If the result is empty, {ixlo} and {size} are both set to
  all zeros. */

bool_t ixbox_is_empty(ix_dim_t d, ix_size_t size[]);
  /* Returns {TRUE} if {size} is {NULL}, or any one of the components
    {size[0..d-1]} is zero. */

bool_t ixbox_has(ix_dim_t d, ix_index_t ix[], ix_index_t ixlo[], ix_size_t size[]);
  /* Returns {TRUE} if and only if {ix[k]} is in the box {B} defined by
    {ixlo} and {size}; specifically, if {ix[k] - ixlo[k]} is in the
    range {0 .. size[k]-1} for every {k} in {0..d-1}. In particular,
    returns {FALSE} if the box {B} is empty. */

void ixbox_intersect
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[],
    ix_index_t ixlo_R[], 
    ix_size_t  size_R[]
  );
  /* Computes their intersection {R} of two {d}-dimensional boxes --
    {A}, defined by the low-corner {ixlo_A[0..d-1]} and the sizes
    {size_A[0..d-1]}, and {B}, defined similarly by {ixlo_B} and
    {size_B}. Stores the result into {ixlo_R} and {size_R}.
    
    The parameters {ixlo_R} and {size_R} may not be {NULL}. */
    
bool_t ixbox_is_contained
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[]
  );
  /* Returns {TRUE} iff the box {A} defined by {ixlo_A,size_A} is
    contained or equal in the box {B} defined by {ixlo_B,size_B}. In
    particular, returns {TRUE} if {A} empty, and {FALSE} if {B} is empty
    but {A} is not. */

bool_t ixbox_equal
  ( ix_dim_t d,
    ix_index_t ixlo_A[], 
    ix_size_t  size_A[],
    ix_index_t ixlo_B[], 
    ix_size_t  size_B[]
  );
  /* Returns {TRUE} iff the box {A} defined by {ixlo_A,size_A} is
    equal to the box {B} defined by {ixlo_B,size_B}. In
    particular, returns {TRUE} if {A} and {B} are both empty. */


void ixbox_print
  ( FILE *wr, 
    ix_dim_t d, 
    char *pref, 
    ix_index_t ixlo[], 
    ix_size_t size[], 
    char *suff
  );
  /* Writes the elements {ixlo[0..d-1]} and {size[0..d-1]} as "0(10) 10(100) ..." 
    preceded by {pref} and followed by {suff}, if these strings are not {NULL}. */

#endif
