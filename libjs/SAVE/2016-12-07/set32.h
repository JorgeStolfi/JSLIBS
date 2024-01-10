#ifndef set32_H
#define set32_H

/* Subsets of the set {0..31}. */
/* Last edited on 2013-10-25 01:12:05 by stolfilocal */

#include <stdint.h>

typedef uint32_t set32_t;
  /* A subset of {0..31}. */

typedef int8_t set32_elem_t;
  /* A set element, in {0..31}, or {set32_NO_ELEM}. */

#define set32_ELEM_MAX (31)
  /* Max valid value of a {set32_elem_t}. */

#define set32_NO_ELEM (-1)
  /* An invalid {set32_t} value. */  

#define set32_belongs(el,A)      (((A) >> (el)) & 1)
#define set32_include(el,A)      ((A) | (1 << (el)))
#define set32_exclude(el,A)      ((A) & (~(1 << (el))))

#define set32_union(A,B)         ((A) | (B))
#define set32_intersection(A,B)  ((A) & (B))
#define set32_difference(A,B)    ((A) & ~(B))

#define set32_range(i,j)         ((1 << ((j)+1)) - (1 << (i)))
  /* The range {i..j}, with {i,j} in {0..31}; empty iff {i > j}. */

int set32_count(set32_t A);
  /* The cardinality of the set {A}.  */

typedef int8_t set32_index_t;
  /* The index (rank) of an element in a set, in order of
    increasing element value. Valid values are {0..31}. */

#define set32_NO_INDEX ((set32_index_t)(-1))
  /* An invalid {set32_index_t} value. */  

set32_elem_t set32_elem_from_index(set32_index_t j, set32_t A);
  /* Returns the {j}th smallest element in the set {A}, counting 
    from {j=0} to {set32_count(A)-1}. Returns {set32_NO_ELEM} if
    {j > set32_count(A)}. */

set32_index_t set32_index_from_elem(set32_elem_t el, set32_t A);
  /* Returns the index (rank) of element {el} in the set of axes
    {A}, sorted by increasing order. Returns {set32_NO_INDEX} if
    {el} is not in the set {A}. */

#endif
