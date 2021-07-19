#ifndef array_H
#define array_H
/* Self-bounded multidimensional arrays of things */

#define spmat_H_COPYRIGHT "Copyright © 2009 by J. Stolfi, UNICAMP"
/* Created on 2009-08-31 by J.Stolfi, UNICAMP */
/* Last edited on 2019-08-30 07:24:56 by jstolfi */

/* SELF-BOUNDED ARRAYS

  This interface defines self-bounded arrays with elements of any
  specified type. Unlike standard C arrays, the number of elements and
  the addressing formula is part of the array's representation.
  
  The representation consists of a record ("dope vector", "descriptor")
  that contains the parameters of the addressing formula, and a pointer
  to a separate memory area where the elements reside.Thus, when passing
  a self-bounded array to a procedure, there is no need to pass
  the array dimensions and indexing coefficientes as separate parameters.
  
  The addressing formula is a completely general affine ("linear") 
  polynomial.  Therefore,  many slicing and reindexing operations can be
  performed by manipulation of the descriptor, without copying the elements.
  
  This interface also defines operations for copying elements between
  arrays and ordinary C vectors. */

#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <ix.h>
#include <ix_descr.h>

/* DECLARING A NEW SPARSE ARRAY TYPE


  Here is how one would define the type {double_array_t} as a multidimensional
  array of doubles:
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <array.h>
    
    array_typedef(double_array_t, double_array, double);
    ------------------------------------------------------------
    
    ------------------------------------------------------------
    / * Implementations  * /

    array_impl(double_array_t, double_array, double);
    ------------------------------------------------------------

  As another example, here is how one would define a type {graph_t}
  as being a directed graph, represented by a boolean adjacency
  matrix:
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <array.h>
    
    array_typedef(graph_t, graph, bool_t);
    ------------------------------------------------------------

    ------------------------------------------------------------
    / * Implementations  * /
    
    array_impl(graph_t, graph, bool_t);
    ------------------------------------------------------------
*/

#define array_typedef(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_TYPEDEF(ARRAY_TYPE,PREFIX,ELEM_TYPE)
/* 
  This macro will also declare prototypes for the functions
  
    ------------------------------------------------------------
    PREFIX##_new
    PREFIX##_make_desc
    PREFIX##_copy
    ------------------------------------------------------------
      
  which are described below.
  
  The names {ARRAY_TYPE} and {PREFIX} can be chosen quite arbitrarily.
  However, if {ELEM_TYPE} is some previously defined type named
  {{XXX}_t} or {XXX}, it is recommended to use {XXX}_array_t}
  as the {ARRAY_TYPE}, and {{XXX}_array} as the {PREFIX}. */

#define array_typeimpl(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_TYPEIMPL(ARRAY_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the implementations of the functions whose
  prototypes were declared with {array_typedef}. It should be
  called after the corresponding call to {array_typedef}, which
  declares the prototypes for those procedures. */
 
/* INDEXING TYPES */
 
#define array_MAX_AXES (ix_descr_MAX_AXES)
  /* Max number of indices of an {array_t}. */

/* ------------------------------------------------------------
  typedef ix_dim_t PREFIX##_dim_t;
------------------------------------------------------------ */
   /* Order (number of indices) of an array or sub-array.  
     Legitimate values are {0..array_MAX_AXES}. */
   
/* ------------------------------------------------------------
  typedef ix_axis_t PREFIX##_axis_t;
------------------------------------------------------------ */
 /* Specifies an axis (index). Legitimate values are {0..array_MAX_AXES-1},
     but {ix_MAX_AXES} may be used as a `null' value. */

/* ------------------------------------------------------------
  typedef ix_size_t PREFIX##_size_t;
------------------------------------------------------------ */
  /* Number of elements along a particular axis. */

/* ------------------------------------------------------------
  typedef ix_index_t PREFIX##_index_t;
------------------------------------------------------------ */
  /* Index of an element. Note that it is signed, even though 
    legitimate values are always non-negative. */

/* ------------------------------------------------------------
  typedef ix_step_t PREFIX##_step_t;
------------------------------------------------------------ */
   /* Element position increment along some axis. */

/* ------------------------------------------------------------
  typedef ix_count_t PREFIX##_count_t;
------------------------------------------------------------ */
  /* Count of elements, for example in whole array. */

/* ------------------------------------------------------------
  typedef ix_pos_t PREFIX##_pos_t;
------------------------------------------------------------ */
  /* Position (index) of an element in the linearized storage vector. */

/* ARRAY ALLOCATION */

/* ------------------------------------------------------------
ARRAY_TYPE PREFIX##_new ( ix_dim_t na, ix_size_t sz[] );
------------------------------------------------------------ */
  /* Allocates a new multidimensional array with size {sz[i]}
    along axis {i}, for {i} in {0..na-1}. 
    
    The elements are stored in consecutive positions, in 
    lexicographic index order (with
    the /last/ index varies fastest, as in C and Pascal). */

/* ------------------------------------------------------------
ARRAY_TYPE PREFIX##_copy( ARRAY_TYPE *A );
------------------------------------------------------------ */
  /* Allocates a new multidimensional array {C} with the same number
    of axes and same size as {A} along each axis, by
    {PREFIX##_new(A.ds.na, A.ds.sz, ixor)}. Then copies all the elements
    of {A} into it with {PREFIX##_assign(C, A)}. The array {A} should
    not have any non-trivial index with zero increment. */

/* DESCRIPTOR ALLOCATION */

/* ------------------------------------------------------------
ARRAY_TYPE *PREFIX##_new_descr ( void );
------------------------------------------------------------ */
  /* Allocates a new array descriptor record, with uninitialized indexing
    formula, and {el = NULL}. */

/* ------------------------------------------------------------
ARRAY_TYPE *PREFIX##_copy_descr ( ARRAY_TYPE *A );
------------------------------------------------------------ */
  /* Allocates a new array descriptor record, and copies into it the
    contents of descriptor {A}, including the {el} pointer.
    Thus, the resulting array with share its elements with {A}. */

/* DEALLOCATION */

/* ------------------------------------------------------------
void PREFIX##_free_elems ( ARRAY_TYPE *A );
------------------------------------------------------------ */
  /* Reclaims the element storage area {A.el} and sets {A.el = NULL}.  
    Does *not* reclaim the descriptor record {A}. */

/* ------------------------------------------------------------
void PREFIX##_free_descr ( ARRAY_TYPE *A );
------------------------------------------------------------ */
  /* Discards the descriptor of array {A} (but not the element storage
    area {A.el}). */

/* ACESSING INDIVIDUAL ELEMENTS */

/* ------------------------------------------------------------
ELEM_TYPE PREFIX##_get_elem ( ARRAY_TYPE *A, ix_index_t ix[] );
------------------------------------------------------------ */
  /* Returns the element of {A} with indices {ix[0..A.ds.na-1]}. */

/* ------------------------------------------------------------
PREFIX##_pos_t *PREFIX##_get_elem_position ( ARRAY_TYPE *A, ix_index_t ix[] );
------------------------------------------------------------ */
  /* Returns the position (linearized index) in {A.el} of the element 
    with indices {ix[0..A.ds.na-1]} of array {A}. */

/* ------------------------------------------------------------
void PREFIX##_set_elem ( ARRAY_TYPE *A, ix_index_t ix[], ELEM_TYPE v );
------------------------------------------------------------ */
  /* Stores {v} into the element of array {A} with indices {ix[0..A.ds.na-1]}. */

/* ASSIGNMENT */
  
/* ------------------------------------------------------------
void PREFIX##_assign ( ARRAY_TYPE *A, ARRAY_TYPE *B );
------------------------------------------------------------ */
  /* Sets every element of {A} to the corresponding element of {B}. The arrays
    must have the same set of valid index tuples.
    
    More precisely, it enumerates all index tuples of {A}, in lexicographic order,
    and assigns the corresponding elements of {A} and {B}.
    Note that elements of {A} that have more than one index tuple are assigned
    multiple times. Also, complicated effects may occur if any element is shared
    by {A} and {B} with different indices. */

/* GET/CHECK SIZE */
  
/* ------------------------------------------------------------
void PREFIX##_get_size ( ARRAY_TYPE *A, ix_size_t sz[] );
------------------------------------------------------------ */
  /* Stores the size of {A} along axis {i} into {sz[i]},
    for {i} in {0..A.ds.na-1}. */
  
/* ------------------------------------------------------------
void PREFIX##_check_size ( ARRAY_TYPE *A, ix_dim_t na, ix_size_t sz[] );
------------------------------------------------------------ */
  /* Fails if the {A} does not have {na} axes, or if its
    size along some axis {i} in {0..A.ds.na-1} is not {sz[i]}. */

/* LIMITS */
 
#define array_MAX_ABS_STEP ix_MAX_ABS_STEP
  /* Max absolute value of a {step} field. */

#define array_MAX_ELEMS ix_MAX_POSITIONS
  /* Max number of *distinct* element positions, not the
    the apparent domain size (the number of valid index tuples). The
    latter can exceed this limit if the array has replicated elements
    (null {step}s). */

#define array_MAX_SIZE (array_MAX_ELEMS)
  /* In theory a memoryless array could have a {size} greater than
    {array_MAX_ELEMS}, but that freedom seems to be of little use, and
    makes the code more complicated. Hence we enforce this limit even
    for memoryless arrays. */

#define array_MAX_INDEX (array_MAX_ELEMS - 1)
  /* We assume that {array_MAX_INDEX + 1} does not overflow an {ix_index_t}.
    The type is signed in order to avoid underflow, e.g. in loops like
    {for (index = N; index >=0; index--).} */

#define array_MAX_POS (array_MAX_ELEMS - 1)
  /* Valid positions run from 0 to the total number of elements minus 1. */

#define array_MAX_BYTES SIZE_MAX
  /* The size argument for {malloc} is a {size_t}. */

#include <array_def.h> 
/*    
    ??? TO DO: Get rid of {ixor}, use {ix_descr_flip_indices} to get F-order; ??? 
    !!! TO DO: in {array_copy}, allow and honor virtually replicated indices. !!!
    
*/
#endif
