/* Infinite multi-dimensional arrays. */
/* Last edited on 2023-03-18 11:31:23 by stolfi */

/*  
  INFINITE ARRAYS
   
    A {ppv_grid_t} defines an /index reduction mapping/, a mapping from
    an infinite {d}-dimensional grid of integers to a finite one, such
    as the valid indices of a {ppv_array_t}. The mapping allows the
    finite grid to be shfted and extended or replicated along each
    dimension.
    
    This interface allows clients to read and write elements of sample
    arrays (such as described by {ppv_array_t}) as if they were
    infinite arrays, without worrying about index wrap-around.
    
  INTENDED USE

    A {ppv_grid_t} structure can be used to represent periodic signals,
    closed curves, repeating image patterns, functions defined on a
    {d}-dimensional toroid by sample interpolaion, etc.. In particular,
    it can be used to store the function samples in finite difference
    methods with mirror or inverting-mirror boundary conditions. */

#ifndef ppv_grid_H
#define ppv_grid_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <ppv_types.h>
#include <ppv_array.h>

/* !!! Give up on all these.  Keep only {NONE} (or {UNIF}, or {CLIP}) and {PREP}. !!! */

typedef enum 
  { ppv_grid_ext_NONE,  /* No extrapolation. */
    ppv_grid_ext_UNIF,  /* Extend with a single sample value. */
    ppv_grid_ext_CLIP,  /* Extend with border samples. */
    ppv_grid_ext_PREP,  /* Replicate fundamental cell. */
    ppv_grid_ext_NREP,  /* Replicate fundamental cell with negation. */
    ppv_grid_ext_PMIR,  /* Mirror the fundamental cell. */
    ppv_grid_ext_NMIR   /* Mirror the fundamental cell with negation. */
  } ppv_grid_ext_t;
  /* A {ppv_grid_ext_t} specifies a mapping from an arbitrary integer index {ix}
    to an index {ixr} with a finite non-empty range {0..n-1}. If {ix} is in {0..n-1}, 
    {ixr} is equal to {ix}.  Otherwise, let {r} be
    {imod(ix,r)} and {q=idiv(ix,r)}.  The options are:
    
      * {ppv_grid_ext_NONE}: do not extend. Attempt to reduce {ix}
      is an error.
      
      * {ppv_grid_ext_UNIF}: extend with unform value. The position of
      {AA[ix]} is {ix_pos_NONE}, fetching its value returns a
      /background/ value, and setting its value is an error.
      
      * {ppv_grid_ext_CLIP}: extend by clipping the index. The element
      {AA[ix]} is aiased to {a[0]} if {ix < 0}, and to {a[n-1]} if {ix
      >= n}.
      
      * {ppv_grid_ext_SREP}: extend by replication. The element
      {AA[ix]} is aliased to {a[r]}.
      
      * {ppv_grid_ext_NREP}: extend by replication and negation. The
      element {AA[ix]} is aliased to {a[r]} if {q} is even, or to the
      negation of {a[r]} if {q} is odd.
      
      * {ppv_grid_ext_SMIR}: extend by mirroring. The element {AA[ix]}
      is aliased to {a[r]} if {q} is even, or to {a[n-1-r]} if {q} is
      odd.
      
      * {ppv_grid_ext_NMIR}: extend by mirroring and negation. The
      element {AA[ix]} is aliased to {a[r]} if {q} is even, or to the
      negation of {a[n-1-r]} if {q} is odd.
      
    Since the elements of a {ppv_array_t} are uninterpreted bit strings, the
    aliasing to negated elements is indicated by a separate sign parameter. See
    {ppv_grid_get_sample} and related procedures. */

typedef 
  struct ppv_grid_t {
    ppv_dim_t d;                 /* Dimension (number of indices) 
    ppv_size_t *size;            /* The size of the finite grid along each axis. */
    ppv_grid_ext_t *ext;         /* Grid extrapolation method for each axis. */
    /* !!! Shifting only works for {SREP} grids. !!! */
    ppv_index_t *shift;          /* Shift to be applied to each axis. */
    /* !!! Skewing requires a triangular matrix of indices. !!! */
    /* !!! That is not compatible with index perms. !!! */
    /* !!! Should be a full {N×N} matrix, but how do we define that? !!! */
    ppv_index_t *skew;           /* Skewing increment for each axis. */
  } ppv_grid_t; 
  /* 
    A {ppv_grid_t} describes a mapping of the {d}-dimensional infinite grid 
    to the finite grid that has range {0..size[ax]-1} for each axis {ax}
    in {0..d-1}. */

/* GRID ALLOCATION */

ppv_grid_desc_t *ppv_grid_new
  ( ppv_dim_t d,
    ppv_size_t size[], 
    ppv_grid_ext_t ext[], 
    ppv_index_t shift[],
    ppv_index_t skew[]
  );
  /* Allocates a descriptor {AA} for a reduction mapping of the infinite {d}-dimensional 
    integer grid with the given parameters. The fini.
    
    The vectors {AA.ext,AA.shift,AA.skew} will be allocated as 
    part of the same record {*AA}.  Thefeore one must NOT call {free}
    on them.  Instead, {free(AA)} will reclaim all the storage used
    by the grid descriptor. */

/* SAMPLE INDEXING */

void ppv_grid_reduce ( ppv_grid_t *AA, ppv_index_t ix[], ppv_index_t ixr[] );
  /* Reduces the index vector [ix[0..d-1]} to finite ranges, as speciifed by {AA},
    and stores the result in {ixr[0..d-1]}; where {d} is {AA.d}. */

???

void ppv_grid_crop ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t skip, ppv_size_t  keep);
  /* Crops the array {A} to the range {skip..skip+keep-1} along axis
    {i}. Thus {size[i]} becomes {keep}, and using value {r} for that
    index in the new descriptor is equivalent to using {r + skip} in
    the old one.
    
    The value of {skip+keep} must not exceed the original {size[i]}.
    {keep} may be zero, in which case the array becomes empty (but
    retains its size along the other axes) and all its steps are
    reset to zero. */

void ppv_grid_subsample ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {i}, so that using {r}
    for that index in the new descriptor is equivalent to using
    {r*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[i]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[i]}. In particular, {size[i]} will be 0 
    if and only if it was 0 originally. */

void ppv_grid_flip ( ppv_grid_t *A, ppv_axis_t i );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {i}
    runs in the opposite direction.  */

void ppv_grid_replicate ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {i}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {step[i]} to 0,
    {size[i]=sz}. On input, {size[i]} must be 1, and {sz} must be
    positive. */

void ppv_grid_swap_indices ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that indices {i} and
    {j} are exchanged. Thus, for instance, {ppv_swap_indices(img, 1,2)} 
    could be used to exchange the rows and columns of a 2D color image,
    assuming index 0 is the color channel.
    If {i == j}, the procedure has no effect. */

void ppv_grid_flip_indices ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the indices {i..j} in reverse order. Thus,
    after the call the call {ppv_flip_indices(A,2,5)}, 
    element {A[p,q,r,s,t,u]} is equivalent to {A[p,q,u,t,s,r]}
    of the original array. */

void ppv_grid_diagonal ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes a diagonal slice
    through the original array, cut parallel to the bisector of axes
    {i} and {j}; and redefines axis {i} to be along that diagonal.
    Thus, using index values {r} and {s} along those axes in the 
    new array is equivalent to using {r} and {s+r} in the original 
    array.
    
    The procedure requires {i != j}, and {0 < A.size[i] <=
    A.size[j]}. It reduces {A.size[j]} by {A.size[i]-1}, and adds
    {A.step[j]} to {A.step[i]} .
    
    Thus, for example, if {A} is a color image with shape
    {3,640,480,1,1,1}, the call {ppv_diagonal(A,2,1)} will cut it down
    to a diagonal band starting at the main diagonal and extending
    {161=640-(480-1)} pixels to the left, and shear it horizontally so
    that it becomes a rectangular image with shape {3,161,480,1,1,1}.
    Then sample {A[c,h,v,0,0,0]} of the new image will be sample
    {A[c,h+v,v,0,0,0]} of the original one. */

void ppv_grid_chop ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t sz, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that index {i} is
    replaced by the combination of indices {i} and {j} (which must
    be different). Upon entry, {size[j]} must be 1. Upon exit,
    {size[i]} will be {sz} and {size[j]} will be {oldsz/sz} where {oldsz} is
    the original value of {size[i]}. 
    
    The original size {oldsz} must be a multiple of {sz}, then the remainder elements
    will become inacessible.  Inparticular, if {oldsz} is less than 
    {sz}, the array will become empty.
    
    For exmple, {ppv_chop(V, 0, 100, 1)} will turn a vector of 5000
    samples (along axis 0) into a greyscale image with 100 samples per
    row (along index 0) and 50 rows (along axis 1). As another
    example, {ppv_chop(A, 1, 64, 3)} could be used to chop a color
    image with shape {3,640,480,1,1,1} into 10 vertical stripes, 64
    pixels wide; and stack those stripes in the depth direction,
    yielding a `color tomogram' with shape {3,64,480,10,1,1}.
    Thus sample {A[c,h,v,d,0,0]} of the new array will be 
    the same as {A[c,64*d+h,v,0,0,0]} of the original. */

/* DESCRIPTOR PRINTOUT */

void ppv_grid_print_descriptor ( FILE *wr, char *pf, ppv_grid_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* DESCRIPTOR VALIDATION */

bool_t ppv_grid_descriptor_is_valid ( ppv_grid_t *A, bool_t die );
  /* Checks whether the descriptor {A} seems to be valid.
    
    To be valid, an array descriptor must satisfy the following
    conditions:

      (1) Two samples with distinct valid index tuples have the same
        position only if {step[i]} is zero for every axis where
        their indices differ.

      (2) For any {i}, if {size[i]} is 1, then {step[i]} is zero.

      (3) If the array is memoryless, then {base} and all {step}
        coefficients are zero, and {el} is NULL.

      (4) The field {npos} is the number of samples with distinct
        positions. Because of (3), {npos} is zero if any {size[i]}
        is zero, otherwise it is the product of {size[i]}
        for every {i} where {step[i]} is nonzero.

      (5) The storage area of every valid sample is contained
        in the allocated area pointed by {el}. (Note that 
        when {bps==0} this is true even if {el==NULL}.)
        
    The procedure returns TRUE if the parameters are valid.
    Otherwise, it either aborts the program with an error message 
    (if {die=TRUE}) or returns FALSE (if {die=FALSE}).
  */

/* DETAILS OF SAMPLE PACKING

  General principles

    The samples of an array {A} are stored in a vector of words,
    pointed to by {A.el}. Each word is an unsigned binary integer with
    {A.bpw} bits.

    A sample is never split across word boundaries; so, if {bps} is
    not an exact divisor or multiple of {bpw}, any incompletely filled
    words will be padded with zeros at the `big' end.
    
    For this implementation, the word size {A.bpw} must be 8, 16, or
    32. The implementation below assumes that words can be directly
    addressed by a pointer of the appropriate type, and that an
    {uint32_t} variable is at least 32 bits long.

  Packing of "small" samples

    If {bps} is less than {bpw}, then {K=bpw/bps} samples are stored
    into each word. Within each word, samples with consecutive
    positions are packed together in big-endian order. Thus, a sample
    whose position is {pos} will be stored in word {iw = pos/K}, and
    its units bit will be bit {ib = (K-1 - pos%K)*bps} of that word.

    For example, consider an image with {bps=6} and 5 samples
    {A,B,C,D,E} with positions {pos = 0,1,2,3,4}. With {bpw=8}, we
    would have {K = 8/6 = 1} samples per word, and the samples would
    be stored in memory as follows:

      | word#   <--0---> <--1---> <--2---> <--3---> <--4--->
      | sample#   <-0-->   <-1-->   <-2-->   <-3-->   <-4-->      
      |         ==AAAAAA ==BBBBBB ==CCCCCC ==DDDDDD ==EEEEEE
      | bit#    7      0 7      0 7      0 7      0 7      0   

    With {bpw = 16}, we would have {K = 16/6 = 2}, and the following layout:

      | word#   <------0-------> <-----1--------> <-----2-------->
      | sample#     <-0--><-1-->     <-2--><-3-->     <-4-->      
      |         ====AAAAAABBBBBB ====CCCCCCDDDDDD ====EEEEEE======
      | bit#    1  1     0     0 1  1     0     0 1  1     0     0   
      |         5  2     6     0 5  2     6     0 5  2     6     0   

    With {bpw = 32}, we would have {K = 32/6 = 5}, and the following layout:

      | word#   <---------------------1-------->
      | sample#   <-0--><-1--><-2--><-3--><-4-->      
      |         ==AAAAAABBBBBBCCCCCCDDDDDDEEEEEE
      | bit#    33     2     1     1     0     0
      |         10     4     8     2     6     0

  Packing of "big" samples

    If {bps} is greater than {bpw}, then each sample will occupy 
    {K=(bps+bpw-1)/bpw} consecutive words, big-end first. Thus, if 
    {bps = 18}, the layout for {bpw = 8} will be

      | word#   <--0---> <--1---> <--2---> <--3---> <--4---> <--5---> ...
      | sample#       <-------0---------->       <--------1---------> ...
      |         ======AA AAAAAAAA AAAAAAAA ======BB BBBBBBBB BBBBBBBB ...
      | bit#    7    2 0 7      0 7      0 7    2 0 7      0 7      0 ... 

    while that for {bpw = 16} will be

      | word#   <------0-------> <-----1--------> <-----2--------> <-----3--------> ...
      | sample#               <---------0------->               <--------1--------> ...
      |         ==============AA AAAAAAAAAAAAAAAA ==============BB BBBBBBBBBBBBBBBB ...
      | bit#    1            0 0 1              0 1            0 0 1              0 ... 
      |         5            2 0 5              0 5            2 0 5              0 ...  

  Sample addressing

    The sample storage area address {A.el} is the address of the word
    that contains the highest-order bit of the sample with position 0.

*/

#endif
