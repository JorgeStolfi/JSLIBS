/* Infinite multi-dimensional arrays. */
/* Last edited on 2021-06-10 13:36:39 by jstolfi */

/*  
  INFINITE ARRAYS
   
    A {ppv_grid_t} implements an /infinite array/ --- a
    multi-dimensional array whose size along some dimensions appears
    to be infinite. It is actually implemented by a finite array whose
    contents is infinitely replicated along those dimensions.
    
    This interface allows clients to read and write elements
    of infinite arrays without worrying about index wrap-around. 
    It also provides various descriptor operations, such as
    subsampling, shifting, flipping and index permutations. 
    
  INTENDED USE

    A {ppv_grid_t} structure can be used to represent periodic 
    signals, closed curves, repeating image patterns, 
    functions defined on the torus, etc..  It can be used
    to store the function samples in finite difference methods
    with mirror or inverting-mirror boundary conditions. */

#ifndef ppv_grid_H
#define ppv_grid_H

#include <ppv_types.h>
#include <ppv_array.h>
#include <stdio.h>

#define ppv_grid_NAXES ppv_array_NAXES
  /* Number of indices in a {ppv_grid_t}. */

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
  /* A {ppv_grid_ext_t} specifies how to extend an array {a[ix]} with a
    finite index range {0..n-1} to an array {A[ix]} defined for all
    integers. If {ix} is in {0..n-1}, the value and position of {A[ix]} are always 
    those of {a[ix]}.  For an index {ix} outside {0..n-1}, let {r} be
    {imod(ix,r)} and {q=idiv(ix,r)}.  The options are:
    
      * {ppv_grid_ext_NONE}: do not extend. The position of {A[ix]} is
      {ix_pos_NONE}, and fetching or setting its value is an error.
      
      * {ppv_grid_ext_UNIF}: extend with unform value. The position of
      {A[ix]} is {ix_pos_NONE}, fetching its value returns a
      /background/ value, and setting its value is an error.
      
      * {ppv_grid_ext_CLIP}: extend by clipping the index. The element
      {A[ix]} is aiased to {a[0]} if {ix < 0}, and to {a[n-1]} if {ix
      >= n}.
      
      * {ppv_grid_ext_SREP}: extend by replication. The element
      {A[ix]} is aliased to {a[r]}.
      
      * {ppv_grid_ext_NREP}: extend by replication and negation. The
      element {A[ix]} is aliased to {a[r]} if {q} is even, or to the
      negation of {a[r]} if {q} is odd.
      
      * {ppv_grid_ext_SMIR}: extend by mirroring. The element {A[ix]}
      is aliased to {a[r]} if {q} is even, or to {a[n-1-r]} if {q} is
      odd.
      
      * {ppv_grid_ext_NMIR}: extend by mirroring and negation. The
      element {A[ix]} is aliased to {a[r]} if {q} is even, or to the
      negation of {a[n-1-r]} if {q} is odd.
      
    Since the elements of a {ppv_array_t} are uninterpreted bit strings, the
    aliasing to negated elements is indicated by a separate sign parameter. See
    {ppv_grid_get_sample} and related procedures. */

typedef 
  struct ppv_grid_t {
    ppv_array_t a;         /* The fundamental array. */
    ppv_grid_ext_t ext[ppv_grid_NAXES];    /* Grid extrapolation methods. */
    /* !!! Shifting only works for {SREP} grids. !!! */
    ppv_index_t    shift[ppv_grid_NAXES];  /* Index shifts. */
    /* !!! Skewing requires a triangular matrix of indices. !!! */
    /* !!! That is not compatible with index perms. !!! */
    /* !!! Should be a full {N×N} matrix, but how do we define that? !!! */
    ppv_index_t skew[ppv_grid_NSKEW];  /* Skewing increment for each index. */
  } ppv_grid_t; 
  /* 
    A {ppv_grid_t} is descriptor for a multi-dimensional infinite of
    small non-negative integers (/samples/), stored in memory 
    in a packed binary format. 
    
    Grid indexing
    
      Each sample of a {ppv_grid_t} {A} is identified by {N=ppv_NAXES}
      integer indices {[ix[0],..ix[N-1]]}.

    Sample indexing

      Each sample of a grid {A} has a /position/ which is an
      affine combination of its {N} indices {ix[0],..ix[N-1]}:

        | pos(A, ix) = A.base + ix[0]*A.step[0] + ·· ix[N-1]*A.step[N-1]

      Each coefficient {A.step[i]} may be positive, negative, or zero.
      In the last case, two samples whose indices differ only in that
      index will have the same position. The field {A.npos} records the
      number of distinct sample positions in the array.

    Sample storage

      The samples of an array {A} are stored in the area pointed by
      {A.el}. The position {pos} of a sample, together with the the
      packing parameters {A.bps} and {A.bpw}, uniquely determines the
      location of the sample within that area. Two samples with distinct
      positions occupy disjoint areas of memory.

      In order to save storage and to simplify image transfer between
      devices and computers with different architectures, the storage
      area is organized as a vector of /words/. The number of bits per
      word is an explicit attribute {A.bpw} of the array, and may be
      different from the machine's native word size. The legal values
      for {bpw}, for this implementation, are 8, 16, and 32. Depending
      on {A.bps} and {A.bpw}, one sample may span two or more words, or
      one word may contain two or more samples. See the section
      "DETAILS OF SAMPLE PACKING" at end of this file.  */

/* ARRAY ALLOCATION */

ppv_grid_t ppv_new_grid ( ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw );
  /* Returns a descriptor {A} for a newly allocated array of samples,
    with {A.size = {sz[0],.. sz[N-1]}} and the specified number of
    bits per sample and per word. The individual sizes {sz[i]} must
    not exceed {ppv_MAX_SIZE},
    
    If the array {A} turns out to be memoryless, {A.el} is set to NULL
    and all samples will have the same position (0).
    
    Otherwise, {A.el} will point to a newly allocated area large
    enough to contain all samples. The samples {A[ix[0],..ix[N-1]]} will
    be packed as compactly as possible, linearized in the obvious way
    with index {ix[0]} varying fastest and {ix[N-1]} the slowest. That
    is, {A.base} will be 0, and each increment {A.step[i]} will be
    the product of all sizes {A.size[j]} with {j<i}. The product
    {sz[0]*..sz[N-1]} must not exceed {ppv_MAX_SAMPLES}, and the total
    storage occupied by them must not exceed {ppv_MAX_BYTES}.
    
    In any case, the call {free(A.el)} will reclaim the storage used
    by the array. */

/* SAMPLE INDEXING */

ppv_pos_t ppv_sample_pos ( ppv_grid_t *A, ppv_index_t *ix );
  /* Computes the position of sample {A[ix[0],.. ix[N-1]]}. Does not check
    whether the element exists or not. */

/* SAMPLE EXTRACTION AND INSERTION AT POSITION */

/* For the procedures below, {pos} is the position of a sample in a
  sample storage area with address {el}. The samples are assumed to be
  packed with {bps} bits per sample and {bpw} bits per word; see the
  section "DETAILS OF SAMPLE PACKING" below. The user must make sure
  that {pos} is valid (i.e. the sample actually exists.)  */

ppv_sample_t ppv_get_sample_at_pos 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos 
  );
  /* Extracts the sample with position {pos}. */

void ppv_set_sample_at_pos 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos, 
    ppv_sample_t qv 
  );
  /* Stores value {qv} into the sample with position {pos}. */

/* INDEX TUPLE MANIPULATION */

void ppv_index_clear ( ppv_index_t *ix );
  /* Sets {ix[i] = 0} for every axis {i}. */

void ppv_index_assign ( ppv_index_t *ix, ppv_index_t *val );
  /* Sets {ix[i] = val[i]} for every axis {i}. */

bool_t ppv_index_first ( ppv_index_t *ix, ppv_grid_t *A );
  /* Sets {ix[i] = 0} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

bool_t ppv_index_last ( ppv_index_t *ix, ppv_grid_t *A );
  /* Sets {ix[i] = A.size[i]-1} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

void ppv_index_shift ( ppv_index_t *ix, ppv_index_t *inc );
  /* Sets {ix[i] += inc[i]} for every axis {i}. */

sign_t ppv_index_compare ( ppv_index_t *ixa, ppv_index_t *ixb );
  /* Returns {NEG}, {ZER}, or {POS} depending on whether {ixa} is less
    than, equal to, or greater than {ixb} in reverse lexicographic
    order (where the first index varies the fastest). */

bool_t ppv_index_next ( ppv_index_t *ix, ppv_grid_t *A, ppv_dim_t d, ppv_pos_t *p );
  /* Set {ix[0..d-1]} to the next combination of values of those
    indices that is valid for {A}, in reverse lexicographic order
    (with the first index in the innermost loop, as in FORTRAN array
    storage).
  
    More precisely, tries to increment {ix[0]} first. If an index
    {ix[i]} to be incremented has reached its limit {A.size[i]-1}, resets
    all indices {ix[0..i]} to 0, and tries to increment index
    {ix[i+1]} instead. If this `carry over' falls off the index tuple
    -- that is, if the tuple {ix} changes from {A.size[0..d-1] - (1,..1)} back to
    {(0,..0)} -- the procedure returns TRUE; otherwise it returns
    FALSE. The effect is undefined if, for some {i}, {A.size[i]} is 0
    or {ix[i]} is not in {0..A.size[i]-1}.
    
    If {p} is not null, assumes that {*p} is the position corresponding
    to the index tuple {ix}, and adjusts it to account for the change,
    using the steps {st}. (If {p} is NULL, {st} is not used.)
    
    This procedure can be used to vary the index tuple {ix} over all
    the elements of a non-empty array {A}, in a single loop, e.g.
    
      | if (ppv_index_first(ix,A))
      |   { ppv_pos_t p = ppv_position(ix,A);
      |     do { ... } while (! ppv_index_next(ix,A,&p));
      |   }

    */

bool_t ppv_index_prev ( ppv_index_t *ix, ppv_grid_t *A, ppv_dim_t d, ppv_pos_t *p );
  /* Like {ppv_index_next}, but in the reverse order. Each index {ix[i]}
    with {i} in {0..d-1} is decremented from {A.size[i]-1} down to 0,
    and reset to {(A.size[i]-1)} when needed; and the returned result
    is TRUE when {ix} changes from {(0,..0)} back to {A.size[0..d-1] -
    (1,..1)}.
    
    Thus, to scan all the elements of an array {A} in reverse order, use 

      | if (ppv_index_last(ix,A))
      |   { ppv_pos_t p = ppv_position(ix,A);
      |     do { ... } while (! ppv_index_prev(ix,A,&p));
      |   }

   */
    
/* INDEX TUPLE VALIDATION */
    
bool_t ppv_index_is_valid ( ppv_index_t *ix, ppv_grid_t *A );
  /* Returns TRUE iff {ix} is a valid index tuple for the array {A},
    i.e., if sample {A[ix[0],..ix[N-1]]} exists. This is true if and
    only if {ix[i]} lies in the range {0..A.step[i]-1}, for every
    axis {i}. */

/* DESCRIPTOR MANIPULATION */

/* The following procedures modify the {size}, {step} and {base} field
  of an array descriptor {A}, so as to change the set of valid indices
  and/or the mapping between indices and sample positions. No samples are
  actually allocated, reclaimed, or changed. 
  
  All these procedures preserve the validity of the descriptor (in the
  sense of {ppv_descriptor_is_valid} below) when given valid
  arguments. */

void ppv_crop ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t skip, ppv_size_t  keep);
  /* Crops the array {A} to the range {skip..skip+keep-1} along axis
    {i}. Thus {size[i]} becomes {keep}, and using value {r} for that
    index in the new descriptor is equivalent to using {r + skip} in
    the old one.
    
    The value of {skip+keep} must not exceed the original {size[i]}.
    {keep} may be zero, in which case the array becomes empty (but
    retains its size along the other axes) and all its steps are
    reset to zero. */

void ppv_subsample ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {i}, so that using {r}
    for that index in the new descriptor is equivalent to using
    {r*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[i]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[i]}. In particular, {size[i]} will be 0 
    if and only if it was 0 originally. */

void ppv_flip ( ppv_grid_t *A, ppv_axis_t i );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {i}
    runs in the opposite direction.  */

void ppv_replicate ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {i}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {step[i]} to 0,
    {size[i]=sz}. On input, {size[i]} must be 1, and {sz} must be
    positive. */

void ppv_swap_indices ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that indices {i} and
    {j} are exchanged. Thus, for instance, {ppv_swap_indices(img, 1,2)} 
    could be used to exchange the rows and columns of a 2D color image,
    assuming index 0 is the color channel.
    If {i == j}, the procedure has no effect. */

void ppv_flip_indices ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the indices {i..j} in reverse order. Thus,
    after the call the call {ppv_flip_indices(A,2,5)}, 
    element {A[p,q,r,s,t,u]} is equivalent to {A[p,q,u,t,s,r]}
    of the original array. */

void ppv_diagonal ( ppv_grid_t *A, ppv_axis_t i, ppv_axis_t j );
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

void ppv_chop ( ppv_grid_t *A, ppv_axis_t i, ppv_size_t sz, ppv_axis_t j );
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

/* ELEMENT ENUMERATION */

typedef void ppv_op_t ( ppv_index_t *ix, ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC );
  /* Client procedure for {ppv_enum}. */

void ppv_enum ( ppv_op_t op, ppv_grid_t *A, ppv_grid_t *B, ppv_grid_t *C );
  /* Enumerates all valid samples of the arrays {A,B,C} in parallel,
    with matching indices. For each index tuple {[ix[0],..ix[N-1]]}, calls
    {op(ix, pA, pB, pC)} where {pA,pB,pC} are the
    positions of the corresponding samples in the three arrays. 
    
    If an array is NULL, it is not scanned and the corresponding {pos}
    is always 0. All {ppv_grid_t} arguments that are not NULL must
    have the same {size} vector. */

/* DESCRIPTOR PRINTOUT */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_grid_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* DESCRIPTOR VALIDATION */

bool_t ppv_descriptor_is_valid ( ppv_grid_t *A, bool_t die );
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
    {unsigned int} variable is at least 32 bits long.

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
