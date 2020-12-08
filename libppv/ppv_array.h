/* Portable multi-dimensional sample arrays. */
/* Last edited on 2016-03-16 15:52:53 by stolfilocal */

#ifndef ppv_array_H
#define ppv_array_H

#include <ppv_types.h>
#include <stdio.h>

/* !!! Add more parameters to {ppv_index_next,ppb_index_prev} as in {indexing.h} !!! */
/* !!! Add an explicit number of indices {A.na} as in {array.h}. !!! */
/* !!! Allow {bpw} to be 64 too. !!! */
/* !!! Make the number of axes variable, with {size} and {step} allocated after the end of record. !!! */

#define ppv_array_NAXES 6
  /* Number of indices in a {ppv_array_t}. */

typedef 
  struct ppv_array_t {
    ppv_size_t size[ppv_array_NAXES];  /* Extent of array's domain along each axis. */
    ppv_step_t step[ppv_array_NAXES];  /* Position increment for each index. */
    ppv_pos_t base;        /* Position of element {[0,..0]}. */
    ppv_nbits_t bps;       /* Bits per sample. */
    ppv_nbits_t bpw;       /* Bits per word (storage unit). */
    void *el;              /* Start of storage area (NULL if no samples). */
  } ppv_array_t; 
  /* 
    A {ppv_array_t} is descriptor for a multi-dimensional array of
    small non-negative integers (/samples/), stored in memory 
    in a packed binary format. 
    
    Array indexing and shape
    
      Each sample of a {ppv_array_t} {A} is identified by
      {N=ppv_NAXES} indices {[ix[0],..ix[N-1]]}.

      The /shape/ of {A} is defined by {N} non-negative integers,
      {A.size[0..N-1]}. Index {ix[i]} ranges from 0 to {A.size[i]-1},
      inclusive

    Intended use

      Although a {ppv_array_t} structure can be used to store
      arbitrary arrays of integers, it was primarily designed to be a
      unified format for uncompressed, uniformly sampled, muti-channel
      and multi-dimensional data -- audio, still image (both 2D and
      3D) and video. 
      
      It was felt that 6 index dimensions were necessary and
      sufficient to satisfy most of those applications.
      Client programs are free to assign these six dimensions as it
      suits them best; the descriptor operations such as
      {ppv_swap_indices} allow the indices to be remapped at negligible
      cost. Clients that need less than {N} dimensions should set the
      {A.size[i]=1} along any unneeded axis {i}.
      
      Most examples below assume that the arrays are used to store
      `animated color tomograms', with index 0 being the color channel
      (R,G,B, etc.) or audio channel (L,R, etc.), indices 1,2,3
      running along the spatial axes (width, height, and depth), index
      4 being time, and index 5 being unused. Thus, a 1-sec color TV
      movie could have shape {3,640,480,1,30,1}, whereas the
      corresponding 44100 Hz, four-channel, frame-chopped audio file
      could have {size={4,29400,1,1,30,1}}.

    Empty arrays

      If any of the sizes {size[i]} is zero, the array has no samples
      and is said to be /empty/. An empty array still has a definite
      shape, which is significant in certain operations.

    Value range of each sample

      The samples of an array {A} have a fixed number of bits per
      sample, {A.bps}. By definition, each sample is an unsigned
      integer in the range {0..2^A.bps - 1}; in particular, if {A.bps}
      is zero, then all samples have value 0.

    Memoryless arrays

      A /memoryless/ array that is one that is empty or has {bps==0},
      and therefore occupies no storage.

    Sample indexing

      Each sample of an array {A} has a /position/ which is an
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

ppv_array_t ppv_new_array ( ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw );
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

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, ppv_index_t ix[] );
  /* Computes the position of sample {A[ix[0],.. ix[N-1]]}. Does not check
    whether the element exists or not. */

/* SAMPLE EXTRACTION AND INSERTION */

/* For the procedures below, {pos} is the position of a sample in a
  sample storage area with address {el}. The samples are assumed to be
  packed with {bps} bits per sample and {bpw} bits per word; see the
  section "DETAILS OF SAMPLE PACKING" below. The user must make sure
  that {pos} is valid (i.e. the sample actually exists.)  */

ppv_sample_t ppv_get_sample 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos 
  );
  /* Extracts the sample with position {pos}. */

void ppv_set_sample 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos, 
    ppv_sample_t qv 
  );
  /* Stores value {qv} into the sample with position {pos}. */

/* INDEX TUPLE MANIPULATION */

void ppv_index_clear ( ppv_index_t ix[] );
  /* Sets {ix[i] = 0} for every axis {i}. */

void ppv_index_assign ( ppv_index_t ix[], ppv_index_t *val );
  /* Sets {ix[i] = val[i]} for every axis {i}. */

bool_t ppv_index_first ( ppv_index_t ix[], ppv_array_t *A );
  /* Sets {ix[i] = 0} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

bool_t ppv_index_last ( ppv_index_t ix[], ppv_array_t *A );
  /* Sets {ix[i] = A.size[i]-1} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

void ppv_index_shift ( ppv_index_t ix[], ppv_index_t *inc );
  /* Sets {ix[i] += inc[i]} for every axis {i}. */

sign_t ppv_index_compare ( ppv_index_t ixa[], ppv_index_t ixb[] );
  /* Returns {NEG}, {ZER}, or {POS} depending on whether {ixa} is less
    than, equal to, or greater than {ixb} in the C index order. */

bool_t ppv_index_next ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p );
  /* Set {ix[0..d-1]} to the next combination of values of those
    indices that is valid for {A}, in the C index order.
  
    More precisely, the procedure scans the indices {ix[na-1]},
    {ix[na-2]}, ... {ix[0]}, looking for an index {ix[i]} that is
    strictly less than its limit {A.size[i]-1}. If the procedure finds
    such an index, it increments that index by one, sets every following
    index {ix[j]} with {j>i} to 0, and returns FALSE. If every index {ix[i]} has
    reached its limit {A.size[i]-1}, the procedure sets all indices back
    to {(0,..0)}, and returns TRUE.
    
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

bool_t ppv_index_prev ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t d, ppv_pos_t *p );
  /* Like {ppv_index_next}, but enumerates the index tuples in the reverse order. 
  
    More precisely, the procedure scans the indices {ix[na-1]},
    {ix[na-2]}, ... {ix[0]}, looking for an index {ix[i]} that is
    strictly positive. If the procedure finds such an index, it
    decrements that index by one, sets every following index {ix[j]} with {j>i} to
    its upper limit {A.size[j]-1}, and returns FALSE. If every index
    {ix[i]} is zero, the procedure sets every index {ix[i]} to its upper
    limit {A.size[i]-1}, and returns TRUE.
    
    Thus, to scan all the elements of an array {A} in reverse order, use 

      | if (ppv_index_last(ix,A))
      |   { ppv_pos_t p = ppv_position(ix,A);
      |     do { ... } while (! ppv_index_prev(ix,A,&p));
      |   }

   */
    
/* INDEX TUPLE VALIDATION */
    
bool_t ppv_index_is_valid ( ppv_index_t ix[], ppv_array_t *A );
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

void ppv_crop ( ppv_array_t *A, ppv_axis_t i, ppv_size_t skip, ppv_size_t  keep);
  /* Crops the array {A} to the range {skip..skip+keep-1} along axis
    {i}. Thus {size[i]} becomes {keep}, and using value {r} for that
    index in the new descriptor is equivalent to using {r + skip} in
    the old one.
    
    The value of {skip+keep} must not exceed the original {size[i]}.
    {keep} may be zero, in which case the array becomes empty (but
    retains its size along the other axes) and all its steps are
    reset to zero. */

void ppv_subsample ( ppv_array_t *A, ppv_axis_t i, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {i}, so that using {r}
    for that index in the new descriptor is equivalent to using
    {r*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[i]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[i]}. In particular, {size[i]} will be 0 
    if and only if it was 0 originally. */

void ppv_flip ( ppv_array_t *A, ppv_axis_t i );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {i}
    runs in the opposite direction.  */

void ppv_replicate ( ppv_array_t *A, ppv_axis_t i, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {i}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {A.step[i]} to 0,
    {A.size[i]=sz}. On input, {size[i]} must be 1, and {sz} must be
    positive. */

void ppv_swap_indices ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j, ppv_dim_t n );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the {n} indices that begin
    with index {i} and the {n} indices that begin with index {j}
    are exchanged.  Thus, for instance, {ppv_swap_indices(img, 1,2, 1)} 
    could be used to exchange the rows and columns of a 2D color image,
    assuming index 0 is the color channel.
    
    The two sets of indices must be either identical ({i==j}, or {n==0}, in
    which case the procedure has no effect) or disjoint. */

void ppv_flip_indices ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the indices {i..j} in reverse order. Thus,
    after the call the call {ppv_flip_indices(A,2,5)}, 
    element {A[p,q,r,s,t,u]} is equivalent to {A[p,q,u,t,s,r]}
    of the original array. */

void ppv_slice ( ppv_array_t *A, ppv_dim_t n, ppv_axis_t ax[], ppv_index_t ix[] );
  /* Modifies the descriptor {A} so that it describes the subset of samples
    where the {n} indices {ax[0..n-1} are fixed to the values
    {ix[0..n-1]}, and then all indices are rearranged so that the affected ones
    (which now have size 1) end up at the end. Thus, for example, after the call
    {ppv_slice(A,2,{2,4},{5,7})}, element {A[p,q,t,u,0,0]} is equivalent to
    {A[p,q,5,t,7,u]} of the original array. 
    
    The indices {ax[0..d-1]} must be in {0..N-1} and strictly increasing. If
    {ax} is null, the procedure assumes {ax[i] == i} for all {i}, that is, sets
    and rearranges the first {n} indices of {A}. The procedure is a no-op if
    {n==0}; otherwise, each index {ix[k]} must be valid, so the array must not
    be empty. */

void ppv_diagonal ( ppv_array_t *A, ppv_axis_t i, ppv_axis_t j );
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

void ppv_chop ( ppv_array_t *A, ppv_axis_t i, ppv_size_t sz, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that index {i} is
    replaced by the combination of indices {i} and {j} (which must
    be different). Upon entry, {A.size[j]} must be 1. Upon exit,
    {A.size[i]} will be {sz} and {A.size[j]} will be {oldsz/sz} where {oldsz} is
    the original value of {A.size[i]}. 
    
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

typedef ix_op_t ppv_op_t;
  /* Client procedure for {ppv_enum}. */

void ppv_enum 
  ( ppv_op_t *op, 
    bool_t reverse, 
    ppv_array_t *A, 
    ppv_array_t *B, 
    ppv_array_t *C 
  );
  /* Enumerates all valid samples of the arrays {A,B,C} in parallel, with
    matching indices. For each index tuple {[ix[0],..ix[N-1]]}, calls {op(ix,
    pA, pB, pC)} where {pA,pB,pC} are the positions of the corresponding samples
    in the three arrays.
    
    The index tuples are scanned from the first one {[0,..0]} increasing if
    {reverse} is false, from the max index tuple decreasing if {reverse} is
    true. In either case, the indices are varied in the C-like order (last index
    varies fastest).
    
    Any array that is NULL is not scanned, and the corresponding {pos}
    is always 0. All {ppv_array_t} arguments that are not NULL must
    have the same {size} vector. */

/* DESCRIPTOR PRINTOUT */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* DESCRIPTOR VALIDATION */

bool_t ppv_descriptor_is_valid ( ppv_array_t *A, bool_t die );
  /* Checks whether the descriptor {A} seems to be valid.
    
    To be valid, an array descriptor must satisfy the following
    conditions:

      (1) Two samples with distinct valid index tuples have the same
        position only if {step[i]} is zero for every axis where
        their indices differ.

      (2) For any {i}, if {A.size[i]} is 1, then {A.step[i]} is zero.

      (3) If the array is memoryless, then {A.base} and all {A.step}
        coefficients are zero, and {A.el} is NULL.

      (4) The field {A.npos} is the number of samples with distinct
        positions. Because of (3), {A.npos} is zero if any {A.size[i]}
        is zero, otherwise it is the product of {A.size[i]}
        for every {i} where {A.step[i]} is nonzero.

      (5) The storage area of every valid sample is contained
        in the allocated area pointed by {A.el}. (Note that 
        when {A.bps==0} this is true even if {A.el==NULL}.)
        
    The procedure returns TRUE if the parameters are valid.
    Otherwise, it either aborts the program with an error message 
    (if {die=TRUE}) or returns FALSE (if {die=FALSE}).
  */

/* DETAILS OF SAMPLE PACKING

  General principles

    The samples of an array {M} are stored in a vector of words,
    pointed to by {M.el}. Each word is an unsigned binary integer with
    {M.bpw} bits.

    A sample is never split across word boundaries; so, if {M.bps} is
    not an exact divisor or multiple of {M.bpw}, any incompletely filled
    words will be padded with zeros at the `big' end.
    
    For this implementation, the word size {M.bpw} must be 8, 16, or
    32. The implementation below assumes that words can be directly
    addressed by a pointer of the appropriate type, and that an
    {unsigned int} variable is at least 32 bits long.

  Packing of "small" samples

    If {M.bps} is less than {M.bpw}, then {K=M.bpw/M.bps} samples are stored
    into each word. Within each word, samples with consecutive
    positions are packed together in big-endian order. Thus, a sample
    whose position is {pos} will be stored in word {iw = pos/K}, and
    its units bit will be bit {ib = (K-1 - pos%K)*M.bps} of that word.

    For example, consider an image with {I.bps=6} and 5 samples
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

    If {M.bps} is greater than {M.bpw}, then each sample will occupy 
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

    The sample storage area address {M.el} is the address of the word
    that contains the highest-order bit of the sample with position 0.

*/

#endif
