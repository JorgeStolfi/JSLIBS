/* Portable multi-dimensional sample arrays. */
/* Last edited on 2021-06-22 13:44:20 by jstolfi */

#ifndef ppv_array_H
#define ppv_array_H

#include <ppv_types.h>
#include <stdio.h>

/* !!! Add more parameters to {ppv_index_next,ppb_index_prev} as in {indexing.h} !!! */
/* !!! Allow {bpw} to be 64 too. !!! */
/* !!! Allow {A.d} to be zero? !!! */

typedef 
  struct ppv_array_t {
    ppv_dim_t d;        /* Number of axes (ndices) in the array. */
    ppv_size_t *size;      /* Extent of array's domain along each axis. */
    ppv_step_t *step;      /* Position increment for each index. */
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
      {d = A.d} indices {[ix[0],..ix[d-1]]}.

      The /shape/ of {A} is defined by {d} non-negative integers,
      {A.size[0..d-1]}. Index {ix[ax]} ranges from 0 to {A.size[ax]-1},
      inclusive.
      
      The descriptor operations such as {ppv_swap_indices}, {ppv_crop}
      and {ppv_slice} allow the indices to be redefined and remapped in
      various ways, at negligible cost.

    Intended use

      Although a {ppv_array_t} structure can be used to store
      arbitrary arrays of integers, it was primarily designed to be a
      unified format for uncompressed, uniformly sampled, muti-channel
      and multi-dimensional data -- audio, still image (both 2D and 3D)
      and video.
      
      For example, a color video file could be represented by a
      {ppv_array_t} {A} with {A.d = 4}, with index 0 being the
      color channel (R,G,B, etc.), indices 1 and 2 running along the
      spatial axes (width and height), and index 3 being time (the frame
      number). Thus, a 10-sec color TV movie at 60 frames per second
      could have {A.size = {3,640,480,600}}.
      
      A tridimensional video, showing the evolution of a 3D tomogram
      over time, could instead have {A.d = 5}, with indices 1,2,3
      being width, height, and depth, and index 4 being time.
      
      As another example, an audio file could have {A.d = 2} where
      index 1 identifies the sound track or stereo channel, and index 2
      is time (the sample index). Then a 4-channel 10-second audio file
      sampled at 44100 Hz would have have {A.size={4,441000}}. If that
      audio file was chopped into overlapping segments associated with
      frames video, then one could have {A.d = 3}, where indices 0 and
      1 are as above, and index 2 is the frame (segment) number. For 100
      millisecond segments of a 10-second 60 fps video, we woould have
      {A.size = {4,4410,600}}.

    Empty arrays

      If any of the sizes {size[ax]} is zero, the array has no samples
      and is said to be /empty/. An empty array still has a definite
      shape, which is significant in certain operations.

    Value range of each sample

      The samples of an array {A} have a fixed number of bits per
      sample, {A.bps}. By definition, each sample is an unsigned
      integer in the range {0..2^A.bps - 1}; in particular, if {A.bps}
      is zero, then all samples have value 0.

    Memoryless arrays

      A /memoryless/ array is one that is empty or has {bps==0},
      and therefore occupies no storage.

    Sample indexing

      Each sample of an array {A} has a /position/ which is an
      affine combination of its {d} indices {ix[0],..ix[d-1]}:

        | pos(A, ix) = A.base + ix[0]*A.step[0] + ·· ix[d-1]*A.step[d-1]

      Each coefficient {A.step[ax]} may be positive, negative, or zero.
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

ppv_array_t *ppv_array_new ( ppv_dim_t d, ppv_size_t sz[], ppv_nbits_t bps, ppv_nbits_t bpw );
  /* Allocates a descriptor {A} for a newly allocated array of samples,
    with {A.d = d}, {A.size = {sz[0],.. sz[d-1]}} and the specified number of bits
    per sample and per word. The individual sizes {sz[ax]} must not
    exceed {ppv_MAX_SIZE},
    
    If the array {A} turns out to be memoryless, {A.el} is set to NULL
    and all samples will have the same position (0).
    
    Otherwise, {A.el} will point to a newly allocated area large
    enough to contain all samples. The samples {A[ix[0],..ix[d-1]]} will
    be packed as compactly as possible, linearized in the obvious way
    with index {ix[0]} varying fastest and {ix[d-1]} the slowest. That
    is, {A.base} will be 0, and each increment {A.step[ax]} will be
    the product of all sizes {A.size[j]} with {j<ax}. 
    
    The number of axes {d} must not exceed {ppv_MAX_AXES}}. The product
    {sz[0]*..sz[d-1]} must not exceed {ppv_MAX_SAMPLES}, and the total
    storage occupied by them must not exceed {ppv_MAX_BYTES}.
    
    In any case, the call {free(A.el)} will reclaim the storage used by
    the array. The vectors {A.size} and {A.step} are allocated in the
    same {malloc} record, so one should never call {free} on them; the
    call {free(A)} will reclaim their storage too. */

ppv_array_t *ppv_array_clone ( ppv_array_t *A );
  /* Creates a copy {S} of the descriptor of the array {A} that shares
    the same sample storage as {A}. 
    
    All fields of {S} a(including {S->step} and {S->size}) will be
    separate memory areas from {A}, and will be initialized with the
    current values in {A}. However, {S} and {A} will share the same
    sample storage area (S->el == A->el). */

ppv_sample_count_t ppv_sample_count( ppv_array_t *A, bool_t reptoo );
  /* Returns the total number of sample positions in {A}. If {reptoo}
    is true, counts replicated elements as distict; that is, counts the 
    number of distinct valid index tuples.  If {reptoo} is false, counts
    those elements only once. */

/* SAMPLE EXTRACTION AND INSERTION */

ppv_sample_t ppv_get_sample ( ppv_array_t *A, ppv_index_t ix[] );
  /* Extracts the sample {A[ix[0],ix[1],.. ix[d-1]]}. */

void ppv_set_sample ( ppv_array_t *A, ppv_index_t ix[], ppv_sample_t qv );
  /* Stores value {qv} into the element {A[ix[0],ix[1],.. ix[d-1]]}. */

/* SAMPLE EXTRACTION AND INSERTION WITH POSITION */

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, ppv_index_t ix[] );
  /* Computes the position of sample {A[ix[0],.. ix[d-1]]}. Does not check
    whether the element exists or not. */

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

void ppv_index_clear ( ppv_dim_t d, ppv_index_t ix[] );
  /* Sets {ix[ax] = 0} for every axis {ax}. */

void ppv_index_assign ( ppv_dim_t d, ppv_index_t ix[], ppv_index_t val[] );
  /* Sets {ix[ax] = val[ax]} for every axis {ax}. */

void ppv_index_shift ( ppv_dim_t d, ppv_index_t ix[], ppv_index_t *inc );
  /* Sets {ix[ax] += inc[ax]} for every axis {ax}. */

sign_t ppv_index_compare ( ppv_dim_t d, ppv_index_t ixa[], ppv_index_t ixb[] );
  /* Returns {NEG}, {ZER}, or {POS} depending on whether {ixa} is less
    than, equal to, or greater than {ixb} in the C index order. */

bool_t ppv_index_first ( ppv_index_t ix[], ppv_array_t *A );
  /* Sets {ix[ax] = 0} for every axis {ax}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[ax]>0} for all {ax}). */

bool_t ppv_index_last ( ppv_index_t ix[], ppv_array_t *A );
  /* Sets {ix[ax] = A.size[ax]-1} for every axis {ax}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[ax]>0} for all {ax}). */

bool_t ppv_index_next ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t na, ppv_pos_t *p );
  /* Set {ix[0..na-1]} to the next combination of values of those
    indices that is valid for {A}, in the C index order.
  
    More precisely, the procedure scans the indices {ix[na-1]},
    {ix[na-2]}, ... {ix[0]}, looking for an index {ix[ax]} that is
    strictly less than its limit {A.size[ax]-1}. If the procedure finds
    such an index, it increments that index by one, sets every following
    index {ix[j]} with {j>ax} to 0, and returns FALSE. If every index {ix[ax]} has
    reached its limit {A.size[ax]-1}, the procedure sets all indices back
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

bool_t ppv_index_prev ( ppv_index_t ix[], ppv_array_t *A, ppv_dim_t na, ppv_pos_t *p );
  /* Like {ppv_index_next}, but enumerates the index tuples in the reverse order. 
  
    More precisely, the procedure scans the indices {ix[na-1]},
    {ix[na-2]}, ... {ix[0]}, looking for an index {ix[ax]} that is
    strictly positive. If the procedure finds such an index, it
    decrements that index by one, sets every following index {ix[j]} with {j>ax} to
    its upper limit {A.size[j]-1}, and returns FALSE. If every index
    {ix[ax]} is zero, the procedure sets every index {ix[ax]} to its upper
    limit {A.size[ax]-1}, and returns TRUE.
    
    Thus, to scan all the elements of an array {A} in reverse order, use 

      | if (ppv_index_last(ix,A))
      |   { ppv_pos_t p = ppv_position(ix,A);
      |     do { ... } while (! ppv_index_prev(ix,A,&p));
      |   }

   */
    
/* INDEX TUPLE VALIDATION */
    
bool_t ppv_index_is_valid ( ppv_index_t ix[], ppv_array_t *A );
  /* Returns TRUE iff {ix} is a valid index tuple for the array {A},
    i.e., if sample {A[ix[0],..ix[na-1]]} exists. This is true if and
    only if {ix[ax]} lies in the range {0..A.step[ax]-1}, for every
    axis {ax}. */

/* DESCRIPTOR MANIPULATION */

/* The following procedures modify the {size}, {step} and {base} field
  of an array descriptor {*A}, so as to change the set of valid indices
  and/or the mapping between indices and sample positions. No samples are
  actually allocated, reclaimed, or changed. 
  
  All these procedures preserve the validity of the descriptor (in the
  sense of {ppv_descriptor_is_valid} below) when given valid
  arguments. */

void ppv_crop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t skip, ppv_size_t  keep);
  /* Crops the array {A} to the range {skip..skip+keep-1} along axis
    {ax}. Thus {size[ax]} becomes {keep}, and using value {r} for that
    index in the new descriptor is equivalent to using {r + skip} in
    the old one.
    
    The value of {skip+keep} must not exceed the original {size[ax]}.
    {keep} may be zero, in which case the array becomes empty (but
    retains its size along the other axes) and all its steps are
    reset to zero. */

void ppv_subsample ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {ax}, so that using {r}
    for that index in the new descriptor is equivalent to using
    {r*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[ax]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[ax]}. In particular, {size[ax]} will be 0 
    if and only if it was 0 originally. */

void ppv_reverse ( ppv_array_t *A, ppv_axis_t ax );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {ax}
    runs in the opposite direction.  */

void ppv_replicate ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {ax}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {A.step[ax]} to 0,
    {A.size[ax]=sz}. On input, {size[ax]} must be 1, and {sz} must be
    positive. */

void ppv_swap_indices ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1, ppv_dim_t n );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the {n} indices that
    begin with index {ax0} and the {n} indices that begin with index
    {ax1} are exchanged. Thus, for instance,
    {ppv_swap_indices(img,1,2,1)} could be used to exchange the rows and
    columns of a 2D color image, assuming index 0 is the color channel.
    
    The two sets of indices must be either identical ({ax0==ax1}), or
    disjoint (|ax0-ax1|>=n). If they are identical, the procedure has no
    effect. */

void ppv_flip_indices ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1 );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the indices {ax0..ax1} in reverse order. Thus,
    after the call the call {ppv_flip_indices(A,2,5)}, 
    element {A[p,q,r,s,t,u]} is equivalent to {A[p,q,u,t,s,r]}
    of the original array. */

void ppv_diagonal ( ppv_array_t *A, ppv_axis_t ax0, ppv_axis_t ax1 );
  /* Modifies the descriptor {A} so that it describes a diagonal slice
    through the original array, cut parallel to the bisector of axes
    {ax0} and {ax1}; and redefines axis {ax0} to be along that diagonal.
    Thus, using index values {r} and {s} along those axes in the 
    new array is equivalent to using {r} and {s+r} in the original 
    array.
    
    The procedure requires {ax0 != ax1}, and {0 < A.size[ax0] <=
    A.size[ax1]}. It reduces {A.size[ax1]} by {A.size[ax0]-1}, and adds
    {A.step[ax1]} to {A.step[ax0]} .
    
    Thus, for example, if {A} is a color image with shape
    {3,640,480,1,1,1}, the call {ppv_diagonal(A,2,1)} will cut it down
    to a diagonal band starting at the main diagonal and extending
    {161=640-(480-1)} pixels to the left, and shear it horizontally so
    that it becomes a rectangular image with shape {3,161,480,1,1,1}.
    Then sample {A[c,h,v,0,0,0]} of the new image will be sample
    {A[c,h+v,v,0,0,0]} of the original one. */

/* OPERATORS THAT CHANGE THE DIMENSION */

ppv_array_t *ppv_slice ( ppv_array_t *A, ppv_axis_t ax, ppv_index_t ix );
  /* Returns a new descriptor {S} that describes the subset of the
    samples of {A} where the index of axis {ax} has value {ix}. The
    array {S} will have one dimension less than {A}, and the index tuple
    of each retained element will have the index corresponding to {ax}
    omitted. Thus, for example, after the call {S = ppv_slice(A,2,5)},
    element {S[p,q,s,t,u]} will be the same as {A[p,q,5,s,t,u]}.
    
    The axis {ax} must be in {0..A->d-1}, and the index {ix} must be valid for it.
    Thus the array {A} must not be empty. */

ppv_array_t *ppv_chop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz );
  /* Returns a new descriptor {S} that describes the same samples
    as {A}, with the same indices, except that axis {ax} is
    replaced by the combination of axes {ax} and a new axis {d}.
    
    Specifically, {S} will have {S.d = d+1}, where {d = A.d};
    {S.size[ax]} will be {sz}; and {S.size[d]} will be {A.size[ax]/sz}. 
    Using index values {ixi} on axis {ax} and {ixj} on axis {d} will
    yield the same element as using index value {ixi*sz +ixd} on
    axis {ax} of {A}
    
    The original size {A.size[ax]} must be a multiple of {sz}, otherwise
    the element will become inacessible.  Inparticular, if {A.size[ax]} is less than 
    {sz}, the array will become empty. ???
    
    For exmple, {K = ppv_chop(V, 0, 100)} will turn one-dimensional array {V} of 5000
    samples (along axis 0) into a two-dimensional monochrome image {K} with {sz=100} samples per
    row (along index 0) and 50 rows (along axis 1). 
    
    As another example, {S = ppv_chop(A, 1, 64)} could be used to chop a
    color image with shape {3,640,480} into 10 vertical stripes, 64
    pixels wide; and stack those stripes in the depth direction,
    yielding a color tomogram with shape {3,64,480,10}. Thus sample
    {S[c,h,v,d]} of the new array will be the same as {A[c,64*d+h,v]} of
    the original. */

/* ELEMENT ENUMERATION */

bool_t ppv_enum 
  ( ppv_index_pos3_op_t *op, 
    bool_t reverse, 
    ppv_array_t *A, 
    ppv_array_t *B, 
    ppv_array_t *C 
  );
  /* Enumerates all valid samples of the arrays {A,B,C} in parallel, with
    matching indices. For each index tuple {[ix[0],..ix[d-1]]}, calls {op(ix,
    pA, pB, pC)} where {pA,pB,pC} are the positions of the corresponding samples
    in the three arrays.
    
    The arrays may have different bits-per-pixel and bits-per-word parameters,
    as well as different position steps in any axes. If {A->bps} is greater
    than {B->bps}, the extra high-order bits of each sample of {A} are set to zeros.
    If {A->bps} is less than {B->bps}, the excess high-order bits of each sample of {B}
    are dropped. 
    
    The procedure stops the enumeraton when the call to {op} returns
    {TRUE}, and returns {TRUE}. Otherwise the enumeration continues
    until all valid index tuples are exhausted, and returns {FALSE}. In
    particular, if the array is empty (that is, one of the sizes
    {sz[ax]} is zero), the procedure returns {FALSE} without ever
    calling {op}. 
    
    The index tuples are scanned from the first one {[0,..0]} increasing if
    {reverse} is false, or from the max index tuple decreasing if {reverse} is
    true. In either case, the indices are varied in the C-like order (last index
    varies fastest).
    
    Any array that is NULL is not scanned, and the corresponding {pos}
    is always 0. All {ppv_array_t} arguments that are not NULL must
    have the same {size} vector. */

void ppv_array_assign ( ppv_array_t *A, ppv_array_t *B );
  /* Copies the samples of {B} into the corresponding samples of {A}.
    The two arrays must have the same size along every axis. 
    
    The assgnments are done in lex order of the indices. This matters
    only if {A} is replicated along one or more axes whereas {B} is not;
    so that multiple elements of {B} with different values may be
    assigned to the same element of {A}. */

/* DESCRIPTOR PRINTOUT */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* DESCRIPTOR VALIDATION */

bool_t ppv_descriptor_is_valid ( ppv_array_t *A, bool_t die );
  /* Checks whether the descriptor {A} seems to be valid.
    
    To be valid, an array descriptor must satisfy the following
    conditions:

      (0) The dimension {A.d} must be in {0..ppv_MAX_DIM}.

      (1) Two samples with distinct valid index tuples have the same
        position only if {step[ax]} is zero for every axis where
        their indices differ.

      (2) For any {ax} in {0..A.d-1}, if {A.size[ax]} is 1, then {A.step[ax]} is zero.

      (3) If the array is memoryless, then {A.base} and all {A.step}
        coefficients are zero, and {A.el} is NULL.

      (4) The field {A.npos} is the number of samples with distinct
        positions. Because of (3), {A.npos} is zero if any {A.size[ax]}
        is zero, otherwise it is the product of {A.size[ax]}
        for every {ax} where {A.step[ax]} is nonzero.

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
