/* Portable 6-dimensional sample arrays. */
/* Last edited on 2005-05-29 10:32:28 by stolfi */

#ifndef ppv_array_H
#define ppv_array_H

#include <ppv_types.h>
#include <stdio.h>

typedef 
  struct ppv_array_t {
    ppv_size_t size[6];  /* Extent of array's domain along each axis. */
    ppv_step_t step[6];  /* Position increment for each index. */
    ppv_pos_t base;      /* Position of element {[0,0,0,0,0,0]}. */
    ppv_tot_size_t npos; /* Number of distinct sample positions. */
    ppv_nbits_t bps;     /* Bits per sample. */
    ppv_nbits_t bpw;     /* Bits per word (storage unit). */
    void *el;            /* Start of storage area (NULL if no samples). */
  } ppv_array_t; 
  /* 
    A {ppv_array_t} is descriptor for a five-dimensional array of
    small non-negative integers (/samples/), in a packed binary
    format.
    
    Intended use

      Although a {ppv_array_t} structure can be used to store
      arbitrary arrays of integers, it was primarily designed to be a
      unified format for uncompressed, uniformly sampled, muti-channel
      and multi-dimensional data -- audio, still image (both 2D and
      3D) and video. In particular, it was felt that six index
      dimensions were necessary and sufficient to satisfy most of
      those applications (e.g. channel + x + y + z + time).

    Array shape

      The /shape/ of an array {A} is defined by six non-negative
      integers, {A.size[0..5]}. Each sample of {A} is identified by six
      indices {[i[0],..i[5]]}, where index {i[ax]} ranges from 0 to
      {A.size[ax]-1}, inclusive.

      Thus, a 1-sec color TV movie could have {size={3,640,480,1,30,1}},
      whereas the corresponding 44100 Hz, four-channel, frame-chopped
      audio file could have {size={4,29400,1,1,30,1}}.

    Empty arrays

      If any {size[ax]} is zero, the array has no samples
      and is said to be /empty/. An empty array still
      has a definite shape, which is significant in certain 
      operations.

    Value range of each sample

      The samples of an array {A} have a definite number of bits per
      sample, {A.bps}. By definition, each sample is an integer in the range
      {0..2^A.bps - 1}; in particular, if {bps} is zero, then all
      samples are zero.

    Memoryless arrays

      A /memoryless/ array that is one that is empty or has {bps==0},
      and therefore occupies no storage.

    Sample indexing

      Each sample of an array {A} also has a /position/ which is an
      affine combination of its six indices:

        | pos(A, p,q,r,s,t,u) = A.base + 
            p*A.step[0] + q*A.step[1] + r*A.step[2] + 
            s*A.step[3] + t*A.step[4] + u*A.step[5]

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
      one word may contain two or more samples. For further details on
      how samples are stored, see the comments in {ppv_packing.h}.  */

/* INDEX TUPLE MANIPULATION */

void ppv_index_clear ( ppv_index_t *ix );
  /* Sets {ix[ax] = 0} for every axis {ax}. */

void ppv_index_assign ( ppv_index_t *ix, ppv_index_t *val );
  /* Sets {ix[ax] = val[ax]} for every axis {ax}. */

void ppv_index_shift ( ppv_index_t *ix, ppv_index_t *inc );
  /* Sets {ix[ax] += inc[ax]} for every axis {ax}. */

sign_t ppv_index_cmp ( ppv_index_t *ixa, ppv_index_t *ixb );
  /* Returns {NEG}, {ZER}, or {POS} depending on whether 
    {ixa} is less than, equal to, or greater than {ixb} 
    in lexicographic order. */
    
/* INDEX TUPLE VALIDATION */
    
bool_t ppv_inside ( ppv_array_t *A, ppv_index_t *ix );
  /* Returns TRUE iff {ix} is a valid index tuple for the array {A},
    i.e., if sample {A[ix[0],..ix[5]]} exists. This is true if and
    only if {ix[ax]} lies in the range {0..A.step[ax]-1}, for every
    axis {ax}. */

/* DESCRIPTOR VALIDATION */

bool_t ppv_valid ( ppv_array_t *A );
  /* Returns TRUE iff the descriptor {A} seems to be valid.
    
    To be valid, an array descriptor must satisfy the following
    conditions:

      (1) Two samples with distinct valid index tuples have the same
        position only if {step[ax]} is zero for every axis where
        their indices differ.

      (2) For any {ax}, if {size[ax]} is 1, then {step[ax]} is zero.

      (3) If the array is memoryless, then {base} and all {step}
        coefficients are zero, and {el} is NULL.

      (4) The field {npos} is the number of samples with distinct
        positions. Because of (3), {npos} is zero if any {size[ax]}
        is zero, otherwise it is the product of {size[ax]}
        for every {ax} where {step[ax]} is nonzero.

      (5) The storage area of every valid sample is contained
        in the allocated area pointed by {el}. (Note that 
        when {bps==0} this is true even if {el==NULL}.)
        
    Condition (1) is hard to check exactly, so this procedure checks
    instead a more stringent condition: (1') there is a permutation of
    the indices such that increasing lexicographic order of the
    permuted index tuples implies increasing sample positions.
  */

/* ARRAY ALLOCATION */

ppv_array_t ppv_new_array ( ppv_size_t *sz, ppv_nbits_t bps, ppv_nbits_t bpw );
  /* Returns a valid descriptor {A} for a newly allocated array of
    samples, with {A.size = {sz[0],.. sz[5]}} and the specified number
    of bits per sample and per word. The individual sizes {sz[ax]}
    must not exceed {ppv_MAX_SIZE},
    
    If the array {A} turns out to be memoryless, {A.el} is set to NULL
    and all samples will have the same position (0).
    
    Otherwise, {A.el} will point to a newly allocated area large
    enough to contain all samples. The samples {A[ix[0],..ix[5]]} will
    be packed as compactly as possible, linearized in the obvious way
    with index {ix[0]} varying fastest and {ix[5]} the slowest. That
    is, {A.base} will be 0, and each increment {A.step[ax]} will be
    the product of all sizes {A.size[bx]} with {bx<ax}. The product
    {sz[0]*..sz[5]} must not exceed {ppv_MAX_SAMPLES}, and the total
    storage occupied by them must not exceed {ppv_MAX_BYTES}.
    
    In any case, the call {free(A.el)} will reclaim the storage used
    by the array. */

/* SAMPLE INDEXING */

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, ppv_index_t *ix );
  /* Computes the position of sample {A[ix[0],.. ix[5]]}. Does not check
    whether the element exists or not. */

/* SAMPLE EXTRACTION AND INSERTION */

/* For the procedures below, {pos} is the position of a sample in a sample
  storage area with address {el}. The samples are assumed to be
  packed with {bps} bits per sample and {bpw} bits per word. The user
  must make sure that {pos} is valid (i.e. the sample actually
  exists.) */

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

/* DESCRIPTOR MANIPULATION */

/* The following procedures modify the {size}, {step} and {base} field
  of a {ppv_array_t}, so as to change the set of valid indices and the
  mapping between indices and sample positions. No samples are
  actually allocated, reclaimed, or changed. All these procedure below
  return valid descriptors when given valid arguments. */

void ppv_crop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t skip, ppv_size_t sz );
  /* Crops the descriptor {A} along axis {ax}, so that {size[ax]}
    becomes {sz}, and using {k} for that index in the new descriptor is
    equivalent to using {skip + k} in the old one.
    
    The value of {skip+sz} must not exceed the original {size[ax]}.
    {sz} may be zero, in which case the array becomes empty (but retains
    its size along the other axes). */

void ppv_subsample ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {ax}, so that using {k}
    for that index in the new descriptor is equivalent to using
    {k*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[ax]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[ax]}. In particular, {size[ax]} will be 0 
    if and only if it was 0 originally. */

void ppv_transpose ( ppv_array_t *A, ppv_axis_t ax, ppv_axis_t bx );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that indices {ax} and
    {bx} are exchanged. Thus, for instance, {transpose(img, 1,2)} 
    could be used to exchange the rows and columns of a 2D color image.
    If {ax == bx}, the procedure has no effect. */

void ppv_flip ( ppv_array_t *A, ppv_axis_t ax );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {ax}
    runs in the opposite direction.  */

void ppv_replicate ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {ax}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {step[ax]} to 0,
    {size[ax]=sz}. On input, {size[ax]} must be 1, and {sz} must be
    positive. */

void ppv_chop ( ppv_array_t *A, ppv_axis_t ax, ppv_size_t sz, ppv_axis_t bx );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that index {ax} is
    replaced by the combination of indices {ax} and {bx} (which must
    be different). Upon entry, {size[bx]} must be 1. Upon exit,
    {size[ax]} will be {sz} and {size[bx]} will be {oldsz/sz} where {oldsz} is
    the original value of {size[ax]}. 
    
    If {oldsz} is not a multiple of {sz}, then the remainder elements
    will become inacessible.  Inparticular, if {oldsz} is less than 
    {sz}, the array will become empty.
    
    For exmple, {ppv_chop(V, 0, 100, 1)} will turn a vector of 5000
    samples (along axis 0) into a greyscale image with 100 samples per
    row (along index 0) and 50 rows (along axis 1). As another
    example, {ppv_chop(A, 1, 64, 3)} could be used to chop a color
    image, 640 pixels wide, into 10 vertical stripes, 64 pixels wide;
    and stack those stripes in the depth direction like the slices of
    a tomograph. */

void ppv_diagonal ( ppv_array_t *A, ppv_axis_t ax, ppv_axis_t bx );
  /* Modifies the descriptor {A} so that it describes a diagonal slice
    through the original array, bisecting axes {ax} and {bx}. Namely,
    {size[ax]} is set to the minimum of {{size[ax],size[bx]}};
    {size[bx]} is set to 1; and using values {k} and {0} for those two
    indices, in the new array, yields the same sample as using {k} for
    both indices in the old array. If {ax == bx}, the procedure has no
    effect.
    
    Thus, for example, if {A} is a color image with shape
    {3,640,480,1,1}, the call {ppv_diagonal(A,1,2)} will turn it into
    an array with shape {3,480,1,1,1} consisting of the 480 color
    pixels along the main diagonal of that image; i.e. sample
    {A[c,i,0,0,0]} of the new array is sample {A[c,i,i,0,0]} of the
    original one. */

/* ELEMENT ENUMERATION */

typedef void ppv_op_t ( ppv_index_t *ix, ppv_pos_t posA, ppv_pos_t posB, ppv_pos_t posC );
  /* Client procedure for {ppv_enum}. */

void ppv_enum ( ppv_array_t *A, ppv_array_t *B, ppv_array_t *C, ppv_op_t op );
  /* Enumerates all valid samples of the arrays {A,B,C} in parallel,
    with matching indices. For each index tuple {[ix[0],..ix[5]]}, calls
    {op(ix, posA, posB, posC)} where {posA,posB,posC} are the
    positions of the corresponding samples in the three arrays. 
    
    If an array is NULL, it is not scanned and the corresponding {pos}
    is always 0. All {ppv_array_t} arguments that are not NULL must
    have the same {size} vector. */

/* DESCRIPTOR PRINTOUT */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_array_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* ERROR MESSAGES */

void ppv_progerror 
  ( const char *msg, 
    const char *file, 
    const unsigned int line, 
    const char* proc 
  );
  /* Prints {file ":" line ": (" *proc ")" *msg} to {stderr} and stops. */

#define ppv_error(msg) \
  ppv_progerror((msg), __FILE__, __LINE__, __FUNCTION__)
  /* Prints {*msg} and the source location.  */

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
