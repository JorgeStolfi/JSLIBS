/* Portable multi-dimensional sample arrays. */
/* Last edited on 2021-07-09 22:31:31 by jstolfi */

#ifndef ppv_array_H
#define ppv_array_H

#include <ppv_types.h>
#include <stdio.h>

/* !!! Add more parameters to {ppv_index_next,ppb_index_prev} as in {indexing.h} !!! */
/* !!! Allow {bpw} to be 64 too. !!! */

typedef 
  struct ppv_array_t {
    ppv_dim_t d;           /* Number of axes (ndices) in the array. */
    ppv_size_t *size;      /* Extent of array's domain along each axis. */
    ppv_step_t *step;      /* Position increment for each index. */
    ppv_pos_t base;        /* Position of element {[0,..0]}. */
    ppv_nbits_t bps;       /* Bits per element. */
    ppv_nbits_t bpw;       /* Bits per word (storage unit). */
    ppv_sample_t maxsmp;   /* Maximum allowed sample value. */
    void *el;              /* Start of storage area (NULL if memoryless). */
  } ppv_array_t; 
  /* 
    A {ppv_array_t} is descriptor for a multi-dimensional array of
    small non-negative integers (/samples/), stored in memory 
    in a packed binary format. 
    
    Array indexing and shape
    
      Each element of a {ppv_array_t} {A} is identified by
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
      is time (the audio sample index). Then a 4-channel 10-second audio file
      sampled at 44100 Hz would have have {A.size={4,441000}}. If that
      audio file was chopped into overlapping segments associated with
      frames video, then one could have {A.d = 3}, where indices 0 and
      1 are as above, and index 2 is the frame (segment) number. For 100
      millisecond segments of a 10-second 60 fps video, we woould have
      {A.size = {4,4410,600}}.

    Empty arrays

      If any of the sizes {size[ax]} is zero, the array has no elements
      and is said to be /empty/. An empty array still has a definite
      shape, which is significant in certain operations.

    Value range of each sample

     Each element of an array {A} is stored in a fixed number of bits,
     the /bits per sample/ attribute {A.bps}. 
     
     By definition, the sample stored in each element is an unsigned
     integer in the range {0..A.maxsmp}, and {A.maxsmp} is at most
     {2^A.bps-1}. In particular, if {A.maxsmp} is zero, then each
     element can assume only one value, namely zero; and if {A.bps} is
     zero, then {A.maxsmp} is zero.

    Memoryless arrays

      A /memoryless/ array is one that is empty or has {bps} equal to
      zero, and therefore occupies no storage.
      
    Element indexing

      Each element of an array {A} has a /position/ which is an
      affine combination of its {d} indices {ix[0],..ix[d-1]}:

        | pos(A, ix) = A.base + ix[0]*A.step[0] + ·· ix[d-1]*A.step[d-1]

      Each coefficient {A.step[ax]} may be positive, negative, or zero.
      In particular, if {A.size[ax]} is 1 or 0, the corresponding
      coefficient {A.step[ax]} is irrelevant and should be zero.
      
    Replicated elements
    
      If a coefficient {A.step[ax]} is zero but {A.size[ax]} is 2 or
      more, two index vectors that differ only in that axis will yield
      the same position, hence will denote the same memory area (if
      any). The array is then said to be /replicated/ along that axis.
      In that case, there will be a dfference between the number of
      valid index tuples and the number of distinct element positions.
      
      On the other hand, two valid index tuples that differ in at least 
      one axis that is not replicated will have distinct positions and
      will be stored in disjoint memory areas.

    Element storage

      If an array {A} is not memoryless, its elemnts are stored in the
      area pointed by {A.el}. The position {pos} of an element, together
      with the the packing parameters {A.bps} and {A.bpw}, uniquely
      determines the location of the element within that area. 
      See {ppv_packing_INFO} below.
      
    Modifying attributes
    
      The attributes of an array descriptor {A} may be set or modified
      by users. However, care must be taken to preserve the properties
      stated above. In particular, unchecked out-of-bounds memory access
      may occur if the fields {A.bps}, {A.bpw}, {A.size} or {A.step} are
      improperly modified. Also, changing {A.maxsmp} from zero to a
      non-zero value requres allocating a storage area for {A.el}.
      
    Significant changes
    
      The attrbute {A.maxsmp} was added on 2021-07-08.  Several
      functions were modified to use {A.maxsmp}  instead of 
      the previously assumed maximum {2^A.bps-1}.
      
      As of 2021-07-08, the procedures {ppv_set_sample} and {ppv_get_sample}
      will check wether the index tuple is valid and the sample value
      is in {0..A.maxsmp}.
      
    */

/* ARRAY ALLOCATION */

ppv_array_t *ppv_array_new ( ppv_dim_t d, ppv_size_t sz[], ppv_sample_t maxsmp );
  /* Allocates a descriptor {A} for a newly allocated array of samples,
    with {A.d = d}, {A.size = {sz[0],.. sz[d-1]}}, and the specified max sample
    size {A.maxsmp}.  
    
    The number of axes {d} cannot exceed {ppv_MAX_AXES}, 
    the individual sizes {sz[ax]} must not exceed {ppv_MAX_SIZE},
    their product {sz[0]*..sz[d-1]} must not exceed {ppv_MAX_SAMPLES},
    and {maxsmp} must not exceed {ppv_MAX_SAMPLE_VAL}.
    
    The number of bits per sample {A.bps} the smallest number that can
    represent any value in {0..maxsmp}; that is, such that {maxsmp <
    2^A.bps}. In particular, if {maxsmp} is zero, {A.bps} will be zero
    and the array will be memoryless. 
    
    The bits per word {A.bpw} will be chosen with attempt to minimize
    the wasted space and/or maximize the packing and unpacking speed.
    
    If the array {A} turns out to be memoryless, {A.el} is set to NULL
    and all elements will have the same position (0).
    
    Otherwise, {A.el} will point to a newly allocated area large enough
    to contain all elements. The total storage occupied by them must not
    exceed {ppv_MAX_BYTES}. The elements {A[ix[0],..ix[d-1]]} will be
    packed as compactly as possible, linearized in the obvious way with
    index {ix[0]} varying fastest and {ix[d-1]} the slowest. That is,
    {A.base} will be 0, and each increment {A.step[ax]} will be the
    product of all sizes {A.size[j]} with {j<ax}.
    
    In any case, the calls {free(A.el); free(A)} will reclaim the
    storage used by the array. The vectors {A.size} and {A.step} are
    allocated in the same {malloc} record, so one should never call
    {free} on them; the call {free(A)} will reclaim their storage
    too. */

ppv_array_t *ppv_array_clone ( ppv_array_t *A );
  /* Creates a copy {S} of the descriptor of the array {A} that shares
    the same element storage as {A}. 
    
    All fields of {S} a(including {S->step} and {S->size}) will be
    separate memory areas from {A}, and will be initialized with the
    current values in {A}. However, {S} and {A} will share the same
    element storage area (S->el == A->el). */

ppv_sample_count_t ppv_compute_npos_steps ( ppv_dim_t d, ppv_size_t size[], ppv_step_t step[] );
  /* If {step} is not {NULL}, computes suitable position increments {step[0..d-1]} for 
    an array with {size[ax]} elements along each axis {ax}. In any case, returns the
    total number of elements in that array.

    Assumes that the array will have no replicated elements; that is, {step[ax]}
    will be zero only if the array is empty or {size[ax]} is 1. */

ppv_sample_count_t ppv_sample_count ( ppv_array_t *A, bool_t reptoo );
  /* Returns the total number of elements in {A}. 
  
    If {reptoo} is true, counts replicated elements as distict; that is,
    returns the count of distinct *valid index tuples*. If {reptoo} is
    false, counts those elements only once; that is, returns the count
    of distinct *element positions*. */

bool_t ppv_is_empty ( ppv_array_t *A );
  /* Returns {TRUE} if and only if the array {A} has no elements;
    that is, iff one of the sizes {A->size[0..d]} is zero.  
    
    Note that this is not the same as {A} being memoryless ({A->el ==
    NULL}), since this also includes the case {A->bps == 0}. */

ppv_sample_t ppv_max_sample( ppv_nbits_t bps );
  /* Returns the max sample value that can be stored in {bps}
    bits; that is, {2^bps-1}. */

ppv_nbits_t ppv_min_bps( ppv_sample_t maxsmp );
  /* Returns a value of bits-per-sample {bps} that is just large
    enough to hold samples in the range {0..maxsmp}.
    Iin particular, returns 0 if {maxsmp} is zero. */

ppv_nbits_t ppv_best_bpw( ppv_nbits_t bps );
  /* Returns a value of bits-per-word {bpw} that minimizes wasted
    space for a sample array with {bps} bits per sample. */

/* SAMPLE EXTRACTION AND INSERTION */

ppv_sample_t ppv_get_sample ( ppv_array_t *A, const ppv_index_t ix[] );
  /* Given an index tuple {ix[0..A.d-1]}, extracts the sample stored in
    element {A[ix]}. Fails if {ix} is not a valid index tuple for {A}, 
    or if {smp} is not in {0..A.maxsmp}.
    
    In particular, if {A.maxsmp} is zero, the index checking is done,
    but a zero sample value is returned without accessing the (non-existent)
    storage area. */

void ppv_set_sample ( ppv_array_t *A, const ppv_index_t ix[], ppv_sample_t smp );
  /* Given an index tuple {ix[0..A.d-1]} and a sample value {smp},
    stores {smp} into the element {A[ix]}.  Fails if the index tuple
    is not valid for {A}, or if {smp} is not in {0..A.maxsmp}.
    
    In particular, if {A.maxsmp} is zero, the index checking is done,
    but {smp} must be zero and it is not acually written to memory. */

bool_t ppv_index_is_valid ( const ppv_index_t ix[], ppv_array_t *A );
  /* Returns TRUE iff {ix[0..A.d-1]} is a valid index tuple for the
    array {A}. This is true if and only if {ix[ax]} lies in the range
    {0..A.step[ax]-1}, for every axis {ax}.
    
    If {ix} is valid, the element {A[ix]} nominally "exists" and has a
    value and position, but will not occupy any memory if {A.maxsmp} is
    zero. */

void ppv_array_assign ( ppv_array_t *A, ppv_array_t *B );
  /* Copies the samples of {B} into the corresponding samples of {A}.
    The two arrays must have the same size along every axis,
    and any samples of {B} must be in the range {0..A.maxsmp}
    
    The assgnments are done in lex order of the indices. This matters
    only if {A} is replicated along one or more axes whereas {B} is not;
    so that multiple elements of {B} with different values may be
    assigned to the same element of {A}. */

/* SAMPLE EXTRACTION AND INSERTION WITH POSITION */

ppv_pos_t ppv_sample_pos ( ppv_array_t *A, const ppv_index_t ix[] );
  /* Given an index tuple {ix[0..A.d-1]}, computes the position 
    of sample {A[ix]}. Does NOT check whether the index tuple is
    valid. */

/* For the procedures below, {pos} is the position of a sample in a
  sample storage area with address {el}. The samples are assumed to be
  packed with {bps} bits per sample and {bpw} bits per word; see the
  section "DETAILS OF SAMPLE PACKING" below. The user must make sure
  that {pos} is valid (i.e. the sample actually exists.) */

ppv_sample_t ppv_get_sample_at_pos 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos 
  );
  /* If {bps} zero, ignores {el} and returns 0. Otherwise extracts the
    value of the element with position {pos} from the area {el}.
    No address range checking is done on {pos}. */

void ppv_set_sample_at_pos 
  ( void *el, 
    ppv_nbits_t bps, 
    ppv_nbits_t bpw, 
    ppv_pos_t pos, 
    ppv_sample_t smp 
  );
  /* If {bps} is zero, ignores {el} and does nothing. Otherwise stores
    value {smp} into the element with position {pos} in the area {el}.
    
    No address range checking is done on {pos}. Fails if {smp} is not in
    {0..2^bps-1}, but does not test it against the array's {maxsmp}. */

void ppv_sample_range(ppv_array_t *A, ppv_sample_t *vminP, ppv_sample_t *vmaxP);
  /* Finds the maximum and mininmum values of all the samples in {A},
    and returns them in {*vminP,*vmaxP}.  If {A} is empty, sets
    {*vminP} to {2^A.bps-1} and {*vmaxP} to zero. */

/* INDEX TUPLE MANIPULATION */

void ppv_index_clear ( ppv_dim_t d, ppv_index_t ix[] );
  /* Sets {ix[ax] = 0} for every axis {ax}. */

void ppv_index_assign ( ppv_dim_t d, ppv_index_t ix[], const ppv_index_t val[] );
  /* Sets {ix[ax] = val[ax]} for every axis {ax}. */

void ppv_index_shift ( ppv_dim_t d, ppv_index_t ix[], ppv_index_t *inc );
  /* Sets {ix[ax] += inc[ax]} for every axis {ax}. */

sign_t ppv_index_compare ( ppv_dim_t d, const ppv_index_t ixa[], const ppv_index_t ixb[] );
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
      
      (3) The field {A.maxsmp} is at most {2^A.bps-1}.

      (4) If the array is memoryless, then {A.base} and all {A.step}
        coefficients are zero, and {A.el} is NULL.

      (5) If the array is not memoryless, the storage area 
        of every element is contained in the allocated 
        area pointed by {A.el}.
        
    The procedure checks conditons (1)--(5). It returns TRUE if the
    attributes of {A} are valid. Otherwise, it either aborts the program with
    an error message (if {die=TRUE}) or returns FALSE (if {die=FALSE}).
  */

/* DETAILS OF ELEMENT PACKING */

#define ppv_packing_INFO \
  "General principles\n" \
  "\n" \
  "  In order to save storage and to simplify image transfer between\n" \
  "  devices and computers with different architectures, the storage\n" \
  "  area {A.el} is organized as a vector of /words/. Each\n" \
  "  word is an unsigned binary integer with {A.bpw} bits. The\n" \
  "  parameter {A.bpw} may be  different from the machine's native word size.\n" \
  "  For this implementation, the word size {A.bpw} must be 8, 16, or\n" \
  "  32. The implementation assumes that words can be directly\n" \
  "  addressed by a pointer of the appropriate type.\n" \
  "\n" \
  "  Each element of an array {A} is stored in {A.bps} bits. Depending\n" \
  "  on {A.bps} and {A.bpw}, one element may span two or more words, or\n" \
  "  one word may contain two or more elements. The sample storage area\n" \
  "  address {A.el} is the address of the word\n" \
  "  that contains the highest-order bit of the sample with position 0.\n" \
  "\n" \
  "Packing of \"small\" samples\n" \
  "\n" \
  "  If {bps} is less than or equal to {bpw}, each element is entirely\n" \
  "  contained in a single word, and {K=A.bpw/A.bps} elements are\n" \
  "  allocated into each word. Within each word, samples with consecutive\n" \
  "  positions are tightly packed together. Specifically a sample whose\n" \
  "  position is {pos} will be stored in word {iw = pos/K}, and its units\n" \
  "  bit will be bit {ib = (K-1 - pos%K)*A.bps} of that word If {A.bps}\n" \
  "  is not an exact divisor of {A.bpw}, any incompletely filled words\n" \
  "  will be padded with zeros at the high-order end.\n" \
  "\n" \
  "  For example, consider an image with {A.bps=6} and 5 samples\n" \
  "  {A,B,C,D,E,F,G} with positions {pos = 0,1,2,3,4,5,6}. With {bpw=8}, we\n" \
  "  would have {K = 8/6 = 1} samples per word, and the samples would\n" \
  "  be stored in memory as follows:\n" \
  "\n" \
  "    | word#   <--0---> <--1---> <--2---> <--3---> <--4---> <--5---> <--6--->\n" \
  "    | sample#   <-0-->   <-1-->   <-2-->   <-3-->   <-4-->   <-5-->   <-6-->\n" \
  "    |         ==AAAAAA ==BBBBBB ==CCCCCC ==DDDDDD ==EEEEEE ==FFFFFF ==GGGGGG\n" \
  "    | bit#    76     0 76     0 76     0 76     0 76     0 76     0 76     0\n" \
  "\n" \
  "  With {bpw = 16}, we would have {K = 16/6 = 2}, and the following layout:\n" \
  "\n" \
  "    | word#   <------0-------> <-----1--------> <-----2--------> <-----3-------->\n" \
  "    | sample#     <-0--><-1-->     <-2--><-3-->     <-4--><-5-->     <--->       \n" \
  "    |         ====AAAAAABBBBBB ====CCCCCCDDDDDD ====EEEEEEFFFFFF ====GGGGGG======\n" \
  "    | bit#    1  1     0     0 1  1     0     0 1  1     0     0 1  1     0     0\n" \
  "    |         5  2     6     0 5  2     6     0 5  2     6     0 5  2     6     0\n" \
  "\n" \
  "  With {bpw = 32}, we would have {K = 32/6 = 5}, and the following layout:\n" \
  "\n" \
  "    | word#   <---------------------1--------> <---------------------1-------->\n" \
  "    | sample#   <-0--><-1--><-2--><-3--><-4-->   <-5--><-6-->                  \n" \
  "    |         ==AAAAAABBBBBBCCCCCCDDDDDDEEEEEE ==FFFFFFGGGGGG==================\n" \
  "    | bit#    33     2     1     1     0     0 33     2     1     1     0     0\n" \
  "    |         10     4     8     2     6     0 10     4     8     2     6     0\n" \
  "\n" \
  "Packing of \"big\" samples\n" \
  "\n" \
  "  If {A.bps} is greater than {A.bpw}, then each sample will occupy\n" \
  "  {K=(bps+bpw-1)/bpw} consecutive words, with word indices {(pos-1)*K}\n" \
  "  to {pos*K-1}.  The units bit of the latter will be the units bit\n" \
  "  of the sample.  If {A.bps} is not an exact multiple of {A.bpw},\n" \
  "  the first word, with index {(pos-1)*K}, will be padded with zero\n" \
  "  bits at the high-order end.\n" \
  "\n" \
  "  For example, if {bps = 18}, the layout for {bpw = 8} will be\n" \
  "\n" \
  "    | word#   <--0---> <--1---> <--2---> <--3---> <--4---> <--5---> ...\n" \
  "    | sample#       <-------0---------->       <--------1---------> ...\n" \
  "    |         ======AA AAAAAAAA AAAAAAAA ======BB BBBBBBBB BBBBBBBB ...\n" \
  "    | bit#    7    2 0 7      0 7      0 7    2 0 7      0 7      0 ...\n" \
  "\n" \
  "  while that for {bpw = 16} will be\n" \
  "\n" \
  "    | word#   <------0-------> <-----1--------> <-----2--------> <-----3--------> ...\n" \
  "    | sample#               <---------0------->               <--------1--------> ...\n" \
  "    |         ==============AA AAAAAAAAAAAAAAAA ==============BB BBBBBBBBBBBBBBBB ...\n" \
  "    | bit#    1            0 0 1              0 1            0 0 1              0 ...\n" \
  "    |         5            2 0 5              0 5            2 0 5              0 ...\n" \
  "\n" \
  "  and that for {bpw = 32} will be\n" \
  "\n" \
  "    | word#   <--------------0---------------> <--------------2---------------> ...\n" \
  "    | sample#               <-------0-------->               <-------1--------> ...\n" \
  "    |         ==============AAAAAAAAAAAAAAAAAA ==============BBBBBBBBBBBBBBBBBB ...\n" \
  "    | bit#    3             1                0 3             1                0 ...\n" \
  "    |         1             7                0 1             7                0 ...\n" \
  ""

size_t ppv_tot_sample_bytes ( ppv_sample_count_t npos, ppv_nbits_t bps, ppv_nbits_t bpw );
  /* Returns the minimum number of bytes needed to store all the voxels
    of an array a total of {npos} elements, given that each sample has
    {bps} bits and the samples are packed in words of {bpw} bits each.
    
    Assumes that the samples are packed as closely as allowed by the
    packing rules above, with no element replication.  In particular,
    returns 0 if {npos == 0} or {bps == 0}. */

/* TESTING AND DEBUGGING */

void ppv_choose_test_size(ppv_dim_t d, ppv_sample_count_t npos, ppv_size_t sz[]);
  /* Chooses varied array sizes {sz[0..d-1]} so that the total number of samples
    will be about {npos}.  
    
    The min size will be at least 1. If {d >= 2}, the max size will be
    about 3 times the min size. */
    
void ppv_dump_storage(FILE *wr, void *el, ppv_pos_t pos, ppv_nbits_t bps, ppv_nbits_t bpw, int nxw);
  /* Prints the word(s) of the storage area {*el} that contain(s) the 
    element with position {pos}, plus {nxw} words of context.
    Assumes that the element is {bps}  bits long and that elements
    are packed into words of {bpw} bits each, as explained in {ppv_packing_INFO}.
    Also prints the decimal value of the sample stored there. */

/* RANDOM FILLING 

  The procedures in this section use the the pseudorandom generators
  from {jsrandom.h}. To ensure that the results are repeatable or
  non-repeatable, as needed, call {srandom} from {stdlib.h} first, with
  an appropriate {seed}. */

void ppv_throw_noise(ppv_array_t *A);
  /* Fills the array {A} with pseudo-random (very pseudo) sample values
    in the range {0..A.maxsmp}. */

void ppv_throw_balls(ppv_array_t *A);
  /* Seta all elements of {A} to sample value 0 except for a few balls
    of random centers and radii, possibly overlapping, whith 
    sample values in {1..A.maxsmp}. */

#endif
