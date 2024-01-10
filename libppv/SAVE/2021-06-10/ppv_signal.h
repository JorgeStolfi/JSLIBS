/* Sampled multi-dimensional signals. */
/* Last edited on 2021-06-10 13:36:46 by jstolfi */

#ifndef ppv_signal_H
#define ppv_signal_H

#define _GNU_SOURCE
#include <ppv_types.h>
#include <ppv_array.h>
#include <interval.h>
#include <stdio.h>

/* 
  Multi-dimensional sampled signals
  
  For this interface, a /signal/ a real function {f()} defined on some
  multi-dimensional domain, represented by an array of /samples/ which
  are interpreted, interpolated and extrapolated according to specific
  methods.

  A PPV signal_t can be used to represent an audio stream, a plane
  or space curve, a two-dimensional image or repeating pattern, a
  three-dimensional tomogram, a video clip, etc.. */

#define ppv_signal_NAXES ppv_array_NAXES
  /* Max dimension of the domain of a signal. */

typedef 
  struct ppv_signal_t {
    ppv_array_t s;                         /* Sample array. */
    double *vtb                            /* Sample decoding table. */
    interval_t dom[ppv_signal_NAXES];      /* Fundamental coordinate interval along each axis. */
    interval_t rng;                        /* Value range of signal. */
    ppv_extrap_t extrap[ppv_signal_NAXES]; /* Extrapolation method along each axis. */
    ppv_interp_t interp[ppv_signal_NAXES]; /* Interpolation method along each axis. */
  } ppv_signal_t; 
  /* A {ppv_signal_t} is descriptor for a signal. */

/* 
  Domain and range

    The domain of any {ppv_signal_t} is a subset of {R^N} where {N =
    ppv_signal_NAXES}. If {f} is such a signal, and {x[0..N-1]} is a vector
    of real numbers, we denote by {f(x)} the value of signal {f}
    at point {x}.

    A signal whose domain has dimension {k<N} can be represented by
    an {N}-dimensional signal which is constant (and sampled only once)
    along the excess {N-k} axes.

    This interface assumes that a signal has a single scalar value
    for each argument. However, vector-valued (multi-channel)
    signals can be accomodated by taking the channel index as an
    additional domain coordinate.
    
  Fundamental cell

    Along each axis {i}, a signal {f} has a /fundamental interval/, a
    real interval {f.dom[i]}. The cartesian product of the fundamental
    intervals is a box {f.dom}, the /fundamental cell/ of {f}.

    Roughly speaking, the signal is sampled and stored only for points
    {x} of {R^N} that lie within the fundamental cell. The signal values 
    inside the cell may be extended to the whole {R^N} by extrapolation,
    as described further on. */
    
double ppv_signal_eval(ppv_signal_t *f, double x[]);
  /* Returns the value of {f(x)}, for a given point {x[0..N-1]} of {R^N}.
    The result may be {±oo} or NaN. */

/*
  Samples and sampling positions

    The value of {f(x)} is determined by an array {f.s} of /samples/.

    Each sample has a /nominal position/, which is a vertex of a regular
    rectangular grid that spans with the fundamental cell.  

    More precisely, for each axis {i} consider the fundamental
    interval {f.dom[i]} divided into a list of {f.s.size[i]} equal
    sub-intervals, and let {ix[0..N-1]} be the indices of some sample
    in {f.s}. The nominal {i}-coordinate of that sample is the
    midpoint of sub-interval number {ix[i]} in this list. */

ppv_signal_t ppv_signal_new(ppv_array_t A);
  /* Creates a {ppv_signal_t} record from a sample array. For each
    axis {i}, the fundamental interval {f.dom[i]} is set to {[0 _
    A.size[i]]}; the interpolation degree {f.interp[i]} is set to 0
    (no interpolation); and the extrapolation methods {f.extrap[i]} is set
    to {ppv_extrap_NONE}. The value table {f.tb} is set to NULL, and
    the value range is set to {[-1 _ +1]}. */

/*     
  Sample value decoding
  
    Each sample is a real number, quantized as a small unsigned integer
    and stored in memory in a packed binary format. Each sample is
    some average of the signal value in the neighborhood of its nominal 
    position. 
    
    Sample values must be /decoded/ (converted to real values) before
    being used. First, the sample is interpreted as an unsigned
    integer {q} in the range {0..2^B-1} where {B = f.s.bps}. The
    integer {q} is then converted to a real value {r} in the range
    {[0_1]}, through the table {f.vtb}. Namely, if {f.vtb} is not
    NULL, it must be a vector with {2^B} elements; and then {r} is
    {f.vtb[q]}, which should be in the range {[0_1]}. If {f.vtb} is
    NULL, then {r} is {(q + 0.5)/2^B}. In any case, {r} is then
    converted to the final sample value by the affine map that takes
    the interval {[0_1]} to the range {f.rng}. */
    
double_vec_t *ppv_signal_table_mu_law(int n, double lo, double hi);
  /* Returns a table that maps integers {0..n-1} to 
    real values spaced according to the mu-law encoding formula,
    shifted and scaled to the range {[lo _ hi]}.  
    
    To use the result as the field {f.vtb} of a {ppv_signal_t},
    one should specify {lo = 0}, {hi = 1}, and {f.rng = [-1 _ +1]}.
    One could also use {lo = -1}, {hi = +1}, and {f.rng = [0 _ 1]}. */
    
double_vec_t *ppv_signal_table_gamma(int n, double lo, double hi, double gamma);    
  /* Returns a table that maps each integer {q} in {0..n-1} to 
    the real value {(q/n)^gamma}, scaled and shifted to
    the range {[lo _ hi]}.
    
    To use the result as the field {f.vtb} of a {ppv_signal_t},
    one should specify {lo = 0}, {hi = 1}, and {f.rng = [-1 _ +1]}.
    One could also use {lo = -1}, {hi = +1}, and {f.rng = [0 _ 1]}. */
    
/* 
  Interpolation
    
    The value of {f(x)}, for any {x} in the fundamental cell, is
    computed by interpolating the decoded values of samples whose
    nominal positions lie near {x}.
    
    The interpolation scheme along the domain axis {i} is determined by the
    attribute {f.interp[i]}, the /interpolation degree/, which is a small unsigned integer.  
    This integer defines the /window size/ {NW = f.interp[i]+1}.

    Let's define an /{i}-slice/ of a {K}-dimensional array as a
    sub-array which has the same index ranges, except that the size
    along axis {i} has been reduced to 1. The interpolation is
    recursive on the domain dimension. Namely, to interpolate a
    {K}-dimensional array, we interpolate {NW} {i}-slices of it, where
    {i = K-1}; and then perform a uni-dimensional interpolation of
    degree {f.interp[i]} on the {NW} resulting values.
    
    The {NW} slices are chosen so that their nominal {i}-coordinates
    are as centered on {x[i]} as possible. Namely, if {NW} is odd, the
    middle slice is the one whose {i}-coordinate is closest to {x[i]}.
    If {NW} is even, the two midle slices are those that bracket
    {x[i]}.
    
    In either case, the decoded sample values are interpolated with a
    B-spline reconstruction filter of degree {f.interp[i]}, yielding a
    signal with continuity order {CT = f.interp[i]-1}.
    
    In particular, if {f.interp[i]} is 0, there is no interpolation,
    and the result depends only on the {i}-slice whose nominal
    {i}-coordinate is nearest to {x[i]}. If {f.interp[i]} is 1, the
    signal value is obtained by linear interpolation of two {i}-slices
    that bracket {x[i]}. */
    
double ppv_signal_univariate_interp(int deg, double s[], double z);
  /* Returns the univariate interpolation of sample values {s[0..deg]}
     at the argument {z}; which should lie in {[0_1]} 
     if {deg} is odd, and in {-0.5 _ +0.5]} if {deg} is even. */
    
/*    
  Extrapolation
    
    For points {x} outside the fundamental cell, the value of {f(x)}
    is either undefined, or is defined by interpolation of an infinite sample
    array that is obtained by extending the array {f.s}.  Extrapolation 
    along axis {i} is determined by the parameter {f.extrap[i]}, which may be:
      
    Domain of definition
    
      If the extrapolation method along axis {i} is NONE, then the
      signal is defined only at values of {x[i]} for which all samples
      needed by the interpolation method exist in the array {f.s}. The
      set of such {x[i]} is a sub-inderval of {f.dom[i]}, called the
      /definition interval/ {f.def[i]} along axis {i}. In all other
      cases, {f(x)} is defined for any value of {x[i]}; we then define
      {f.def[i] == f.dom[i]} by convention.
    
    Nominal sample positions
    
      The nominal sample positions along axis {i} depend on the
      interpolation order {f.interp[i]} and the extrapolation order
      {f.extrap[i]}.
      
      Consider a sample whose index tuple is {ix[0..N-1]}. The 
      nominal {i}-coordinate of that sample is {LO +
      (ix[i]+0.5)*(HI-LO)/SZ}, where {LO = f.dom[i].end[0]}, {HI =
      f.dom[i].end[1]}, and {SZ = f.s.size[i]}.  */

void ppv_signal_set_interpolation(ppv_signal_t *f, ppv_axis_t i, ppv_interp_t interp);
  /* Sets the interpolation method of {f} along axis {i} to {interp}.
    If {f.extrap[i]} is NONE, this procedure will also adjust {f.dom[i]}
    so that the nominal positions of the samples are not changed. */

void ppv_signal_set_extrapolation(ppv_signal_t *f, ppv_axis_t i, ppv_extrap_t extrap);
  /* Sets the extrapolation method of {f} along axis {i} to {extrap}.
    If {f.extrap[i]} changes to or from NONE, this procedure will also adjust {f.dom[i]}
    so that the nominal positions of the samples are not changed. */

/* EVALUATION */

ppv_pos_t ppv_signal_eval(ppv_signal_t *f, double []);
  /* Evaluates the signal {f} at point {x}, using its current
    interpolation and extrapolation methods. */

/* SAMPLE EXTRACTION AND INSERTION */

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

bool_t ppv_index_first ( ppv_index_t *ix, ppv_signal_t *A );
  /* Sets {ix[i] = 0} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

bool_t ppv_index_last ( ppv_index_t *ix, ppv_signal_t *A );
  /* Sets {ix[i] = A.size[i]-1} for every axis {i}, and returns TRUE iff that
    is a valid index tuple for {A} (i.e. {A.size[i]>0} for all {i}). */

void ppv_index_shift ( ppv_index_t *ix, ppv_index_t *inc );
  /* Sets {ix[i] += inc[i]} for every axis {i}. */

sign_t ppv_index_compare ( ppv_index_t *ixa, ppv_index_t *ixb );
  /* Returns {NEG}, {ZER}, or {POS} depending on whether {ixa} is less
    than, equal to, or greater than {ixb} in reverse lexicographic
    order (where the first index varies the fastest). */

bool_t ppv_index_next ( ppv_index_t *ix, ppv_signal_t *A, ppv_dim_t d, ppv_pos_t *p );
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

bool_t ppv_index_prev ( ppv_index_t *ix, ppv_signal_t *A, ppv_dim_t d, ppv_pos_t *p );
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
    
bool_t ppv_index_is_valid ( ppv_index_t *ix, ppv_signal_t *A );
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

void ppv_crop ( ppv_signal_t *A, ppv_axis_t i, ppv_size_t skip, ppv_size_t  keep);
  /* Crops the array {A} to the range {skip..skip+keep-1} along axis
    {i}. Thus {size[i]} becomes {keep}, and using value {r} for that
    index in the new descriptor is equivalent to using {r + skip} in
    the old one.
    
    The value of {skip+keep} must not exceed the original {size[i]}.
    {keep} may be zero, in which case the array becomes empty (but
    retains its size along the other axes) and all its steps are
    reset to zero. */

void ppv_subsample ( ppv_signal_t *A, ppv_axis_t i, ppv_size_t stride );
  /* Subsamples the descriptor {A} along axis {i}, so that using {r}
    for that index in the new descriptor is equivalent to using
    {r*stride} in the old one.
    
    The increment {stride} must be positive. The extent {size[i]} is
    reduced to the maximum {sz} such that {(sz-1)*stride} is less than 
    the original {size[i]}. In particular, {size[i]} will be 0 
    if and only if it was 0 originally. */

void ppv_flip ( ppv_signal_t *A, ppv_axis_t i );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that the index {i}
    runs in the opposite direction.  */

void ppv_replicate ( ppv_signal_t *A, ppv_axis_t i, ppv_size_t sz );
  /* Modifies the descriptor {A} so that it has size {sz} along axis
    {i}, but using any value for that index yields the same sample as
    using the value 0. Namely, this procedure sets {step[i]} to 0,
    {size[i]=sz}. On input, {size[i]} must be 1, and {sz} must be
    positive. */

void ppv_swap_indices ( ppv_signal_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the same indices, except that indices {i} and
    {j} are exchanged. Thus, for instance, {ppv_swap_indices(img, 1,2)} 
    could be used to exchange the rows and columns of a 2D color image,
    assuming index 0 is the color channel.
    If {i == j}, the procedure has no effect. */

void ppv_flip_indices ( ppv_signal_t *A, ppv_axis_t i, ppv_axis_t j );
  /* Modifies the descriptor {A} so that it describes the same samples
    as before, with the indices {i..j} in reverse order. Thus,
    after the call the call {ppv_flip_indices(A,2,5)}, 
    element {A[p,q,r,s,t,u]} is equivalent to {A[p,q,u,t,s,r]}
    of the original array. */

void ppv_diagonal ( ppv_signal_t *A, ppv_axis_t i, ppv_axis_t j );
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

void ppv_chop ( ppv_signal_t *A, ppv_axis_t i, ppv_size_t sz, ppv_axis_t j );
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

void ppv_enum ( ppv_op_t op, ppv_signal_t *A, ppv_signal_t *B, ppv_signal_t *C );
  /* Enumerates all valid samples of the arrays {A,B,C} in parallel,
    with matching indices. For each index tuple {[ix[0],..ix[N-1]]}, calls
    {op(ix, pA, pB, pC)} where {pA,pB,pC} are the
    positions of the corresponding samples in the three arrays. 
    
    If an array is NULL, it is not scanned and the corresponding {pos}
    is always 0. All {ppv_signal_t} arguments that are not NULL must
    have the same {size} vector. */

/* DESCRIPTOR PRINTOUT */

void ppv_print_descriptor ( FILE *wr, char *pf, ppv_signal_t *A, char *sf );
  /* Prints the attributes of descriptor {A} to file {wr},
    bracketed by the strings "pf" and "sf". */
  
/* DESCRIPTOR VALIDATION */

bool_t ppv_descriptor_is_valid ( ppv_signal_t *A, bool_t die );
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
