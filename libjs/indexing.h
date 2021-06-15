#ifndef indexing_H
#define indexing_H

/* Multidimensional array indexing tools */
/* Last edited on 2021-06-14 08:47:12 by jstolfi */

#include <bool.h>
#include <sign.h>
#include <stdint.h>

/* !!! Eliminate {reverse} since one can use {ix_flip}. !!! */

/* These procedures provide tools for simple multi-dimensional array
  indexing. They assume arrays with box-like shape, indices that
  start at 0, and array elements stored at positions that depend
  linearly on the indices.
  
  INDEXING FORMULA
  
  More precisely, let {M} be an array with {d} indices {ix[0..d-1]},
  where each {ix[i]} ranges from 0 to some maximum {sz[i]-1}. This
  module assumes that element {M[ix[0],..ix[d-1]} will be stored in
  element {Mv[p]} of some vector {Mv}; the element's /position/
  {p} being given by
  
    {  p = bp + SUM { ix[i]*st[i] : i = 0,..d-1 } }         (1)
    
  for some /base position/ {bp}, and some /steps/ (position
  increments) {st[i]}.
  
  STEP MANIPULATION
  
  The advantage of having an explicit step vector {st[]}, instead of
  computing the steps from the sizes {sz[]}, is that one can perform
  many array indexing transformations in {O(1)} time, without actually
  moving or copying the elements.
  
  In particular, by manipulating the parameters {st}, {sz} and {bp}
  one can extract arbitrary sub-arrays (rows, columns, planes, blocks,
  etc) of any array; reverse the order of vectors; transpose, rotate,
  and flip arrays along any axis; subsample the elements along
  any axis; extract diagonal vectors and planes; and more.
  
  VIRTUAL REPLICATION
  
  Any step {st[i]} can be zero; in that case, the corresponding index
  {ix[i]} is irrelevant. This feature can be used to virtually fill an
  array with copies of a single element, or of a lower-dimensional
  slice, without actually duplicating the samples.
  
  TRIVIAL INDICES
  
  When the array {M} has size {sz[i] = 1} along some axis {i}, the only
  valid value for index {ix[i]} is zero. We then say that {M} is
  /trivial along/ that axis, os that {i} is a /trivial axis/ for {M}.
  
  EMPTY ARRAYS
  
  If {sz[i]==0} for at some axis {i}, the array has zero elements ---
  it is an `empty array'. Such arrays still may have non-zero size
  along other axes, and that information is significant in many
  contexts. With the noted exceptions, all procedures below can be
  used on empty arrays, usually with sensible results.
  
  ELEMENT STORAGE AND MANIPULATION (NOT!)
  
  This module is concerned only with computing the mapping from
  indices {ix[0..d-1]} to the position {p}, and vice-versa; and with
  index manipulations like array scanning and such. Storage allocation
  for the elements and manipulation of their values are left to the
  client. In particular, the `position' may be a memory (byte)
  address, a a bit address, an index into some vector, or something
  else entirely. */

/* DATA TYPES */

typedef uint8_t ix_dim_t;
  /* Type for the number of indices {d}. */

typedef uint8_t ix_axis_t;
  /* Type for axis number, or index of an index, usually from 0 to {d-1}. */

typedef uint64_t ix_size_t;
  /* Type for the size {sz[i]} of an array along some axis. */

typedef int64_t ix_index_t;
  /* Type for individual array element indices {ix[i]}.
    They may be negative on occasion (e.g. in reverse 'for' loops).  */

typedef int64_t ix_step_t; 
  /* Type for position steps {st[i]} along any axis. It is signed
    to allow for array flipping. */

typedef uint64_t ix_pos_t; 
  /* Type for element positions. */

#define ix_pos_NONE (UINT64_MAX)
  /* An {ix_pos_t} value that means "no such element". */

typedef uint64_t ix_count_t;
  /* A count of the number of elements in an array. */

/* ELEMENT INDEXING */

ix_pos_t ix_position 
  ( ix_dim_t d, 
    const ix_index_t ix[], 
    ix_pos_t bp, 
    const ix_step_t st[]
  );
  /* Computes the position of element {M[ix[0],..ix[d-1]}, that is,
    {bp + SUM{ix[i]*st[i] : i=0,..d-1}}. Does not test the validity of
    the indices, so the result is undefined if any {ix[i]} is
    negative, or exceeds its maximum legal value. */

ix_pos_t ix_position_safe 
  ( ix_dim_t d, 
    const ix_index_t ix[], 
    const ix_size_t sz[], 
    ix_pos_t bp, 
    const ix_step_t st[]
  );
  /* Checks whether the index tuple {ix[0..d-1]} is valid, using
    {ix_is_valid(d,ix,sz)}. If it is valid, returns the position of
    element {M[ix[0],..ix[d-1]}, computed with
    {ix_position(d,ix,bp,st)}. If it is invalid, returns
    {ix_pos_NONE}. */

void ix_indices 
  ( ix_pos_t p, 
    ix_dim_t d, 
    const ix_size_t sz[], 
    const ix_step_t st[], 
    ix_pos_t bp, 
    ix_index_t ix[] 
  );
  /* Computes the indices {ix[0,..d-1]} of an element from its given
    position {p}, the base position {bp}, the sizes {sz[0..d-1]},
    and the position steps {st[0..d-1]}. Undefined if {p} is not 
    a valid position computable from {sz}.
    
    Currently it assumes that there is an ordering of the steps such
    that the slices in that order have disjoint position ranges. That
    is true, in particular, if the elements are obtained by clipping,
    subsampling, and index swapping from a packed element arrangement
    (in C or Fortran order). */

/* INDEX TUPLE MANIPULATION */

void ix_fill ( ix_dim_t d, ix_index_t ix[], ix_index_t val );
  /* Sets {ix[i] = val} for {i=0,..d-1}. */

void ix_assign ( ix_dim_t d, ix_index_t ix[], const ix_index_t val[] );
  /* Sets {ix[i] = val[i]} for {i=0,..d-1}. */

bool_t ix_assign_min ( ix_dim_t d, ix_index_t ix[], const ix_size_t sz[] );
  /* Sets {ix[i] = 0} for {i=0,..d-1}, and returns TRUE iff that
    is a valid index tuple (i.e. {sz[i]>0} for all {i}). */

bool_t ix_assign_max ( ix_dim_t d, ix_index_t ix[], const ix_size_t sz[] );
  /* Sets {ix[i] = sz[i]-1} for {i=0,..d-1}, and returns TRUE iff that
    is a valid index tuple (i.e. {sz[i]>0} for all {i}). */

void ix_shift ( ix_dim_t d, ix_index_t ix[], const ix_index_t inc[] );
  /* Sets {ix[i] += inc[i]} for {i=0,..d-1}. */

bool_t ix_is_valid ( ix_dim_t d, const ix_index_t ix[], const ix_size_t sz[] );
  /* Returns TRUE iff {ix[i]} lies in the range {0..sz[i]-1}, 
    for {i=0,..d-1}. */

/* SIZE TUPLE MANIPULATION */

void ix_sizes_assign ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = szb[i]} for {i=0,..d-1}. */

bool_t ix_sizes_shrink ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = min(sza[i], szb[i])} for {i=0,..d-1}.
    Returns TRUE iff any {sza[i]} was actually changed. */

bool_t ix_sizes_expand ( ix_dim_t d, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = max(sza[i], szb[i])} for {i=0,..d-1}.
    Returns TRUE iff any {sza[i]} was actually changed. */

ix_size_t ix_max_size ( ix_dim_t d, const ix_size_t sz[] );
  /* The maximum of {{0, sz[0],.. sz[d-1]}}. */

ix_size_t ix_min_size ( ix_dim_t d, const ix_size_t sz[] );
  /* The minimum of{{ix_MAX_SIZE, sz[0],.. sz[d-1]}}. */

/* STEP TUPLE MANIPULATION */

void ix_steps_assign ( ix_dim_t d, ix_step_t sta[], const ix_step_t stb[] );
  /* Sets {sta[i] = stb[i]} for {i=0,..d-1}. */

/* ELEMENT COUNTING */

bool_t ix_is_empty ( ix_dim_t d, const ix_size_t sz[] );
  /* TRUE iff the array is empty, i.e. {sz[i]==0} for some {i} in {0..d-1}. */

ix_count_t ix_num_tuples ( ix_dim_t d, const ix_size_t sz[] );
  /* Number of valid index tuples in array; 0 if the array is empty.
    Beware that virtual element replication may create huge arrays that exceed
    the capacity of an {ix_count_t}. */

ix_count_t ix_num_positions ( ix_dim_t d, const ix_size_t sz[], const ix_step_t st[] );
  /* Number of *distinct* positions, not counting virtual replications.
    Assumes that the parameters are valid (see {ix_parms_are_valid} below). */

ix_pos_t ix_min_pos ( ix_dim_t d, const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[] );
ix_pos_t ix_max_pos ( ix_dim_t d, const ix_size_t sz[], ix_pos_t bp, const ix_step_t st[] );
  /* Minimum and maximum position for any valid index tuple.
    Returns {ix_pos_NONE} if the array is empty. */

bool_t ix_same_size ( ix_dim_t d, const ix_size_t sza[], const ix_size_t szb[], bool_t die );
  /* The procedure returns TRUE if {sza[i] == szb[i]} for all {i} in
    {0..d-1}.  Otherwise, it either aborts the program with an error
    message (if {die=TRUE}) or returns FALSE (if {die=FALSE}). */ 

bool_t ix_contained ( ix_dim_t d, const ix_size_t sza[], const ix_size_t szb[], bool_t die );
  /* The procedure returns TRUE if {sza[i] <= szb[i]} for all {i} in
    {0..d-1}. Otherwise, it either aborts the program with an error
    message (if {die=TRUE}) or returns FALSE (if {die=FALSE}). */ 

/* STEP MANIPULATION */

/* The following procedures modify the indexing parameters {sz}, {*bp}
  and {st} so as to change the set of valid indices and/or the mapping
  between indices and positions, as described. Obviously, no array elements are
  actually allocated, reclaimed, or changed.
  
  In the describing comments, {ixA[0..d-1]} is a generic tuple valid
  AFTER the call, and {ixB[0..d-1]} is the tuple that, BEFORE the call,
  would yield the same position.
  
  All the procedures below preserve the validity of the parameters
  {sz,*bp,st} (as defined by {ix_parms_are_valid} below).  Unless said otherwise,
  the operations preserve the empty/non-empty status of the array. */

void ix_crop 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t skip, 
    ix_size_t keep
  );
  /* Crops the array along axis {i} to the index range {ini..fin}. 
  
    More precisely, sets {sz[i]} to {keep}, and modifies {*bp,st}
    so that {ixB[i] = ixA[i]+skip} and {ixB[k] = ixA[k]} for all {k\neq i}.
    If {keep==0}, the array becomes empty (but retains its size along
    the other axes); in that case {*bp} and {st[0..d-1]} are reset to 0.
    The result is undefined if {skip+keep > sz[i]}.
    
    Thus, to clip a two-dimensional {m × n} matrix to its upper right
    quarter, use {ix_crop(d,sz,&bp,st, 0, 0,m/2)} and then
    {ix_crop(d,sz,&bp,st, 1, n/2,n-n/2)}, where {m=sz[0]} and {n=sz[1]}. 
    To reduce that matrix to its column number {c}, use
    {ix_crop(d,sz,&bp,st, 1, c,1)}. */

void ix_subsample 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t stride 
  );
  /* Subsamples the array along axis {i} by taking one every {stride} slices.
    The sampling step {stride} must be positive. 
    
    More precisely, reduces {sz[i]} to {floor((sz[i]-1)/stride)+1}, and modifies
    {*bp} so that {ixB[i] = stride*ixA[i]} and {ixB[k] = ixA[k]} for all {k\neq i}. */

void ix_flip 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i
  );
  /* Flips the array along axis {i}. 
  
    More precisely, modifies {*bp} and {st[i]} so that {ixB[i]=sz[i]-1-ixA[i]}, 
    and {ixB[k]=ixA[k]} for all {k\neq i}. The size {sz[i]} is not changed.  */

void ix_replicate 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t size 
  );
  /* Virtually replicates the array along axis {i}, {size} times.
  
    On input, the array must be trivial along axis {i} ({sz[i] == 1}).
    The procedure sets {sz[i]=size} and {st[i]=0}, so that {ixB[i]=0} and
    {ixB[k]=ixA[k]} for all {k\neq i}. The replication count {size} must be
    positive. */

void ix_swap_indices 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_axis_t j,
    ix_dim_t n
  );
  /* Transposes the array, by exchanging the {n} indices starting with index {i}
    with the {n} indices {j} starting with index {j}. The two sets of indices
    must be either identical ({i==j} or {n==0}, in which case the procedure has
    no effect) or disjoint.
  
    More precisely, swaps {sz[i+k]} with {sz[j+k]} and {st[i+k]} with {st[j+k]},
    for {k} in {0..n-1}; so that {ixB[i+k]=ixA[j+k]}, {ixB[j+k] = ixA[i+k]}, for
    {k} in {0..n-1}, and {ixB[r]=ixA[r]} for every other axis {r}.
    
    Thus, to transpose a two-dimensional matrix, in the ordinary sense,
    use {ix_swap_indices(d,sz,&bp,st, 0,1, 1)}.  */

void ix_flip_indices
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_axis_t j 
  );
  /* Reverses the order of all array indices between axis {i} and axis {j}, 
    inclusive.
  
    More precisely: if {i<j}, the procedure swaps {sz[i+k]} with
    {sz[j-k]} and {st[i+k]} with {st[j-k]} for all {k} in
    {0..(j-i-1)/2}; thus {ixB[i+k] = ixA[j-k]} for all those {k}. 
    It is a no-op if {i>=j}. */
 
void ix_slice
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_dim_t n, 
    const ix_axis_t ax[], 
    const ix_index_t ix[] 
  );
  /* Sets the {n} indices along axes {ax[0..n-1]} to the values
    {ix[0..n-1]}, and then performs an index permutation that moves those axes
    to the end of the index tuple and brings the rest of the indices to the
    front, preserving their order.
    
    The axes {ax[0..n-1]} must be a subset of {0..d-1} and must be given in
    strictly increasing order. If {ax} is NULL, assumes {ax[k] = k} for any {k}
    in {0..n-1}, i.e. sets the first {n} indices and moves them to the end. The
    procedure is a no-op if {n==0}; otherwise, each index value {ix[k]} must be
    valid for the axis {ax[k]}, so the array must not be empty.
    
    Thus, for example, to extract a two-dimensional slice 
    {B} from a six-dimensional array {A}, such that {B[r,u] = A[15,r,10,20,u,17]}
    for all {r,u} valid in the right-hand side, use {ix_slice(d,sz,&bp,st,4,ax,ix)}
    where {ax[0..3]=(0,2,3,5)} and {ix[0..3] = (15,10,20,17)}. */
   
void ix_diagonal 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i,
    ix_axis_t j 
  );
  /* Extracts a diagonal slice of the array, cut parallel to the
    bisector of axes {i} and {j}; and redefines axis {i} to be along
    that diagonal. Thus, using index values {r} and {s} for those axes,
    after the call, is equivalent to using {r} and {r+s} before the call.
  
    Upon entry, requires {i != j} and {sz[i] <= sz[j]}. The
    procedure reduces {sz[j]} by {sz[i]-1}, and adds {st[j]} to
    {st[i]}; so that {ixB[i]=ixA[i]}, {ixB[j]=ixA[j]+ixA[i]}, and
    {ixB[k]=ixA[k]} for all other axes {k}. If the array was non-empty,
    it will remain so.
    
    Thus, for example, if {M} is a {100×120×5} matrix, the call
    {ix_diagonal(d,sz,&bp,st,0,1)} will chop it down to a diagonal
    band, starting with the main diagonal and extending 21 elements to
    the right, indexed as a {100×21×5} matrix. */

void ix_chop 
  ( ix_dim_t d, 
    ix_size_t sz[], 
    ix_pos_t *bp, 
    ix_step_t st[], 
    ix_axis_t i, 
    ix_size_t stride, 
    ix_axis_t j 
  );
  /* Chops the array perpendicularly to axis {i} into fat slices
    of thickness {stride}, that are stacked along axis {j}. 
  
    Upon entry, the array size {sz[i]} along axis {i} must be a
    multiple of {stride}, and the array must be trivial along axis {j}
    ({sz[j] == 1}). The procedure sets {sz[j] = sz[i]/stride}, {sz[i] =
    stride}, and modifies {st} so that {ixB[i] = (ixA[i] + stride*ixA[j])},
    {ixB[j] = 0}. The effect is undefined if {stride == 0} or {j == i}
    
    For exmple, {ix_chop(d,sz,&bp,st, 0, 100, 1)} will turn a vector of
    5000 elements (along axis 0) into a two-dimensional arryay with
    100 elements per column (along index 0) and 50 columns (along axis 1). */

/* INDEX TUPLE STEPPING */
 
typedef enum ix_order_t
  { ix_order_L = 0, /* Last index is changed first (`C/Pascal' order). */
    ix_order_F = 1  /* The first index is changed first (`Fortran' order). */
  } ix_order_t;
  /* An {ix_order_t} value specifies the order in which index tuples
    are enumerated by some procedures, such as {ix_next} and {ix_enum}. 
    
    The value {ix_order_L} specifies the standard lexicographic order,
    with the LAST index in the innermost loop, as in C and Pascal. If
    {sz} is {{2,3}}, the order is {(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)}.
    
    The value {ix_order_F} specifies flipped lexicographic order,
    with the FIRST index in the innermost loop, as in FORTRAN. If
    {sz} is {{2,3}}, the order is {(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)}. */

bool_t ix_next
  ( ix_dim_t d, 
    ix_index_t ix[],
    const ix_size_t sz[],
    ix_order_t ixor,
    const ix_step_t stA[],
    ix_pos_t *pA, 
    const ix_step_t stB[],
    ix_pos_t *pB, 
    const ix_step_t stC[],
    ix_pos_t *pC
  );
  /* Set {ix[0..d-1]} to the next valid index tuple for arrays with
    size {sz[0..d-1]}, in the order specified by {ixor}.
    
    The procedure also optionally updates the positions {pA,pB,pC} of
    the elements corresponding to the tuple {ix} in up to three arrays
    {A,B,C}, with same size {sz} but independent step vectors
    {stA,stB,stC}.  
    
    For a typical example of its use, see the implementation of
    {ix_enum} below.
    
    More precisely, the procedure scans the indices {ix[0]} to
    {ix[d-1]}, in the fast-to-slow order, looking for an index {ix[i]}
    that is strictly less than its limit {sz[i]-1}. If the procedure
    finds such an index, it increments that index by one, sets every
    previously scanned index {ix[j]} to 0, and returns FALSE. If every
    index {ix[i]} has readed its limit {sz[i]-1}, the procedure sets
    all indices back to {(0,..0)}, and returns TRUE.
    
    The effect is undefined if any size {sz[i]} is 0, or if any index
    {ix[i]} lies outside its valid range {0..sz[i]-1}.
    
    If {pA} and {stA} are not null, the procedure assumes that {*pA}
    is the position corresponding to the index tuple {ix} in some
    array {A} with step vector {stA[0..d-1]}. The procedure then
    adjusts {*pA} to account for the change in the index tuple {ix}.
    In particular, when the procedure returns TRUE, it resets {*pA}
    the position corresponding to {(0,..0)} in {A}.
    
    Similarly, if {pB} and/or {pC} are non-NULL, the procedure assumes
    that {*pB,*pC} are positions corresponding to {ix} into to two
    other arrays {B,C}, with the same indices but independent steps
    vectors {stB,stC}; and updates them in the same way. */
    
bool_t ix_prev
  ( ix_dim_t d, 
    ix_index_t ix[],
    const ix_size_t sz[],
    ix_order_t ixor,
    const ix_step_t stA[],
    ix_pos_t *pA, 
    const ix_step_t stB[],
    ix_pos_t *pB, 
    const ix_step_t stC[],
    ix_pos_t *pC 
  );
  /* Set {ix[0..d-1]} to the previous valid index tuple for arrays with
    size {sz[0..d-1]}, in the order specified by {ixor}.
    
    The procedure also optionally updates the positions {pA,pB,pC} of
    the elements corresponding to the tuple {ix} in up to three arrays
    {A,B,C}, with same size {sz} but independent step vectors
    {stA,stB,stC}.  
    
    This function is basically the inverse of {ix_next}. For a typical
    example of its use, see the implementation of {ix_enum} below.
    
    More precisely, the procedure scans the indices {ix[0]} to
    {ix[d-1]}, in the fast-to-slow order, looking for an index
    {ix[i]} that is strictly positive. If the procedure finds such an
    index, it decrements that index by one, sets every previously
    scanned index {ix[j]} to its upper limit {sz[j]-1}, and returns
    FALSE. If every index {ix[i]} is zero, the procedure sets every
    index {ix[i]} to its upper limit {sz[i]-1}, and returns TRUE.
    
    The effect is undefined if any size {sz[i]} is 0, or if any index
    {ix[i]} lies outside its valid range {0..sz[i]-1}.
    
    If {pA} and {stA} are not null, the procedure assumes that {*pA} is the
    position corresponding to the index tuple {ix} in some array with
    step vector {stA[0..d-1]}. The procedure then adjusts {*pA} to account for the
    change in the index tuple {ix}. In particular, when the procedure
    returns TRUE, it resets {*pA} to the position corresponding to the
    last tuple {(sz[0]-1,sz[1]-1, ..., sz[d-1]-1)}.
    
    Similarly, if {pB} and/or {pC} are non-NULL, the procedure assumes
    that {*pB,*pC} are positions corresponding to {ix} into to two
    other arrays {B,C}, with the same indices but independent steps
    vectors {stB,stC}; and updates them in the same way. */

sign_t ix_compare ( ix_dim_t d, const ix_index_t ixa[], const ix_index_t ixb[], ix_order_t ixor );
  /* Returns -1, 0, or +1 depending on whether the index tuple {ixa}
    comes before, coincides, or follows {ixb} in the order defined by {ixor}.
    The result is undefined if any index is negative. */
    
int ix_sync_level 
  ( ix_dim_t d, 
    const ix_size_t sz[], 
    const ix_index_t ix[], 
    ix_order_t ixor, 
    bool_t reverse
  );
  /* Returns the number of indices in {ix[0..d-1]} that would be set to their
    respective starting values if the index tuple were advanced in
    order defined by {ixor}; either in increasing order (if {reverse} is false) 
    or decreasing order (if [reverse} is true).  
    
    This procedure is useful to detect when the enumeration has reached the
    first or last element of a new slice (row, plane, etc.) of an array. */

/* ELEMENT ENUMERATION */

typedef bool_t ix_index_op_t (const ix_index_t ix[]);
  /* Type of a procedure that operates on a given index tuple {ix[0..d-1]}, returning
    a boolean value.  */

typedef bool_t ix_index_pos_op_t ( const ix_index_t ix[], ix_pos_t pA );
typedef bool_t ix_index_pos2_op_t ( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB );
typedef bool_t ix_index_pos3_op_t ( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC );
  /* Type of a procedure that operates on a given index tuple {ix[0..d-1]}, and corresponding
    array position values {posA}, {posB}, {posC}, returning a boolean value.  */

bool_t ix_enum 
  ( ix_index_pos3_op_t *op,
    ix_dim_t d,
    const ix_size_t sz[],
    ix_order_t ixor,
    bool_t reverse,
    ix_pos_t bpA, 
    const ix_step_t stA[], 
    ix_pos_t bpB, 
    const ix_step_t stB[], 
    ix_pos_t bpC, 
    const ix_step_t stC[]
  );
  /* Applies the operation {op} to all valid index tuples of three arrays 
    {A,B,C} in parallel. 
    
    More precisely, enumerates all valid index tuples {ix[0..d-1]}
    varying the indices in the order defined by {ixor}; from
    {(0,..0)} to {sz-(1,..1)} if {reverse} is FALSE, and from
    {sz-(1,..1)} to {(0,..0)} if reverse is TRUE. For each tuple {ix},
    computes the relative positions {pA,pB,pC} using offsets {bpA,bpB,bpC}
    and steps {stA,stB,stC}, respectively; and calls
    {op(ix,pA,pB,pC)}.
    
    The procedure stops the enumeraton when the call to {op} returns
    {TRUE}, and returns {TRUE}. Otherwise the enumeration continues
    until all valid index tuples are exhausted, and returns {FALSE}. In
    particular, if the array is empty (that is, one of the sizes
    {sz[ax]} is zero), the procedure returns {FALSE} without ever
    calling {op}. 
    
    If {stA} is NULL, {pA} is always {bpA}. Ditto for the other two
    arrays.
    
    Note that any virtually replicated array elements will be visited
    multiple times --- once for each valid index tuple that refers to
    them. A repeated visit to an element of {A} can be recognized by
    {(ix[i]>0) && ((stA==NULL)||(stA[i]==0))} for some {i}; and ditto for
    the other arrays. */

/* INDEXING OF COMPACT ARRAYS

  The elements of an array are often /packed/ in consecutive position of
  a vector, without any gaps, in a specific ordering of the index
  tuples.  In those situations, the steps {st[0..d-1]} are
  determined by the sizes {sz[0..d-1]}.
  
  In the the C/Pascal convention, where the LAST index has
  the smallest increment, the steps are
    
    { st[i] = PROD { sz[j] : j = i+1,..d-1 } }              (2)
    
  In this case, the position of an element can be computed 
  by the `recursive' formula
  
    { p = ix[0]; 
      for (i = 1; i < d; i++) { p = p*sz[i-1] + ix[i]; }
      p += bp;
    }                                                       (3)
  
  Under the FORTRAN convention, where the FIRST index has the smallest
  address increment, the steps are
  
    { st[i] = PROD { sz[j] : j = 0,..i-1 } }                (4)
  
  The position of an element can be computed as 
  
    { p = ix[d-1]; 
      for (i = d-2; i >= 0; i--) { p = p*sz[i+1] + ix[i]; } 
      p += bp;
    }                                                       (5)
  
*/ 

ix_pos_t ix_packed_position 
  ( ix_dim_t d, 
    const ix_index_t ix[], 
    ix_pos_t bp, 
    const ix_size_t sz[], 
    ix_order_t ixor
  );
  /* Computes the position of element {M[ix[0],..ix[d-1]} in an array
    {M} whose elements are stored in consecutive positions, starting at index {bp},
    in the packed indexing order specified by {ixor}.
    
    The procedure returns {ix_pos_NONE} if any {ix[i]} is negative or
    exceeds {sz[i]-1}. */

void ix_packed_steps 
  ( ix_dim_t d, 
    const ix_size_t sz[], 
    ix_order_t ixor, 
    ix_step_t st[]
  );
  /* This procedure computes the position increments {st[0..d-1]} for
    an array with sizes {sz[0..d-1]}, assuming that the positions are
    all distinct and consecutive, in the packed indexing order specified by {ixor}. In
    particular, {st[d-1] == 1} for the C/Pascal order, and {st[0] ==
    1} for the Fortran order.
    
    If any size {sz[i]} is 1, the corresponding increment {st[i]} (which 
    is irrelevant) is set to 0. If any size {sz[i]} is 0 (empty array), 
    all increments {st[0..d-1]} are set to 0. */

void ix_packed_indices 
  ( ix_dim_t d, 
    ix_pos_t p,
    ix_pos_t bp,
    const ix_size_t sz[],
    ix_order_t ixor, 
    ix_index_t ix[]
  );
  /* Computes the index tuple {ix[0..d-1]} that corresponds to element
    position {p} in the packed array with base position {bp} and
    sizes {sz[0..d-1]}. The procedure returns {ix_pos_NONE} if
    there is no such index tuple; in particular, if the array is empty. */

/* PARAMETER VALIDATION */

bool_t ix_parms_are_valid 
  ( ix_dim_t d, 
    const ix_size_t sz[], 
    ix_pos_t bp, 
    const ix_step_t st[], 
    bool_t die
  );
  /* Checks whether the indexing parameters {sz[0..d-1]}, bp, and {st[0..d-1]} 
    satisfy the following conditions:

      (0) Every valid index tuple should yield a non-negative
        position.
        
      (1) For any {i}, if {sz[i]==1}, then {st[i]==0}.

      (2) If any {sz[i]} is zero, then {bp} and all {st}
        coefficients are zero.

      (3) The parameters {bp}, {sz[i]}, {st[i]}, and the number of 
        positions do not exeed their limits (see below).

    The procedure returns TRUE if the parameters are valid. Otherwise,
    it either aborts the program with an error message (if {die=TRUE})
    or returns FALSE (if {die=FALSE}). */

bool_t ix_positions_are_distinct 
  ( ix_dim_t d, 
    const ix_size_t sz[], 
    const ix_step_t st[], 
    bool_t die 
  );
  /* Checks whether the indexing parameters {sz[0..d-1]} and
    {st[0..d-1]} imply distinct positions for distinct valid index
    tuples, except for virtual replication. That is, whether two
    tuples {ixa[0..d-1]} and {ixb[d-1]} give the same position iff and
    only iff {ixa[i]!=ixb[i]} implies {st[i]=0}.
    
    The procedure assumes that the parameters are valid in the sense
    of {ix_parms_are_valid}. It returns TRUE if the distinctness
    condition is satisfied; otherwise, it either aborts the program
    with an error message (if {die=TRUE}) or returns FALSE (if
    {die=FALSE}). */

/* INDEX RANGE REDUCTION */

typedef enum
  { ix_reduction_SINGLE,  /* ... *,*,*,0,1,2,3,4,5,*,*,*,*,*,*,*,*,*,* ... */
    ix_reduction_EXTEND,  /* ... 0,0,0,0,1,2,3,4,5,5,5,5,5,5,5,5,5,5,5 ... */
    ix_reduction_REPEAT,  /* ... 3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3 ... */
    ix_reduction_MIRROR,  /* ... 2,1,0,0,1,2,3,4,5,5,4,3,2,1,0,0,1,2,3 ... */
    ix_reduction_PXMIRR   /* ... 3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,4,5 ... */
  } ix_reduction_t;
#define ix_reduction_FIRST ix_reduction_SINGLE
#define ix_reduction_LAST  ix_reduction_PXMIRR
  /* A mapping from unrestricted {ix_index_t} values (positive or negative) to 
    integers in some range {0..N-1}.  The comments above illustrate each mapping
    for {N=6} and various indices from {-3} to {+15}.  The code '*' denotes {-1}. */

ix_index_t ix_reduce ( ix_index_t i, ix_size_t N, ix_reduction_t red );
  /* Returns the index {i} reduced to {0..N-1} as specified by {red} */

void ix_reduce_range ( ix_index_t i0, ix_size_t m, ix_size_t N, ix_reduction_t red, ix_index_t i[] );
  /* Stores into {i[0..m-1]} the indices {i0..i0+m-1} reduced to {0..N-1} as specified by {red}. */

/* LIMITS */

#define ix_MAX_DIM 40
  /* Should be enough. Note that if no indices are trivial or replicated, 
    an array with {N} indices has at least {2^N} elements. */
     
#define ix_MAX_AXIS (ix_MAX_DIM - 1)
  /* Follows from {ix_MAX_DIM}. */

#define ix_MAX_ABS_STEP (1099511627775LLU)
  /* Since we are using 64-bit integers, we can be generous. We
    arbitrarily set {ix_MAX_ABS_STEP = 2^40-1}, enough to step over 1
    trillion elements but still safely away from overflow. */

#define ix_MAX_POSITIONS (ix_MAX_ABS_STEP + 1)
  /* By combining subsampling with diagonal extraction a client could
    get the step of an array to be as high as {ix_MAX_POSITIONS-1},
    which leads to this limit on ix_MAX_POSITIONS (2^40 > 1 trillion
    positions).
    
    Note that this is a limit on the number of *distinct* element
    positions, not on the number of valid index tuples. The latter can
    exceed this limit if the array has virtually replicated samples
    (i.e. one or more steps {st[i]} equal to zero). */

#define ix_MAX_SIZE (ix_MAX_POSITIONS)
  /* In principle the size of an empty array along some axis could be
    greater than {ix_MAX_POSITIONS}, since that limit applies to the
    total number of element positions (which is zero). However, that
    feature seems to be of little use, and complicates the code. 
    Thus we enforce this limit on individual sizes, too, even for empty
    arrays. */

#define ix_MAX_INDEX (ix_MAX_SIZE - 1)
  /* We need at least this much in order to index every element, and no more
    than this since we cannot create arrays bigger than that.
    Note that {ix_MAX_INDEX + 1} does not overflow.} */

#define ix_MAX_POS (ix_MAX_POSITIONS - 1)
  /* Valid element positions run from 0 to the max number of 
    positions minus 1. */

#endif
