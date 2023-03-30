#ifndef indexing_descr_H
#define indexing_descr_H

/* Multidimensional array descriptors. */
/* Last edited on 2023-03-18 11:00:50 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <ix.h>

/* BASIC TYPES */

/* NUMBER OF AXES IN DESCRIPTOR */

#define ix_descr_MAX_AXES 6
  /* Number of indices (dimensions, axes) in any array descriptor. Lower-dimensional
    arrays are obtained by setting the size to 1 along the unwanted axes. */

/* DESCRIPTORS */

typedef struct ix_descr_t 
  { ix_dim_t na;                      /* Effective number of indices (axes). */
    ix_size_t sz[ix_descr_MAX_AXES];  /* Extent of array's domain along each axis. */
    ix_step_t st[ix_descr_MAX_AXES];  /* Position increment for each index. */
    ix_pos_t bp;                      /* Position of element {[0,..0]}. */
  } ix_descr_t; 
  /* An indexing formula for arrays with {na} indices. The position of
    an element with index tuple {ix[0..na-1]} is 
    
      {D.bp + ix[0]*D.st[0] +·· ix[D.na-1]*D.st[D.na-1]]}
      
    See {indexing.h} for more details. */

ix_descr_t ix_descr_make ( ix_dim_t na, const ix_size_t sz[], const ix_step_t st[], ix_pos_t bp );
  /* Returns a descrptor with the given fields.  Performs no checking. */

/* STEP TUPLE MANIPULATION */

ix_descr_t ix_descr_from_sizes ( ix_dim_t na, const ix_size_t sz[] );
  /* Returns a descrptor record {D} for an array that has size {sz[i]}
    along axis {i}, for {i} in {0..na-1}.
    
    The base position {D.bp} is set to zero. The position increments
    {D.st[0..na-1]} are computed so as to imply a compact set of element
    positions, in lexicographic order. */

void ix_descr_steps_assign ( ix_descr_t *D, ix_step_t sta[], const ix_step_t stb[] );
  /* Sets {sta[i] = stb[i]} for {i=0,...D.na-1}. */

/* ACESSING INDIVIDUAL ELEMENTS */

ix_pos_t ix_descr_position ( ix_descr_t *D, const ix_index_t ix[] );
  /* Returns the position of the element of {D} with indices 
    {ix[0..D.na-1]}. */

/* INDEX TUPLE MANIPULATION */


void ix_descr_indices_fill ( ix_descr_t *D, ix_index_t ix[], ix_index_t val );
  /* Sets {ix[i] = val} for {i=0,..D.na-1}. */

void ix_descr_indices_assign ( ix_descr_t *D, ix_index_t ix[], const ix_index_t val[] );
  /* Sets {ix[i] = val[i]} for {i=0,..D.na-1}. */

void ix_descr_indices_shift ( ix_descr_t *D, ix_index_t ix[], const ix_index_t inc[] );
  /* Sets {ix[i] += inc[i]} for {i=0,..D.na-1}. */

bool_t ix_descr_indices_first ( ix_descr_t *D, ix_index_t ix[] );
  /* Sets {ix[i] = 0} for {i=0,..D.na-1}, and returns TRUE iff that
    is a valid index tuple for {D} (i.e. {D.sz[i]>0} for all {i}). */

bool_t ix_descr_indices_last ( ix_descr_t *D, ix_index_t ix[] );
  /* Sets {ix[i] = D.sz[i]-1} for {i=0,..D.na-1}, and returns TRUE iff that
    is a valid index tuple (i.e. {D.sz[i]>0} for all {i}). */

bool_t ix_descr_indices_are_valid ( ix_descr_t *D, const ix_index_t ix[] );
  /* Returns TRUE iff {ix[0..D.na-1]} is a valid index tuple for {D};
    i.e. if {ix[i]} lies in the range {0..D.sz[i]-1}, 
    for {i=0,..D.na-1}. */

/* SIZE TUPLE MANIPULATION */

void ix_descr_sizes_assign ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = szb[i]} for {i=0,..D.na-1}. */

bool_t ix_descr_sizes_shrink ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = min(sza[i], szb[i])} for {i=0,..D.na-1}.
    Returns TRUE iff any {sza[i]} was actually changed. */

bool_t ix_descr_sizes_expand ( ix_descr_t *D, ix_size_t sza[], const ix_size_t szb[] );
  /* Sets {sza[i] = max(sza[i], szb[i])} for {i=0,..D.na-1}.
    Returns TRUE iff any {sza[i]} was actually changed. */

ix_size_t ix_descr_max_size ( ix_descr_t *D );
  /* Returns the max dimension of {D} along any axis. */

ix_size_t ix_descr_min_size ( ix_descr_t *D );
  /* Returns the min dimension of {D} along any axis in {0..D.na-1}. */

/* ELEMENT COUNTING */

bool_t ix_descr_is_empty ( ix_descr_t *D );
  /* TRUE iff the array {D} is empty, i.e. {D.sz[i]==0} for some {i} in {0..D.na-1}. */

ix_count_t ix_descr_num_tuples ( ix_descr_t *D );
  /* Number of valid index tuples in the array {D}; 0 if the array is empty.
    Beware that virtual element replication may create huge arrays that exceed
    the capacity of an {ix_count_t}. */

ix_count_t ix_descr_num_positions ( ix_descr_t *D );
  /* Number of *distinct* element positions in {D}, not counting virtual replications.
    Assumes that the descriptor is valid (see {ix_descr_is_valid} below). */

ix_pos_t ix_descr_min_pos ( ix_descr_t *D );
ix_pos_t ix_descr_max_pos ( ix_descr_t *D );
  /* Minimum and maximum position of all elements of {D} with a valid
    index tuple. Undefined if the array is empty. */

bool_t ix_descr_same_size ( ix_descr_t *A, ix_descr_t *B, bool_t die );
  /* Returns TRUE if {A.na == B.na} and {A.sz[i] == B.sz[i]} for all
    {i} in {0..A.na-1}, returns TRUE. Otherwise, returns FALSE (if
    {die == FALSE}) or fails with error (if {die == TRUE}). */ 

bool_t ix_descr_contained ( ix_descr_t *D, ix_descr_t *B, bool_t die );
  /* Returns TRUE if the domain of {A} is contained in the domain of
    {B} --- that is, if if {A.na == B.na} {A.sz[i] <= B.sz[i]} for all
    {i} in {0..A.na-1}. Otherwise, returns FALSE (if {die == FALSE})
    or fails with error (if {die == TRUE}). */ 

/* DESCRIPTOR MANIPULATION */

/* The following procedures modify the fields {D.sz}, {*D.bp} and
  {D.st} of the given descriptor {D} so as to change the set of valid
  indices and/or the mapping between indices and positions. Obviously,
  no array elements are actually allocated, reclaimed, or changed.
  
  In the describing comments, {ixA[0..D.na-1]} is a generic tuple valid
  for {D} AFTER the call, referring to some position {p}; and
  {ixB[0..D.na-1]} is the tuple that, BEFORE the call, would yield the
  same position.
  
  All the procedures below preserve the validity of the descriptor
  {D} (as defined by {ix_descr_is_valid} below).  Unless said otherwise,
  the operations preserve the empty/non-empty status of the array. */

void ix_descr_crop ( ix_descr_t *D, ix_axis_t i, ix_size_t skip, ix_size_t keep );
  /* Crops the array described by {D} along axis {i} to the 
    index range {ini..fin}. 
  
    More precisely, sets {D.sz[i]} to {keep}, and modifies {D.bp,D.st}
    so that {ixB[i] = ixA[i]+skip} and {ixB[k] = ixA[k]} for all {k\neq i}.
    If {keep==0}, the array becomes empty (but retains its size along
    the other axes); in that case {D.bp} and {D.st[0..D.na-1]} are reset to 0.
    The result is undefined if {skip+keep > D.sz[i]}.
    
    Thus, to clip a two-dimensional {m × n} matrix to its upper right
    quarter, use {ix_descr_crop(D, 0, 0,m/2)} and then
    {ix_descr_crop(D, 1, n/2,n-n/2)}, where {m=D.sz[0]} and {n=D.sz[1]}. 
    To reduce that matrix to its column number {c}, use
    {ix_descr_crop(D, 1, c,1)}. */

void ix_descr_subsample ( ix_descr_t *D, ix_axis_t i, ix_size_t stride );
  /* Subsamples the array described by {D} along axis {i} by taking
    one every {stride} slices. The sampling step {stride} must be positive. 
    
    More precisely, reduces {D.sz[i]} to {floor((D.sz[i]-1)/stride)+1}, and modifies
    {D.bp} so that {ixB[i] = stride*ixA[i]} and {ixB[k] = ixA[k]} for all {k\neq i}. */

void ix_descr_flip ( ix_descr_t *D, ix_axis_t i );
  /* Flips the array described by {D} along axis {i}. 
  
    More precisely, modifies {D.bp} and {D.st[i]} so that
    {ixB[i] = D.sz[i]-1-ixA[i]}, and {ixB[k] = ixA[k]} for all {k\neq i}. The
    size {D.sz[i]} is not changed. */

void ix_descr_replicate ( ix_descr_t *D, ix_axis_t i, ix_size_t size );
  /* Virtually replicates the array described by {D} along axis {i}, {size} times.
  
    On input, the array must be trivial along axis {i} (i.e., {D.sz[i]} must be 1).
    The procedure sets {D.sz[i]=size} and {D.st[i]=0}, so that {ixB[i]=0} and
    {ixB[k]=ixA[k]} for all {k\neq i}. The replication count {size} must be
    positive. */

void ix_descr_swap_indices ( ix_descr_t *D, ix_axis_t i, ix_axis_t j, ix_dim_t n );
  /* Transposes the array described by {D}, by exchanging the {n} axes starting
    with axis {i} with the {n} axes starting with index {j}. The two sets of indices
    must be either identical ({i==j} or {n==0}, in which case the procedure has
    no effect) or disjoint.
  
    More precisely, swaps {D.sz[i+k]} with {D.sz[j+k]} and {D.st[i+k]} with
    {D.st[j+k]}, for {k} in {0..n-1}; so that {ixB[i+k]=ixA[j+k]}, {ixB[j+k] =
    ixA[i+k]}, for {k} in {0..n-1}, and {ixB[r]=ixA[r]} for every other axis
    {r}. If {i == j}, the procedure has no effect.
    
    Thus, to transpose a two-dimensional matrix, in the ordinary sense, use
    {ix_descr_swap_indices(D, 0,1, 1)}. To reduce it to its row number {r}, viewed
    as a column matrix, use {ix_descr_crop(D, 0, r,1)} and then
    {ix_descr_swap_indices(D, 0,1, 1)}. */

void ix_descr_flip_indices ( ix_descr_t *D, ix_axis_t i, ix_axis_t j );
  /* Reverses the order of all indices of array described by {D} 
    between axis {i} and axis {j}, inclusive.
  
    More precisely: if {i<j}, the procedure swaps {D.sz[i+k]} with
    {D.sz[j-k]} and {D.st[i+k]} with {D.st[j-k]} for all {k} in
    {0..(j-i-1)/2}; thus {ixB[i+k] = ixA[j-k]} for all those {k}. 
    It is a no-op if {i>=j}. */

void ix_descr_slice( ix_descr_t *D, ix_dim_t n, const ix_axis_t ax[], const ix_index_t ix[] );
  /* Reduces {D} to the sub-array of {D} where the {n} indices along axes
    {ax[0..n-1]} are set to the values {ix[0..n-1]}, and then eliminated from
    the array. 
    
    The axes {ax[0..n-1]} must be in strictly increasing order, and must be less
    than {D.na}. If {ax} is NULL, assumes {ax[k] = k} for any {k} in {0..n-1},
    i.e. sets and eliminates the first {n} indices of {D}. It is a no-op if
    {n==0}; otherwise, each index value {ix[k]} must be valid for the axis
    {ax[k]}, so the array must not be empty. The number of indices {D.na} gets
    reduced by {n}.
  
    More precisely: {ixB[kB] = ixA[kB-k]} for {kB} in {0..D.na-1} but not in
    {ax[0..n-1]}, {ixB[kB] = ix[k]} for {kB} in {ax[0..n-1]}, where {k} is the
    number of axes less than {kB} that are in {ax[0..n-1]}. */

void ix_descr_diagonal ( ix_descr_t *D, ix_axis_t i, ix_axis_t j );
  /* Extracts a diagonal slice of the array described by {D}, parallel to the
    bisector of axes {i} and {j}; and redefines axis {i} to be along
    that diagonal. Thus, using index values {r} and {s} for those axes,
    after the call, is equivalent to using {r} and {r+s} before the call.
  
    Upon entry, requires {i != j} and {0 < D.sz[i] <= D.sz[j]}. The
    procedure reduces {D.sz[j]} by {D.sz[i]-1}, and adds {D.st[j]} to
    {D.st[i]}; so that {ixB[i]=ixA[i]}, {ixB[j]=ixA[j]+ixA[i]}, and
    {ixB[k]=ixA[k]} for all other axes {k}. If the array was non-empty,
    it will remain so.
    
    Thus, for example, if {D} is a {100×120×5} matrix, the call
    {ix_descr_diagonal(D,0,1)} will chop it down to a diagonal
    band, starting with the main diagonal and extending 21 elements to
    the right, indexed as a {100×21×5} matrix. */

void ix_descr_chop ( ix_descr_t *D, ix_axis_t i, ix_size_t stride, ix_axis_t j );
  /* Chops the array described by {D} perpendicularly to axis {i} into
    fat slices of thickness {stride}, that are stacked along axis {j}. 
  
    Upon entry, the array size {D.sz[i]} along axis {i} must be a
    multiple of {stride}, and the array must be trivial along axis {j}
    ({D.sz[j] == 1}). The procedure sets {D.sz[j] = D.sz[i]/stride},
    {D.sz[i] = stride}, and modifies {st} so that {ixB[i] = (ixA[i] +
    stride*ixA[j])}, {ixB[j] = 0}. The effect is undefined if {stride
    == 0} or {j == i}
    
    For exmple, {ix_descr_chop(D, 0, 100, 1)} will turn a vector of
    5000 elements (along axis 0) into a two-dimensional array with
    100 elements per column (along index 0) and 50 columns (along axis 1). */

/* INDEX TUPLE STEPPING */

bool_t ix_descr_next ( ix_descr_t *D, ix_index_t ix[], ix_pos_t *p );
  /* Increments the index tuple {ix} to the next valid tuple in
    lexicographic indexing order.
    
    More precisely, the procedure scans the indices {ix[D.na-1]} to
    {ix[0]}, in that order, looking for an index {ix[i]}
    that is strictly less than its limit {D.sz[i]-1}. If the procedure
    finds such an index, it increments that index by one, sets every
    previously scanned index {ix[j]} to 0, and returns FALSE. If every index
    {ix[i]} has readed its limit {D.sz[i]-1}, the procedure sets all
    indices back to {(0,..0)}, and returns TRUE.
    
    If {p} is not null, assumes that {*p} is the position corresponding
    to the index tuple {ix}, and adjusts it to account for the change.
  
    This procedure can be used to vary the index tuple {ix[0..D.na-1]}
    over all the elements of an array with descriptor {D}, in a single
    loop, e.g.
    
      | if (ix_descr_assign_min(D,ix))
      |   { ix_pos_t p = ix_descr_position(D,ix);
      |     do { ... } while (! ix_descr_next(D,ix,ixor,&p));
      |   }

    */
    
bool_t ix_descr_prev ( ix_descr_t *D, ix_index_t ix[], ix_pos_t *p );
  /* Like {ix_descr_next}, but in the reverse order. 
    
    More precisely, the procedure scans the indices {ix[D.na-1]} to
    {ix[0]}, in that order, looking for an index
    {ix[i]} that is strictly positive. If the procedure finds such an
    index, it decrements that index by one, sets every previously
    scanned index {ix[j]} to its upper limit {D.sz[j]-1}, and returns
    FALSE. If every index {ix[i]} is zero, the procedure sets every
    index {ix[i]} to its upper limit {D.sz[i]-1}, and returns TRUE.
    
    If {p} is not null, assumes that {*p} is the position corresponding
    to the index tuple {ix}, and adjusts it to account for the change.
    
    Thus, to scan an array in reverse order, use
    
      | if (ix_descr_assign_max(D,ix,ixor))
      |   { ix_pos_t p = ix_descr_position(D,ix);
      |     do { ... } while (! ix_descr_prev(D,ix,ixor,&p));
      |   }

    */

sign_t ix_descr_compare ( ix_descr_t *D, const ix_index_t ixa[], const ix_index_t ixb[] );
  /* Returns -1, 0, or +1 depending on whether the index tuple {ixa}
    comes before, coincides, or follows {ixb} in lexicographic order
    The result is undefined if any index is negative. */
    
/* ELEMENT ENUMERATION */

typedef ix_index_pos3_op_t ix_descr_op_t;
  /* Client procedure for {ix_descr_enum}. */

bool_t ix_descr_enum 
  ( ix_descr_op_t *op,
    ix_order_t ixor,
    bool_t reverse,
    ix_descr_t *A,
    ix_descr_t *B,
    ix_descr_t *C
  );
  /* Applies the operation {op} to all valid index tuples 
    (and corresponding positions) of three arrays {A,B,C} in parallel.
    The three arrays must have the same number of indices, that is,
    {na == A.na == B.na == C.na}.
    
    More precisely, let {sz[i]} be the minimum of
    {A.sz[i],B.sz[i],C.sz[i]}, for {i} in {0..na-1}. The procedure
    enumerates all valid index tuples {ix[0..na-1]} between {(0,..0)}
    and {sz-(1,..1)}. For each tuple {ix}, it computes the the
    corresponding positions {pA,pB,pC} in {A,B,C}, respectively; and
    calls {op(ix,pA,pB,pC)}.
    
    The procedure stops the enumeraton when the call to {op} returns
    {TRUE}, and returns {TRUE}. Otherwise the enumeration continues
    until all valid index tuples are exhausted, and returns {FALSE}. In
    particular, if the array is empty (that is, one of the sizes
    {sz[ax]} is zero), the procedure returns {FALSE} without ever
    calling {op}.
    
    If {A} is NULL, {A.sz[i]} is assumed to be {+oo}, and {pA} is
    always 0. Ditto for the other two descriptors. If all three
    descriptors are NULL, the procedure is a no-op.
    
    The index tuples are enumerated by varying the indices
    in the order specified by {ixor}, as in {ix_descr_next}
    (if {reverse} is FALSE) or {ix_descr_prev} (if {reverse} is TRUE).
    Use {ix_descr_flip} to scan some or all indices in reverse.
    
    Note that virtually replicated array elements will be visited
    multiple times (once for each valid index tuple that refers to
    them). A repeated visit to an element of {A} can be recognized by
    {(ix[i]>0) && (A != NULL) && (A.st[i]==0)} for some {i}; and ditto for
    the other two descriptors. */

/* PARAMETER VALIDATION */

bool_t ix_descr_is_valid ( ix_descr_t *D, bool_t die );
  /* Checks whether the descriptor {D} satisfies the following conditions:

      (0) Every valid index tuple yields a non-negative  position.
        
      (1) For any {i}, if {D.sz[i]==1}, then {D.st[i]==0}.

      (2) If any {D.sz[i]} is zero, then {D.bp} and all 
        increments {D.st[0..D.na-1]} are zero.

    The procedure returns TRUE if these conditions hold. Otherwise,
    it either aborts the program with an error message (if {die=TRUE})
    or returns FALSE (if {die=FALSE}). */

bool_t ix_descr_positions_are_distinct ( ix_descr_t *D, bool_t die );
  /* Checks whether the indexing parameters of {D} imply distinct
    positions for distinct valid index tuples, except for virtual
    replication. That is, whether two tuples {ixa[0..D.na-1]} and
    {ixb[D.na-1]} give the same position iff and only iff
    {ixa[i]!=ixb[i]} implies {D.st[i]=0}.
    
    The procedure assumes that the descriptor valid in the sense
    of {ix_descr_is_valid}. It returns TRUE if the distinctness
    condition is satisfied; otherwise, it either aborts the program
    with an error message (if {die=TRUE}) or returns FALSE (if
    {die=FALSE}). */

/* INPUT/OUTPUT */

void ix_descr_print_descr ( FILE *wr, char *pre, ix_descr_t *D, int32_t wd, char *suf );
  /* Writes the descriptor {D} to {wr}, in readable ASCII format. 
    The output consists of four lines 
    
      | "axes = {naeff}" 
      | "base = {D.bp}"
      | "step = {D.st[0]} {D.st[1]} ... {D.st[naeff-1]}"
      | "size = {D.sz[0]} {D.sz[1]} ... {D.sz[naeff-1]}"
    
    where {naeff} is the largest integer such that {D.sz[naeff-1] != 1}.
    
    Each line is prefixed by the string {pre} (which defaults to "" if NULL) 
    and suffixed with the string {suf} (which defaults to "\n" if NULL).  Each size or
    step is padded on the left to total width {wd} or more. Elements are separaed by 
    blanks. Positive {step} elements will have an explicit "+" sign. */

#endif
