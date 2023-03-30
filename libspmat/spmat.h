#ifndef spmat_H
#define spmat_H
/* Generic sparse matrix types, and operations thereon */

#define spmat_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:46:14 by stolfi */

/* 
  !!! change sort_entries so that (+1,+2) means row-by-row ? !!!
  !!! define {PREFIX##_extract_element} !!!

  !!! unify {PREFIX##_add_{row,column}} into {PREFIX##_add_elements} with {drow,dcol} bits ? !!!
  !!! Is there a way to define automatically {PREFIX##_MAX_SIZE}, etc.? !!!

  !!! {PREFIX##_add_row} should store {val[j]} into {M[row,col+j]} !!!
  !!! {PREFIX##_add_row} should accept {nv != M.cols} !!!
  !!! define also {PREFIX##_fill_row} !!!

  !!! {PREFIX##_add_col} should store {val[i]} into {M[row+i,col]} !!!
  !!! {PREFIX##_add_col} should accept {nv != M.rows} !!!
  !!! define also {PREFIX##_fill_col} !!!

  !!! define {PREFIX##_add_block} that adds a submatrix. !!! 

  !!! define also {PREFIX##_find_element_sorted(M, row, col, srow, scol)} !!!
  !!! define also {PREFIX##_flip(M, frow, fcol)} !!!
  !!! define also {PREFIX##_rotate(M, qturns)} !!!
  !!! define also {PREFIX##_shift(M, row, col)} !!!
  !!! define also {PREFIX##_scatter(M, srow, scol)}, etc. !!!
  !!! define also {PREFIX##_scan(M, element_visit_function_t)} !!!
*/

/* SPARSE MATRIX DESCRIPTORS

  In applied mathematics, a /sparse matrix/ is a two-dimensional
  matrix that contains a large proportion of /trivial elements/ ---
  elements that are equal to zero, or to some other
  application-specific value.
  
  This interface defines a /packed representation/ for sparse
  matrices, namely a list of the element values together with their
  row and column indices in the original (/unpacked/) representation.
  Any element not on this list is assumed to be trivial. The
  representation is handled through a /descriptor/, a record
  containing the nominal dimensions (row and column counts) of the
  matrix, the count of elements stored in the list, and a pointer to
  an array of (row,col,value) triplets.

  USING SPARSE MATRICES
  
  A new sparse matrix type, with elements of a specific type, is
  declared with the {spmat_typedef} macro, and implemented
  with {spmat_impl} macro. These macros also declare and
  implement the basic functions for handling such sparse matrices.
  Additional functions are available through the interfaces
  {spmat_io.h} and {spmat_linalg.h}.

  The incremental creation of sparse matrices is both easy and
  efficient with the tools provided here. For example, suppose that
  the type {dspmat_t} was previously defined as a sparse matrix of
  {double} values. The following bit of code creates such a matrix,
  with 50 rows and 200 columns, and fills it one element at a time:
  
    ------------------------------------------------------------
    // Create a {50 × 200} sparse matrix {M}
    dspmat_t M = dspmat_new(50,200,400); // Guesing ~400 entries.
    
    // Store some elements into {M}, in arbitrary order:
    int32_t count = 0;
    while(! finished(...))
      { int32_t i = ...;
        int32_t j = ...;
        double Mij = ...; 
        count = dspmat_add_element(&M, count, i, j, Mij);
      }
    
    // Reclaim unused entries:
    dspmat_trim(&M, count); 
    
    // Sort the entries by rows:
    dspmat_sort_entries(&M, +1, 0);
    
    // Print the elements, one row per line:
    int32_t r = -1; // Current row index.
    for (k = 0; k < M.ents; k++)
      { dspmat_entry_t *ek = &(M.e[k]);
        if (ek->row != r) { r = ek->row; printf("\n [%d] :", r); }
        printf(" [%d]=%f", ek->col, ek->val);
      }
    printf("\n");
    ------------------------------------------------------------
    
  The function {dspmat_add_element} above will automatically reallocate
  {M}'s entry list, as needed, to accomodate the new entries. Each
  expansion doubles the size of the entry list, so the total cost of
  multiple expansions is small and proportional to the number of
  elements added. The call {dspmat_trim} truncates the list to the
  specified size, freeing the unused space.
*/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>

/* DECLARING A NEW SPARSE MATRIX TYPE

  Here is how one would define the type {dspmat_t} as a sparse matrix of
  doubles:
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <spmat.h>
    
    spmat_typedef(dspmat_t, dspmat, double);
    ------------------------------------------------------------
    
    ------------------------------------------------------------
    / * Implementations  * /

    #define dspmat_trivial_elem (0.0)
    #define dspmat_elem_is_trivial(X) ((X)==0.0)
    
    spmat_impl(dspmat_t, dspmat, double);
    ------------------------------------------------------------

  As another example, here is how one would define a type {graph_t}
  as being a sparse directed graph, represented by a boolean adjacency
  matrix where most of the entries are FALSE (meaning `no edge'):
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <spmat.h>
    
    spmat_typedef(graph_t, graph, bool_t);
    ------------------------------------------------------------

    ------------------------------------------------------------
    / * Implementations  * /

    #define graph_trivial_elem (FALSE)
    #define graph_elem_is_trivial(X) ((X)==FALSE)
    
    spmat_impl(graph_t, graph, bool_t);
    ------------------------------------------------------------ */

#define spmat_typedef(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_TYPEDEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/* 
  This macro defines a data type called {MATRIX_TYPE}, which is a
  record that describes a sparse matrix with element values are of
  type {ELEM_TYPE}. The operations on those matrices will have
  names beginning with {PREFIX}.
  
  This macro will also define the type {PREFIX##_entry_t} (a triplet
  with fields {.row,.col,.value}), and will generate prototype
  declarations for the following procedures:  
  
    ------------------------------------------------------------
    PREFIX##_new
    PREFIX##_expand
    PREFIX##_trim
    PREFIX##_make_desc
    PREFIX##_copy
    
    PREFIX##_find_element
    PREFIX##_add_element
    
    PREFIX##_add_row
    PREFIX##_extract_row
    PREFIX##_scan_row
    
    PREFIX##_add_col
    PREFIX##_extract_col 
    PREFIX##_scan_col
    
    PREFIX##_add_diagonal
    PREFIX##_fill_diagonal
    
    PREFIX##_sort_entries
    PREFIX##_transpose
    ------------------------------------------------------------
      
  These procedures are described below.
     
  The names {MATRIX_TYPE} and {PREFIX} can be chosen quite arbitrarily.
  However, it is recommented to use {PREFIX##_t} as the {MATRIX_TYPE}.
  For example, 
  
    ------------------------------------------------------------
    spmat_typedef(mymatrix_t, mymatrix, myelem_t)
    spmat_typedef(cmat_t, cmat, complex)
    spmat_typedef(graph_t, graph, bool_t)
    ------------------------------------------------------------
 */
 
#define spmat_impl(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the implementations of the functions whose
  prototypes were declared with {spmat_typedef}.
  
  BEFORE calling this macro, the client must also declare the following 
  macros:
    
    ------------------------------------------------------------
    #define PREFIX##_trival_elem (...)
      / * The trivial element value. * / 

    #define PREFIX##_elem_is_trivial(X) (...)
      / * A boolean expression that yields TRUE iff (X) is a
        trivial element value. * / 
    ------------------------------------------------------------  
    
  The macro {PREFIX##_elem_is_trivial(X)} is used by the sparse
  matrix functions to detect elements that need not be stored
  explicitly, and often has the obvious definition
  
    ------------------------------------------------------------
    #define mymatrix_elem_is_trivial(X) ((X)== mymatrix_trivial_elem)
    ------------------------------------------------------------
  
  (Note the extra '()'s around the 'X' to ensure proper C parsing.)
  However, for matrices with floating-point elements, one might want
  to use instead
  
    ------------------------------------------------------------
    #define mymatrix_elem_is_trivial(X) (fabs(X) < TINY)
    ------------------------------------------------------------
  
  where {TINY} is the expected magnitude of roundoff 
  errors.  For matrices of strings, one may want instead
  
    ------------------------------------------------------------
    #define mymatrix_elem_is_trivial(X) (strlen(X) == 0)
    ------------------------------------------------------------
   
    and so on. */
    
/* ======================================================================

  The remainder of this interface describes the functions provided
  by the macros {spmat_typedef} and {spmat_impl}.  */

/* DATA TYPES */

typedef uint32_t spmat_size_t;
  /* Type of a count of columns or rows in a sparse matrix. */
  
typedef uint32_t spmat_index_t;
  /* Type of column or row index into a sparse matrix. */

typedef uint32_t spmat_count_t;
  /* Type of a count of non-trivial entries in a sparse matrix. */

typedef uint32_t spmat_pos_t;
  /* Type of an index of an entry in a sparse matrix. */

/* ------------------------------------------------------------
  typedef spmat_size_t PREFIX##_size_t;
  typedef spmat_index_t PREFIX##_index_t;
  typedef spmat_count_t PREFIX##_count_t;
  typedef spmat_pos_t PREFIX##_pos_t;
------------------------------------------------------------ */
  /* These types are alternative names for {spmat_size_t},
    {spmat_index_t}, etc. */

/* ------------------------------------------------------------
  typedef ELEM_TYPE PREFIX##_elem_t;
------------------------------------------------------------ */
  /* Another name for the type of the matrix elements. */

/* ------------------------------------------------------------
typedef struct MATRIX_TYPE
  { spmat_size_t rows;
    spmat_size_t cols;
    PREFIX##_entry_t *e;
    spmat_count_t ents;
  } MATRIX_TYPE;
------------------------------------------------------------ */
  /* A variable {M} of this type is a descriptor for a sparse matrix
    of {ELEM_TYPE} values that has {M.rows} rows and {M.cols} columns
    in the unpacked form, but has only {M.ents} explicitly stored
    elements, namely the {.val} fields of the entries
    {M.e[0..M.ents-1]}.

    Every entry {t} in {M.e[]} must have {t.row < M->rows} and {t.col
    < M->cols}. The entries may be in any order, but some uses may
    require that there must not be two entries with the same {row,col}
    indices. The address {.e} may be NULL if (and only if) {M.ents == 0}.*/

/* ------------------------------------------------------------
typedef struct PREFIX##_entry_t
  { spmat_index_t row;
    spmat_index_t col;
    ELEM_TYPE val
  } PREFIX##_ENTRY_t;
------------------------------------------------------------ */
  /* A variable {e} of this type represents an element of the sparse
    matrix. It contains the indices {e.row} and {e.col} of the 
    element, and its value {e.val} (which is usually non-trivial). */

/* MATRIX ALLOCATION */

/* ------------------------------------------------------------
MATRIX_TYPE PREFIX##_new
  ( spmat_size_t rows, 
    spmat_size_t cols, 
    spmat_count_t ents
  );
------------------------------------------------------------ */
  /* Returns a record of type {MATRIX_TYPE}, the descriptor of a newly
    allocated sparse matrix of {ELEM_TYPE} elements. The matrix will
    have {rows} rows and {cols} columns if fully expanded, and its
    entry vector will have space for {ents} non-trivial elements. (All
    three paramters can be modified afterwards; see below.)

    A procedure that modifies a sparse matrix {M} usually expects {M}
    to be properly initialized, with {M->e} set to a valid address (or
    NULL) and {M->ents} set to the size of {M->e[]}. This is true even
    for procedures that completely overwrite their arguments, like
    {PREFIX##_copy}. For that reason, it is a good practice to
    initialize every sparse matrix variable as soon as it is declared
    --- if nothing else, with {PREFIX##_new(0,0,0)}.*/

/* ------------------------------------------------------------
void PREFIX##_expand(MATRIX_TYPE *M, spmat_pos_t pos);
------------------------------------------------------------ */
  /* Makes sure that the entry {M->e[pos]} exists, by reallocating {M->e}
    if necessary.

    Namely, if {pos < M->ents}, the procedure does nothing. If {pos >=
    M->ents}, the procedure allocates a larger entry vector from the
    heap, copies the old entries of {M} into it, reclaims the old
    vector {M->e}, makes {M->e} point to the new area, and sets
    {M->ents} to its size. In that case, the new vector will have at
    least {pos+1} entries, and will be at least twice the size of the
    old one.*/

/* ------------------------------------------------------------
void PREFIX##_trim(MATRIX_TYPE *M, spmat_count_t ents);
------------------------------------------------------------ */
  /* Reallocates the entry area {M->e}, if necessary, so that
    it has exactly {ents} entries.

    Namely, if {ents == M->ents}, the procedure does nothing. If {ents
    != M->ents}, the procedure allocates a new entry vector from the
    heap, with precisely {ents} entries, copies into it the entries
    {M.e[0..size-1]}, reclaims {M->e}, makes {M->e} point to the new
    vector, and sets {M.ents = ents}.

    If the new size {ents} is zero, the address {M->e} is set
    to NULL.*/

/* ------------------------------------------------------------
MATRIX_TYPE PREFIX##_make_desc
  ( spmat_size_t rows,
    spmat_size_t cols,
    ELEM_TYPE *e,
    spmat_count_t ents
  );
------------------------------------------------------------ */
  /* Assembles a sparse matrix descriptor from the row count {rows}, the
    column count {cols}, the entry count {ents}, and the address {e} of
    the first explicit entry. The client must make sure that the entries
    {e[0..ents-1]} actually exist. The procedures {PREFIX##_expand} and
    {PREFIX##_trim} can be applied to the resulting descriptor only if
    the address {e} was obtained through {malloc}. */

/* COPYING */

/* ------------------------------------------------------------
void PREFIX##_copy(MATRIX_TYPE *M, MATRIX_TYPE *N);
------------------------------------------------------------ */
  /* Stores into {N} a copy of {M}.

    More precisely, it makes sure that {N} hs exactly {M->ents} entries,
    (with {PREFIX##_trim(N, M->ents)}), then copies all entries from
    {M->e[]} to {N->e[]}, and sets {N->rows=M->rows}, {N->cols=M->cols}.
    Note that {N} must have been properly initialized. */

/* FINDING ELEMENTS */

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_find_element
  ( MATRIX_TYPE *M,
    spmat_pos_t posIni,
    spmat_pos_t posLim,
    spmat_index_t row,
    spmat_index_t col
  );
------------------------------------------------------------ */
  /* Scans {M->e[posIni..MIN(posLim,M->ne)-1]} looking for an entry {M->e[k]}
    with indices {row,col}.  If the entry exists, returns its position
    {k}; otherwise returns {MIN(posLim,M->ne)}.   */

/* ADDING AND EXTRACTING SINGLE ELEMENTS */

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_add_element
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    spmat_index_t col,
    ELEM_TYPE val
  );
------------------------------------------------------------ */
  /* Stores the element value {val} with indices {row,col} into entry
    {M->e[pos]} of matrix {M}, expanding {M} as needed.

    More precisely, if {row >= M->rows}, the procedure sets {M->rows =
    row+1}; if {col >= M->cols}, it sets {M->cols = cols+1}. Then, if
    the element value {val} is not trivial, the procedure sets
    {M->e[pos].row = row}, {M->e[pos].col = col}, and {M->e[pos].val =
    val}.  In that case, if {pos >= M->ents}, the vector {M.e} is
    automatically reallocated (with {PREFIX##_expand})
    so that {M->e[pos]} exists.

    The function will return {pos} itself if {val} is trivial, or
    {pos+1} otherwise. The client must make sure that the matrix {M}
    does not contain any other entry with indices {row,col}. */

/* ADDING, EXTRACTING, AND SACNNING ELEMENTS BY ROWS OR COLUMNS */

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_add_row
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    ELEM_TYPE val[],
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Stores {val[0..nv-1]} into {M->e}, starting at position
    {M->e[pos]}, as elements of row {row} of the matrix.

    More precisely, if {row >= M.rows}, then {M.rows} is set to {row+1};
    if {nv > M.cols}, then {M.cols} is set to {nv}. Then, each
    non-trivial element {val[j]} is stored into {M->e}, as the entry in
    row {row} and column {j} of the matrix.

    The vector {M->e} will be reallocated as needed (with
    {PREFIX##_expand}) to accomodate the new entries.

    The function returns the index of the entry of {M.e} immediately
    after the last entry that was assigned; or {pos} itself, if all
    elements of {val} were trivial. The procedure does not check
    whether {M} already contains any entries with the same indices as
    the stored elements.*/

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_add_col
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t col,
    ELEM_TYPE val[],
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Stores {val[0..nv-1]} into {M->e}, starting at position
    {M->e[pos]}, as elements of column {col} of the matrix.

    More precisely, if {col >= M.cols}, then {M.cols} is set to {col+1};
    if {nv > M.rows}, then {M.rows} is set to {nv}. Then, each
    non-trivial element {val[i]} is stored into {M->e}, as the entry in
    row {i} and column {col} of the matrix.

    The vector {M->e} will be reallocated as needed (with
    {PREFIX##_expand}) to accomodate the new entries.

    The function returns the index of the entry of {M.e} immediately
    after the last entry that was assigned; or {pos} itself, if all
    elements of {val} were trivial. The procedure does not check
    whether {M} contains any entries with the same indices as the
    stored elements. */
    
/* ------------------------------------------------------------
spmat_pos_t PREFIX##_extract_row
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    ELEM_TYPE val[],
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Expands row {row} of {M}, which starts at {M->e[pos]}, and stores
    it into the ordinary vector {val[0..nv-1]}.

    The function assumes that all entries of {M} with row index {row}
    are stored in consecutive positions of {M->e}, starting with entry
    {M->e[pos]}. Thus, if {pos >= M->ents} or {M->e[pos].row != row},
    the function assumes that all elements of that row are trivial.

    The size {nv} of {val} must be equal to {M->ncols}. The function
    returns the index of the entry in {M->e} immediately after the last
    entry that was actually extracted; or {pos} itself, if that row was
    completely trivial.
    
    Warning: if the extracted row contains two or more entries with
    the same column index, only the last one will be stored into {val}. */

/* ------------------------------------------------------------
void PREFIX##_extract_col
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t col,
    ELEM_TYPE val[],
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Expands column {col} of {M}, which starts at {M->e[pos]}, and stores
    it into the ordinary vector {val[0..nv-1]}.

    The function assumes that all entries of {M} with column index
    {col} are stored in consecutive positions of {M->e}, starting with
    entry {M->e[pos]}. Thus, if {pos >= M->ents} or {M->e[pos].col !=
    col}, the function assumes that all elements of that column are
    trivial.

    The size {nv} of {val} must be equal to {M->nrows}. The function
    returns the index of the entry in {M->e} immediately after the last
    entry that was actually extracted; or {pos} itself, if that column was
    completely trivial.
    
    Warning: if the extracted column contains two or more entries with
    the same row index, only the last one will be stored into {val}.

    Here is a typical example of usage of these procedures:

      ------------------------------------------------------------
      dspmat_sort_entries(&M, +1, 0);  // Sort entries of {M} by row.
      dspmat_t R = dspmat_new(M->rows,M->cols,0);  // Matrix for result.
      double v[M.cols];
      int32_t pM = 0, pR = 0; // Scan the entries of {M} and {R}.
      int32_t i;              // Scans the rows of {M} and {R}.
      for (i = 0; i < M.rows; i++)
        { pM = dspmat_extract_row(&M, pM, i, v);
          for (j = 0; j < M.cols; j++) { modify(v[j]); }
          pR = dspmat_add_row(&R, pR, i, v);
        }
      dspmat_trim(&R, pR);
      ------------------------------------------------------------
  */

/* ------------------------------------------------------------
typedef ELEM_TYPE PREFIX##_entry_scan_proc 
  ( spmat_index_t row, 
    spmat_index_t col, 
    ELEM_TYPE val
  );
  
spmat_pos_t PREFIX##_scan_row
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    PREFIX##_entry_scan_proc proc
  );
------------------------------------------------------------ */
  /* Calls {proc(row, col, val)} for every non-trivial element in
    row {row} of {M}, where {val} is the element's value,
    and {col} its column ndex.

    The function assumes that all entries of {M} with row index {row} are
    stored in consecutive positions of {M->e}, starting with entry
    {M->e[pos]}. Thus, if {pos >= M->ents} or {M->e[pos].row != row},
    the function assumes that all elements of that row are trivial.
    
    The procedure returns the index of the entry in {M->e} immediately
    after the last entry scanned (that is, the end of that row in {M.e}). */

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_scan_col
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t col,
    PREFIX##_entry_scan_proc proc
  );
------------------------------------------------------------ */
  /* Calls {proc(row, col, val)} for every non-trivial element in
    column {col} of {M}, where {val} is the element's value,
    and {row} its row index.

    The function assumes that all entries of {M} with column index {col} are
    stored in consecutive positions of {M->e}, starting with entry
    {M->e[pos]}. Thus, if {pos >= M->ents} or {M->e[pos].col != col},
    the function assumes that all elements of that column are trivial.
    
    The procedure returns the index of the entry in {M->e} immediately
    after the last entry scanned (that is, the end of that column in {M.e}). */
    
/* ADDING ELEMENTS ALONG DIAGONALS*/

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_add_diagonal
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    spmat_index_t col,
    ELEM_TYPE val[],
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Stores {val[0..nv-1]} into {M->e}, starting at position
    {M->e[pos]}, as a diagonal sequence of elements of {M}. The sequence
    begins at row {row} and column {col} of the matrix, and moves
    diagonally down and to the right, wrapping around the edges of {M}
    if necessary.

    More precisely each non-trivial value {val[k]} is stored into {M->e},
    as the entry in row {(row+k)%M->rows} and column {(col+k)%M->cols}
    of the matrix.

    The vector {M->e} will be reallocated as needed (with
    {PREFIX##_expand}) to accomodate the new entries. The matrix
    dimensions {M->rows} and {M->cols} are not affected.

    The function returns the index of the entry of {M.e} immediately
    after the last entry that was assigned; or {pos} itself, if all
    elements of {val} were trivial.
    
    The procedure does not check whether this poeration leaves {M}
    with two or more entries with the same indices. Note that this may
    happen if {nv} exceeds the least common multiple of {M.rows} and
    {M.cols}. */

/* ------------------------------------------------------------
spmat_pos_t PREFIX##_fill_diagonal
  ( MATRIX_TYPE *M,
    spmat_pos_t pos,
    spmat_index_t row,
    spmat_index_t col,
    ELEM_TYPE val,
    spmat_size_t nv
  );
------------------------------------------------------------ */
  /* Similar to {}, but uses the single value {val},
    replicated {nv} times. */

/* SORTING ENTRIES BY INDICES*/

/* ------------------------------------------------------------
void PREFIX##_sort_entries(MATRIX_TYPE *M, int32_t orow, int32_t ocol);
------------------------------------------------------------ */
  /* Rearranges the entries {M->e[0..M->ents-1]} in an order that depends
    on the integers {orow} and {ocol}. The function swaps whole entries
    only; it does not modify any individual fields.

    If {orow} is nonzero, the entries will be sorted according to their
    {.row} indices: increasing if {orow>0}, or decreasing if {orow<0}.
    If {orow} is zero, the row index is not relevant for the ordering.
    The {ocol} parameter refers to the column index {.col}, in the same
    way.

    The absolute values of {orow} and {ocol} (which must be different)
    specify the order of the comparisons: a larger value means that the
    corresponding index is more important.

    Thus, for example, {mymatrix_sort_entries(&M, +1, 0)} will sort the
    entries by inreasing row index, but will leave the entries within
    each row in arbitrary order; while {mymatrix_sort_entries(&M, +1,
    -2)} will sort them by decreasing column index, and, within each column,
    by ascending row index. */

/* ------------------------------------------------------------
void PREFIX##_sort_entries_ins
  ( MATRIX_TYPE *M, 
    int32_t orow, 
    int32_t ocol, 
    spmat_pos_t posIni, 
    spmat_pos_t posLim
  );
------------------------------------------------------------ */
  /* Similar to {PREFIX##_sort_entries}, but uses insertion sort, and
    only applies to the entries {M->e[posIni..posLim-1]}. */

/* TRANSPOSITION */

/* ------------------------------------------------------------
void PREFIX##_transpose(MATRIX_TYPE *M);
------------------------------------------------------------ */
  /* Exchanges the {.row} and {.col} indices of all
    entries of {M}, and also exchanges {M->rows} with {M->cols}.

    It does not change the order of the entries within {M.e}; so, for
    example, if {M} was sorted by rows before the call, it will be
    sorted by columns after it. */

/* CONDENSING AND MERGING ENTRIES */

/* ------------------------------------------------------------
typedef void PREFIX##_entry_merge_proc_t 
  ( spmat_index_t row, 
    spmat_index_t col, 
    ELEM_TYPE valA, 
    ELEM_TYPE valB
  );

void PREFIX##_merge(MATRIX_TYPE *A, MATRIX_TYPE *B, PREFIX##_entry_merge_proc_t proc);
------------------------------------------------------------ */
  /* Calls {proc(i,j,A[i,j],B[i,j])} for all corresponding pairs
    of elements of {A} and {B} where at least one of them 
    is explicitly represented.

    The entries of {A} MUST BE SORTED by increasing row index, then by
    increasing column index within each row (as produced by
    {PREFIX##_sort_entries(&A, +2, +1)}); and ditto for {B}.
    The pairs are visited in the same order.

    The matrices {A} and {B} must also have the same dimensions
    {A.rows==B.rows} and {A.cols==B.cols}.
    
    If the {A} matrix contains {na} entries with the same {.row} and
    {.col} indices, and the {B} matrix contains {nb} entries with
    those same indices, then {proc} is called {max(na,nb)} entries on
    those entries. Namely, each of those {A} entries is paired with
    the corresponding {B} entry, in order of occurrence, the shortest
    list being implicitly completed with trivial values. */

/* ------------------------------------------------------------
typedef ELEM_TYPE PREFIX##_entry_condense_proc_t 
  ( spmat_index_t row, 
    spmat_index_t col, 
    ELEM_TYPE v0, 
    ELEM_TYPE v1
  );

void PREFIX##_condense(MATRIX_TYPE *A, PREFIX##_entry_condense_proc_t proc);
------------------------------------------------------------ */
  /* Combines any set of two or more consecutive entries of {A}
    which have the same column and row indices into a single entry,
    by successive calls to {proc}.  Discards trivial entries.
    
    Namely, if the values of those entries are {v[0], v[1], ... v[n-1]},
    then the value {vt} of the combined entry will be obtained by 
    {vt = v0; vt = proc(vt, v[1]); vt = proc(vt, v[2]); ... vt = proc(vt,v[n-1])}.
    The combined entry is retained in {A} only if {vt} is not trivial.
    
    The matrix is automatically trimmed if the number of entries
    decreases.
  
    If, upon entry, the entries of {A} are arranged so that all entries with the
    same indices appear in consecutive positions, then, upon exit, {A}
    will have at most one entry (with non-trivial value) for each 
    index pair. */

/* AUXILIARY FUNCTIONS */

int32_t spmat_compare_indices
  ( spmat_index_t arow, 
    spmat_index_t acol, 
    spmat_index_t brow, 
    spmat_index_t bcol, 
    int32_t orow, 
    int32_t ocol
  );
  /* Compares two index pairs {arow,acol} and {brow,bcol}
    according to the partial order specified by {orow} and {ocol}.

    If {orow} is positive, the row indices {arow,brow} must be increasing; If
    {orow} is negative, they must be decreasing. If {orow} is zero,
    the row indices are irrelevant relevant for the ordering. The
    {ocol} parameter refers to the column indices {acol,bcol}, in the same way.

    The absolute values of {orow} and {ocol} (which must be different)
    specify the order of the comparisons: a larger value means that the
    corresponding index is more important.

    The result is negative if the pair {arow,acol} precedes the pair
    {brow,bcol} in the partial order; it is positive if {arow,acol}
    follows {brow,bcol}; and is zero if the two pairs are equivalent
    in the partial order. */

void *spmat_alloc(spmat_count_t ents, size_t elsz);
  /* Allocates an area with space for {ents} entries of size {elsz}.
    Bombs out if there is no space for the request.
    If {ents == 0}, the result is {NULL}. */

void spmat_expand
  ( void **eP, 
    spmat_count_t *entsP, 
    spmat_pos_t pos, 
    size_t elsz
  );
  /* Makes sure that element {(*eP)[pos]} exists, reallocating and
    copying the array {**eP} if {pos >= (*entsP)}. If that happens,
    the new {(*entsP)} will be about twice as big as the old one.
    Assumes that each entry has size {elsz}. */

void spmat_trim
  ( void **eP, 
    spmat_count_t *entsP, 
    spmat_count_t ents, 
    size_t elsz
  );
  /* Makes sure that {(*entsP) == ents}, reallocating and copying
    the  array {**eP} if that is not the case. If {ents == 0}, will
    set {(*eP) == NULL}.  Assumes that each entry has size {elsz}. */

/* MAXIMUM SAFE VALUES */

#define spmat_MAX_SIZE (UINT32_MAX/2)
#define spmat_MAX_ROWS (spmat_MAX_SIZE)
#define spmat_MAX_COLS (spmat_MAX_SIZE)
  /* Max safe value for a row or column count ({spmat_size_t}). */

#define spmat_MAX_INDEX (spmat_MAX_SIZE-1)
  /* Max safe value for a valid row or column index ({spmat_index_t}). */

#define spmat_NO_INDEX (spmat_MAX_SIZE)
  /* A {spmat_index_t} value used to mean `no such column' or `no such row'. */

#define spmat_MAX_COUNT (UINT32_MAX/2)
#define spmat_MAX_ENTS (spmat_MAX_COUNT)
  /* Max safe value for a valid stored entry count ({spmat_count_t}). */

#define spmat_MAX_POS (spmat_MAX_COUNT-1)
  /* Max safe value for a valid stored entry index ({spmat_pos_t}). */

#define spmat_NO_POS (spmat_MAX_POS)
  /* A {spmat_pos_t} value used to mean `no such entry'. */

#include <spmat_def.h>

#endif
