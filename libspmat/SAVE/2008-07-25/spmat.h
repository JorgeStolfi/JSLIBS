//@{sparse_matrix.h}
#ifndef sparse_matrix_H
#define sparse_matrix_H
/* Generic sparse matrix types, and operations thereon */
/* Last edited on 2008-07-25 00:56:24 by stolfi */

/* 
  !!! change sort_entries so that (+1,+2) means row-by-row ? !!!
  !!! define {PREFIX##_extract_element} !!!

  !!! unify _add_{row,column} into {PREFIX##_add_elements} with {drow,dcol} bits ? !!!

  !!! {PREFIX##_add_row} should store {val[j]} into {M[row,col+j]} !!!
  !!! {PREFIX##_add_row} should accept {vals != M.cols} !!!
  !!! define also {PREFIX##_fill_row} !!!

  !!! {PREFIX##_add_column} should store {val[i]} into {M[row+i,col]} !!!
  !!! {PREFIX##_add_column} should accept {vals != M.rows} !!!
  !!! define also {PREFIX##_fill_column} !!!

  !!! define {PREFIX##_add_block} that adds a submatrix. !!! 

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
  
  This interface defines a /packed representation/ for sparse matrices,
  namely a list of the element values together with their row and
  column indices in the original (/unpacked/) representation. Any
  element not on this list is assumed to be trivial. The
  representation is handled through a /descriptor/, a record
  containing the nominal dimensions (row and column counts) of the
  matrix, the count of elements stored in the list, and a pointer to
  an array of (row,col,value) triplets.

  USING SPARSE MATRICES
  
  A new sparse matrix type, with elements of a specific type, 
  is declared with the {sparse_matrix_typedef}
  macro, and implemented with {sparse_matrix_impl} macro.  These
  macros also declare and implement the basic functions for handling
  such sparse matrices.  Additional functions are available through
  the interfaces  {sparse_matrix_io.h} and {sparse_matrix_linalg.h}.

  The incremental creation of sparse matrices is both easy and
  efficient with the tools provided here. For example, suppose that
  the type {dmat_t} was previously defined as a sparse matrix of
  {double} values. The following bit of code creates such a matrix,
  with 50 rows and 200 columns, and fills it one element at a time:
  
    ------------------------------------------------------------
    // Create a {50 × 200} sparse matrix {M}
    dmat_t M = dmat_new(50,200,400); // Guesing ~400 entries.
    
    // Store some elements into {M}, in arbitrary order:
    int count = 0;
    while(! finished(...))
      { int i = ...;
        int j = ...;
        double Mij = ...; 
        count = dmat_add_element(&M, count, i, j, Mij);
      }
    
    // Reclaim unused entries:
    dmat_trim(&M, count); 
    
    // Sort the entries by rows:
    dmat_sort_entries(&M, +1, 0);
    
    // Print the elements, one row per line:
    int r = -1; // Current row index.
    for (k = 0; k < M.ents; k++)
      { dmat_entry_t *ek = &(M.e[k]);
        if (ek->row != r) { r = ek->row; printf("\n [%d] :", r); }
        printf(" [%d]=%f", ek->col, ek->val);
      }
    printf("\n");
    ------------------------------------------------------------
    
  The function {dmat_add_element} above will automatically reallocate
  {M}'s entry list, as needed, to accomodate the new entries. Each
  expansion doubles the size of the entry list, so the total cost of
  multiple expansions is small and proportional to the number of
  elements added. The call {dmat_trim} truncates the list to the
  specified size, freeing the unused space.
*/

#include <stdlib.h>
#include <bool.h>
#include <ref.h>
#include <stdint.h>

/* DECLARING A NEW SPARSE MATRIX TYPE

  Here is how one would define the type {dmat_t} as a sparse matrix of
  doubles:
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <sparse_matrix.h>
    
    sparse_matrix_typedef(dmat_t, dmat, double);
    ------------------------------------------------------------
    
    ------------------------------------------------------------
    / * Implementations  * /

    #define dmat_trivial_elem (0.0)
    #define dmat_elem_is_trivial(X) ((X)==0.0)
    
    sparse_matrix_impl(dmat_t, dmat, double);
    ------------------------------------------------------------

  As another example, here is how one would define a type {graph_t}
  as being a sparse directed graph, represented by a boolean adjacency
  matrix where most of the entries are FALSE (meaning `no edge'):
  
    ------------------------------------------------------------
    / * Declarations * /
    #include <sparse_matrix.h>
    
    sparse_matrix_typedef(graph_t, graph, bool_t);
    ------------------------------------------------------------

    ------------------------------------------------------------
    / * Implementations  * /

    #define graph_trivial_elem (FALSE)
    #define graph_elem_is_trivial(X) ((X)==FALSE)
    
    sparse_matrix_impl(graph_t, graph, bool_t);
    ------------------------------------------------------------ */

#define sparse_matrix_typedef(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_TYPEDEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/* 
  This macro defines a data type called {MATRIX_TYPE}, which is a
  record that describes a sparse matrix with element values are of
  type {ELEM_TYPE}. The operations on those matrices will have
  names beginning with {PREFIX
  
  This macro will also define the type {{PREFIX}_entry_t} (a triplet
  with fields {.row,.col,.value}), and will generate prototype
  declarations for the following procedures:  
  
    ------------------------------------------------------------
    {PREFIX}_new
    {PREFIX}_expand
    {PREFIX}_trim
    {PREFIX}_make_desc
    {PREFIX}_copy
    {PREFIX}_add_element
    {PREFIX}_add_row
    {PREFIX}_extract_row
    {PREFIX}_add_column
    {PREFIX}_extract_column 
    {PREFIX}_add_diagonal
    {PREFIX}_fill_diagonal
    {PREFIX}_sort_entries
    {PREFIX}_traspose
    ------------------------------------------------------------
      
  These procedures are described below.
     
  The names {MATRIX_TYPE} and {PREFIX} can be chosen quite arbitrarily.
  However, it is recommented to use {{PREFIX}_t} as the {MATRIX_TYPE}.
  For example, 
  
    ------------------------------------------------------------
    sparse_matrix_typedef(mymatrix_t, mymatrix, myelem_t)
    sparse_matrix_typedef(cmat_t, cmat, complex)
    sparse_matrix_typedef(graph_t, graph, bool_t)
    ------------------------------------------------------------
 */
 
#define sparse_matrix_impl(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the implementations of the functions whose
  prototypes were declared with {sparse_matrix_typedef}.
  
  BEFORE calling this macro, the client must also declare the following 
  macros:
    
    ------------------------------------------------------------
    #define {PREFIX}_trival_elem (...)
      / * The trivial element value. * / 

    #define {PREFIX}_elem_is_trivial(X) (...)
      / * A boolean expression that yields TRUE iff (X) is a
        trivial element value. * / 
    ------------------------------------------------------------  
    
  The macro {{PREFIX}_elem_is_trivial(X)} is used by the sparse
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
  by the macros {sparse_matrix_typedef} and {sparse_matrix_impl}.  */

/* DATA TYPES */

/* ------------------------------------------------------------
typedef struct MATRIX_TYPE
  { uint32_t rows;
    uint32_t cols;
    uint32_t ents;
    PREFIX##_entry_t *e;
  } MATRIX_TYPE;
------------------------------------------------------------ */
  /* A variable {M} of this type is a descriptor for a sparse matrix
    of {ELEM_TYPE} values that has {M.rows} rows and {M.cols} columns
    in the unpacked form, but has only {M.ents} explicitly stored
    elements, namely the {.val} fields of the entries
    {M.e[0..M.ents-1]}.

    Every entry {t} in {M.e[]} must have {t.row < M->rows} and {t.col <
    M->cols}. The entries may be in any order, but there must not be two
    entries with the same {row,col} indices. The address {.e} may be
    NULL if (and only if) {M.ents == 0}.*/

/* ------------------------------------------------------------
typedef struct PREFIX##_entry_t
  { uint32_t row;
    uint32_t col;
    ELEM_TYPE val
  } PREFIX##_ENTRY_t;
------------------------------------------------------------ */
  /* A variable {e} of this type represents an element of the sparse
    matrix. It contains the indices {e.row} and {e.col} of the 
    element, and its value {e.val} (which is usually non-trivial). */

/* MATRIX ALLOCATION */

/* ------------------------------------------------------------
MATRIX_TYPE PREFIX##_new(uint32_t rows, uint32_t cols, uint32_t ents);
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
void PREFIX##_expand(MATRIX_TYPE *M, uint32_t pos);
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
void PREFIX##_trim(MATRIX_TYPE *M, uint32_t ents);
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
  ( uint32_t rows,
    uint32_t cols,
    ELEM_TYPE *e,
    uint32_t ents
  );
------------------------------------------------------------ */
  /* Assembles a sparse matrix descriptor from the row count {rows}, the
    column count {cols}, the entry count {ents}, and the address {e} of
    the first explicit entry. The client must make sure that the entries
    {e[0..ents-1]} actually exist. The procedures {{PREFIX}_expand} and
    {{PREFIX}_trim} can be applied to the resulting descriptor only if
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

/* ADDING AND EXTRACTING SINGLE ELEMENTS */

/* ------------------------------------------------------------
uint32_t PREFIX##_add_element
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t row,
    uint32_t col,
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

/* ADDING AND EXTRACTING ELEMENTS BY ROWS OR COLUMNS */

/* ------------------------------------------------------------
uint32_t PREFIX##_add_row
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t row,
    ELEM_TYPE val[],
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Stores {val[0..vals-1]} into {M->e}, starting at position
    {M->e[pos]}, as elements of row {row} of the matrix.

    More precisely, if {row >= M.rows}, then {M.rows} is set to {row+1};
    if {vals > M.cols}, then {M.cols} is set to {valss}. Then, each
    non-trivial element {val[j]} is stored into {M->e}, as the entry in
    row {row} and column {j} of the matrix.

    The vector {M->e} will be reallocated as needed (with
    {PREFIX##_expand}) to accomodate the new entries.

    The function returns the index of the entry of {M.e} immediately
    after the last entry that was assigned; or {pos} itself, if all
    elements of {val} were trivial. The client must make sure that {M}
    does not contain any entry with the same indices as the stored
    elements.*/

/* ------------------------------------------------------------
uint32_t PREFIX##_add_column
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t col,
    ELEM_TYPE val[],
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Stores {val[0..vals-1]} into {M->e}, starting at position
    {M->e[pos]}, as elements of column {col} of the matrix.

    More precisely, if {col >= M.cols}, then {M.cols} is set to {col+1};
    if {vals > M.rows}, then {M.rows} is set to {vals}. Then, each
    non-trivial element {val[i]} is stored into {M->e}, as the entry in
    row {i} and column {col} of the matrix.

    The vector {M->e} will be reallocated as needed (with
    {PREFIX##_expand}) to accomodate the new entries.

    The function returns the index of the entry of {M.e} immediately
    after the last entry that was assigned; or {pos} itself, if all
    elements of {val} were trivial. The client must make sure that {M}
    does not contain any entry with the same indices as the stored
    elements. */
    
/* ------------------------------------------------------------
uint32_t PREFIX##_extract_row
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t row,
    ELEM_TYPE val[],
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Expands row {row} of {M}, which starts at {M->e[pos]}, and stores
    it into the ordinary vector {val[0..vals-1]}.

    The function assumes that all entries of {M} with row index {row} are
    stored in consecutive positions of {M->e}, starting with entry
    {M->e[pos]}. Thus, if {M->e[pos].row != row}, the function assumes
    that all elements of that row are trivial.

    The size {vals} of {val} must be equal to {M->ncols}. The function
    returns the index of the entry in {M->e} immediately after the last
    entry that was actually extracted; or {pos} itself, if that row was
    completely trivial.*/

/* ------------------------------------------------------------
void PREFIX##_extract_column
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t col,
    ELEM_TYPE val[],
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Expands column {col} of {M}, which starts at {M->e[pos]}, and stores
    it into the ordinary vector {val[0..vals-1]}.

    The function assumes that all entries of {M} with column index {col}
    are stored in consecutive positions of {M->e}, starting with entry
    {M->e[pos]}. Thus, if {M->e[pos].col != col}, the function assumes
    that all elements of that column are trivial.

    The size {vals} of {val} must be equal to {M->nrows}. The function
    returns the index of the entry in {M->e} immediately after the last
    entry that was actually extracted; or {pos} itself, if that column was
    completely trivial.

    Here is a typical example of usage of these procedures:

      ------------------------------------------------------------
      dmat_sort_entries(&M, +1, 0);  // Sort entries of {M} by row.
      dmat_t R = dmat_new(M->rows,M->cols,0);  // Matrix for result.
      double v[M.cols];
      int pM = 0, pR = 0; // Scan the entries of {M} and {R}.
      int i;              // Scans the rows of {M} and {R}.
      for (i = 0; i < M.rows; i++)
        { pM = dmat_extract_row(&M, pM, i, v);
          for (j = 0; j < M.cols; j++) { modify(v[j]); }
          pR = dmat_add_row(&R, pR, i, v);
        }
      dmat-trim(&R, pR);
      ------------------------------------------------------------
  */

/* ADDING ELEMENTS ALONG DIAGONALS*/

/* ------------------------------------------------------------
uint32_t PREFIX##_add_diagonal
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t row,
    uint32_t col,
    ELEM_TYPE val[],
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Stores {val[0..vals-1]} into {M->e}, starting at position
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
    elements of {val} were trivial. The client must make sure that {M}
    will not end up with two entries with the same indices. Note that
    this may happen if {vals} exceeds the least common multiple of
    {M.rows} and {M.cols}. */

/* ------------------------------------------------------------
uint32_t PREFIX##_fill_diagonal
  ( MATRIX_TYPE *M,
    uint32_t pos,
    uint32_t row,
    uint32_t col,
    ELEM_TYPE val,
    uint32_t vals
  );
------------------------------------------------------------ */
  /* Similar to {}, but uses the single value {val},
    replicated {vals} times. */

/* SORTING ENTRIES BY INDICES*/

/* ------------------------------------------------------------
void PREFIX##_sort_entries(MATRIX_TYPE *M, int orow, int ocol);
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

/* TRANSPOSITION */

/* ------------------------------------------------------------
void PREFIX##_transpose(MATRIX_TYPE *M);
------------------------------------------------------------ */
  /* Exchanges the {.row} and {.col} indices of all
    entries of {M}, and also exchanges {M->rows} with {M->cols}.

    It does not change the order of the entries within {M.e}; so, if {M}
    was sorted by rows before the call, it will be sorted by columns
    after it.  */

#include <sparse_matrix_def.h>

#endif
//@{sparse_matrix_io.h}
#ifndef sparse_matrix_io_H
#define sparse_matrix_io_H
/* Reading and writing generic sparse matrices */
/* Last edited on 2008-07-24 18:50:40 by stolfi */

/*
  This interface provides basic I/O operations on generic sparse
  matrices.

  For example, if {dmat_t} was defined as a sparse matrix of {double}
  values, one can write

    ------------------------------------------------------------
    //Read a sparse matrix {M} from {stdin}:
    dmat_t M = dmat_new(0,0,0);
    dmat_read(stdin, &M);

    //Write the result to {stdout}:
    dmat_write(stdout, &M);
    ------------------------------------------------------------ */
  /* The I/O functions are made available by the macros
  {sparse_matrix_io_def} and {sparse_matrix_io_impl} (see below).
  For example,

    ------------------------------------------------------------
    / * Declarations * /
    #include <sparse_matrix_io.h>

    sparse_matrix_typedef(dmat_t, dmat, double);

    sparse_matrix_io_def(dmat_t, dmat, double);
    ------------------------------------------------------------ */
  /*   ------------------------------------------------------------
    / * Implementations  * /
    #include <sparse_matrix_io.h>

    sparse_matrix_typedef(dmat_t, dmat, double);

    void dmat_elem_read(FILE *rd, double *valP)
      { fscanf(rd, "%f", valP); }

    void dmat_elem_write(FILE *wr, double *valP)
      { fprintf(wr, "%6.3f", *valP); }

    sparse_matrix_io_impl(dmat_t, dmat, double);
    ------------------------------------------------------------ */

#define sparse_matrix_io_def(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_IO_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the prototype declarations of
  four procedures, explained later:

     ------------------------------------------------------------
     {PREFIX}_write
     {PREFIX}_read
     {PREFIX}_elem_write
     {PREFIX}_elem_read
     ------------------------------------------------------------ */

#define sparse_matrix_io_impl(MATRIX_TYPE,PREFIX,ELEMT_TYPE) \
  sparse_matrix_IO_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  These macros expand into implementations of the procedures

     ------------------------------------------------------------
     {PREFIX}_write
     {PREFIX}_read
     ------------------------------------------------------------
  
  BEFORE calling this macro, the client must provide adequate
  implementations for the functions {{PREFIX}_elem_write} and
  {{PREFIX}_elem_read}, as explained below. */

/* ======================================================================

  The remainder of this interface describes the functions provided
  by the macros {sparse_matrix_io_def} and {sparse_matrix_io_impl}.  */

/* MATRIX I/O */

/* ------------------------------------------------------------
void PREFIX##_write(FILE *wr, MATRIX_TYPE *M);
------------------------------------------------------------ */
  /* Writes the sparse matrix {M} to the output stream {wr}.

    The output contains a distinctive header line, then three lines with
    {M}'s dimensional attributes ({M->rows}, {M->cols}, {M->ents}), then
    the entries {M->e[0..M->ents-1]}, one per line (with row, column,
    and value); and finally a distinctive `footer line.

    Uses the client-defined function {PREFIX##-elem_write} to write each
    value; which will be preceded by a space, and followed by an end-of-line.*/

/* ------------------------------------------------------------
void PREFIX##_read(FILE *rd, MATRIX_TYPE *M);
------------------------------------------------------------ */
  /* Reads a sparse matrix from file {rd} and stores it into {M}.

    The fields {M.rows} and {M.cols} will be reset as specified in the
    file. The vector {M.e} will be re-alocated if necessary to contain
    exactly the entries read from the file.

    Uses the client-defined function PREFIX##_elem_read} to parse
    each value. */

/* ELEMENT I/O FUNCTIONS

  The macro {sparse_matrix_io_def} also declares the following
  procedures:*/

/* ------------------------------------------------------------
void PREFIX##_elem_write(FILE *wr, ELEM_TYPE *valP);
------------------------------------------------------------ */
  /* Writes to the output stream {wr} the value {*valP}.

    This procedure must be implemented by the client. It is often
    advisable to use a human- and machine-readable plain text
    representation that preserves all the relevant information contained
    in {*valP}. For example, a prudent choice for {double} values would
    be {fprintf(wr, "%24.16e", *valP)}.*/

/* ------------------------------------------------------------
void PREFIX##_elem_read(FILE *rd, ELEM_TYPE *valP);
------------------------------------------------------------ */
  /* Parses an {ELEM_TYPE} value from the input stream {rd}, and stores
    it into the variable {val}.

    This function must be implemented by
    the client. It is often advisable to define it so that it can parse
    any output that may be produced by {{PREFIX}_elem_write}. */

#include <sparse_matrix_io_def.h>

#endif
//@{sparse_matrix_linalg.h}
#ifndef sparse_matrix_linalg_H
#define sparse_matrix_linalg_H
/* Linear algebra operations on generic sparse matrices. */
/* Last edited on 2008-07-24 23:47:17 by stolfi */

#include <stdlib.h>
#include <bool.h>
#include <ref.h>
#include <stdint.h>
#include <sparse_matrix.h>

/* LINEAR ALGEBRA OPERATIONS

  This interface provides some basic linear algebra functions
  on generic sparse matrices.  For example, if {dmat_t}
  was defined as a sparse matrix with {double} elements,
  one could write:

    ------------------------------------------------------------
    //Create an {n × n} identity matrix {I}:
    dmat_t I = dmat_new(n,n,n);
    dmat_identity(&I,n);

    //Set {R = M - 2*I}:
    dmat_t R = dmat_new(n,n,0);
    dmat_mix(1.0, &M, -2.0, &I, &R);

    //Set {S = R*R}:
    dmat_t S = dmat_new(n,n,0);
    dmat_mul(&R, &R, &S);

    //Write the result to {stdout}:
    dmat_write(stdout, &S);
    ------------------------------------------------------------ */
  /* To obtain these linear algebra functions, one should use

    ------------------------------------------------------------
    / * Declarations * /
    #include <sparse_matrix_linalg.h>

    sparse_matrix_typedef(dmat_t, dmat, double);

    sparse_matrix_linalg_def(dmat_t, dmat, double);
    ------------------------------------------------------------ */
  /*   ------------------------------------------------------------
    / * Implementations  * /
    #include <sparse_matrix_linalg.h>

    sparse_matrix_impl(dmat_t, dmat, double);

    #define dmat_trivial_elem() (0.0)
    #define dmat_elem_is_trivial(X) ((X)==0.0)

    #define dmat_elem_zero (0.0)
    #define dmat_elem_one (1.0)
    #define dmat_elem_add(X,Y) ((X)+(Y))
    #define dmat_elem_mul(X,Y) ((X)*(Y))

    sparse_matrix_impl(dmat_t, dmat, double);
    sparse_matrix_linalg_impl(dmat_t, dmat, double);
    ------------------------------------------------------------

 */

#define sparse_matrix_linalg_def(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_LINALG_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into prototype declarations
  of the basic linear algebra procedures:

    ------------------------------------------------------------
    {PREFIX}_identity
    {PREFIX}_mix
    {PREFIX}_mul
    {PREFIX}_map_col
    {PREFIX}_map_row
    ------------------------------------------------------------
  
  They are approriate when the addition and multiplication of
  {ELEM_TYPE} values are defined. */

#define sparse_matrix_linalg_impl(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_LINALG_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the implementations of the linear algebra
  procedures declared by {sparse_matrix_linalg_def}.

  BEFORE calling this macro, the client must declare the following
  macros or procedures:

    ------------------------------------------------------------
    #define {PREFIX}_elem_add(X,Y) (...)
       / * The addition of two {ELEM_TYPE} values {X,Y}. * /

    #define {PREFIX}_elem_mul(X,Y) (...)
       / * The multiplication of two {ELEM_TYPE} values {X,Y}. * /

    #define {PREFIX}_elem_one (...)
      / * The {ELEM_TYPE} value that is the multiplicative unit. * /
    ------------------------------------------------------------
  
  This interface assumes that {{PREFIX}_trivial_elem} (the element
  value which is not stored explicitly in the matrix) is the
  zero value for {{PREFIX}_elem_add}.

/* ======================================================================

  The remainder of this interface describes the functions provided by
  the macros {sparse_matrix_linalg_def} and {sparse_matrix_linalg_impl}. */

/* IDENTITY MATRIX*/

/* ------------------------------------------------------------
void PREFIX##_identity(MATRIX_TYPE *M, uint32_t size);
------------------------------------------------------------ */
  /* Sets the matrix {M} to the identity matrix with
    {size} rows and {size} columns.*/

/* ------------------------------------------------------------
void PREFIX##_mix
  ( ELEM_TYPE a,
    MATRIX_TYPE *A,
    ELEM_TYPE b,
    MATRIX_TYPE *B,
    MATRIX_TYPE *C
  );
------------------------------------------------------------ */
  /* Stores into {C} the linear combination {a*A + b*B}.

    The entries of {A} MUST BE SORTED by increasing row index, then by
    increasing column index within each row (as produced by
    {{PREFIX}_sort_entries(&A, +2, +1)}); and ditto for {B}. The entries
    of {C} will be in the same order. The matrices {A} and {B} may be
    the same, but the matrix {C} must be storage-disjoint from both.

    The matrices {A} and {B} must also have the same dimensions
    {A.rows==B.rows} and {A.cols==B.cols}, which will be assigned to
    {C.rows} and {C.cols}. The vector {C.e} will be re-alocated, if
    necessary, to hold exactly the non-trivial entries of the result.*/

/* ------------------------------------------------------------
void PREFIX##_mul(MATRIX_TYPE *A, MATRIX_TYPE *B, MATRIX_TYPE *C);
------------------------------------------------------------ */
  /* Stores into {C} the matrix product {A*B}.

    The entries of {A} (only) MUST BE SORTED by increasing row index
    (as obtained with {{PREFIX}_sort_entries(&A, +1, 0)}.
    The entries of {C} will be sorted by increasing row index,
    then by increasing column index within each row. The matrices {A} and {B}
    may be the same, but the matrix {C} must be
    storage-dsjoint from both.

    The procedure requires {A.cols == B.rows}, and will set
    {R.rows=A.rows} and {R.cols=B.cols}. The vector {R.e} will be
    re-alocated if necessary to hold exactly the non-trivial entries
    of the result. */

/* MATRIX-VECTOR MULTIPLICATION*/

/* ------------------------------------------------------------
void PREFIX##_map_col
  ( MATRIX_TYPE *M,
    ELEM_TYPE a[],
    uint32_t na,
    ELEM_TYPE b[],
    uint32_t nb
  );
------------------------------------------------------------ */
  /* This macro declares the function {{PREFIX}_map_col}. The call
    {{PREFIX}_map_col(&M, a, na, b, nb)} will compute
    the product of the matrix {M} by the column vector
    {a[0..na-1]} and store the result into {b[0..nb-1]}.

    The procedure requires {na == M.cols} and {nb == M.rows}. The
    vectors {a} and {b} must be storage-disjoint. The entries of {M}
    need not be sorted. */*/

/* ------------------------------------------------------------
void PREFIX##_map_row
  ( ELEM_TYPE a[],
    uint32_t na,
    MATRIX_TYPE *M,
    ELEM_TYPE b[],
    uint32_t nb
  );
------------------------------------------------------------ */
  /* This macro declares the function {{PREFIX}_map_row}. The call
    {{PREFIX}_map_row(a, na, &M, b, nb)} will compute
    the product of the row vector
    {a[0..na-1]} by the matrix {M}, and store the result into {b[0..nb-1]}.

    The procedure requires {na == M.rows} and {nb == M.cols}. The
    vectors {a} and {b} must be storage-disjoint. The entries of {M}
    need not be sorted. */

#include <sparse_matrix_linalg_def.h>

#endif
//@{sparse_matrix_def.h}
#ifndef sparse_matrix_def_H
#define sparse_matrix_def_H
/* The ugly entrails of {sparse_matrix.h}. */
/* Last edited on 2008-07-23 19:21:41 by stolfi */

/* DEFINITION MACROS */

#define sparse_matrix_TYPEDEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_DEFINE_ENTRY_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DEFINE_MATRIX_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_NEW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_EXPAND(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_TRIM(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_MAKE_DESC(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_COPY(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_IS_TRIVIAL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_TRIVIAL_ELEM(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_ADD_ELEMENT(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_ADD_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_ADD_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_ADD_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_FILL_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_EXTRACT_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_EXTRACT_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_SORT_ENTRIES(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_TRANSPOSE(MATRIX_TYPE,PREFIX,ELEM_TYPE)

#define sparse_matrix_DEFINE_ENTRY_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct PREFIX##_entry_t \
    { uint32_t row; \
      uint32_t col; \
      ELEM_TYPE val \
    } PREFIX##_ENTRY_t

#define sparse_matrix_DEFINE_SPARSE_MATRIX_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct MATRIX_TYPE \
    { uint32_t rows; \
      uint32_t cols; \
      uint32_t ents; \
      PREFIX##_entry_t *e; \
    } MATRIX_TYPE

#define sparse_matrix_DECLARE_NEW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_new(uint32_t rows, uint32_t cols, uint32_t ents)

#define sparse_matrix_DECLARE_EXPAND(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(MATRIX_TYPE *M, uint32_t pos)

#define sparse_matrix_DECLARE_TRIM(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(MATRIX_TYPE *M, uint32_t ents)

#define sparse_matrix_DECLARE_MAKE_DESC(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_make_desc(uint32_t rows, uint32_t cols, ELEM_TYPE *e, uint32_t ents)

#define sparse_matrix_DECLARE_COPY(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_copy(MATRIX_TYPE *M, MATRIX_TYPE *N)

#define sparse_matrix_DECLARE_ADD_ELEMENT(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_element(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val)

#define sparse_matrix_DECLARE_ADD_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_row(MATRIX_TYPE *M, uint32_t pos, uint32_t row, ELEM_TYPE val[], uint32_t vals)

#define sparse_matrix_DECLARE_ADD_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_column(MATRIX_TYPE *M, uint32_t pos, uint32_t col, ELEM_TYPE val[], uint32_t vals)

#define sparse_matrix_DECLARE_ADD_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_diagonal(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val[], uint32_t vals)

#define sparse_matrix_DECLARE_FILL_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_fill_diagonal(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val, uint32_t vals)

#define sparse_matrix_DECLARE_EXTRACT_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_extract_row(MATRIX_TYPE *M, uint32_t pos, uint32_t row, ELEM_TYPE val[], uint32_t vals)

#define sparse_matrix_DECLARE_EXTRACT_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_extract_column(MATRIX_TYPE *M, uint32_t pos, uint32_t col, ELEM_TYPE val[], uint32_t vals)

#define sparse_matrix_DECLARE_SORT_ENTRIES(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_sort_entries(MATRIX_TYPE *M, int orow, int ocol)

#define sparse_matrix_DECLARE_TRANSPOSE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_transpose(MATRIX_TYPE *M)

/* IMPLEMENTATION MACROS */

#define sparse_matrix_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_IMPLEMENT_NEW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_EXPAND(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_TRIM(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_MAKE_DESC(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_COPY(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_ADD_ELEMENT(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_ADD_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_ADD_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_ADD_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_FILL_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_EXTRACT_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_EXTRACT_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_SORT_ENTRIES(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_TRANSPOSE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_bOgUs /* To eat the semicolon. */

#define sparse_matrix_IMPLEMENT_NEW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_new(uint32_t rows, uint32_t cols, uint32_t ents) \
    { sparse_matrix_t M = sparse_matrix_alloc(rows, cols, ents, sizeof(PREFIX##_entry_t)); \
      return (MATRIX_TYPE){M.rows, M.cols, M.ents, (ELEM_TYPE *)M.e}; \
    } \

#define sparse_matrix_IMPLEMENT_EXPAND(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(MATRIX_TYPE *M, int pos) \
    { if (pos >= M->ents) \
        { sparse_matrix_expand((sparse_matrix_t *)M, pos, sizeof(PREFIX##_entry_t)); } \
    } \

#define sparse_matrix_IMPLEMENT_TRIM(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(MATRIX_TYPE *M, uint32_t ents) \
    { sparse_matrix_trim((sparse_matrix_t *)M, ents, sizeof(PREFIX##_entry_t)); } \

#define sparse_matrix_IMPLEMENT_MAKE_DESC(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_make_desc(uint32_t rows, uint32_t cols, PREFIX##_entry_t *e, uint32_t ents) \
    { return (MATRIX_TYPE){rows, cols, ents, e}; } \

#define sparse_matrix_IMPLEMENT_COPY(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_copy(PREFIX *M, PREFIX *N) \
    { N->rows = M->rows; N->cols = M->cols; \
      PREFIX##_trim(N, M->ents); \
      int k; \
      for (k = 0; k < M->ents; k++) { N->e[k] = M->e[k]; } \
    }

#define sparse_matrix_IMPLEMENT_ADD_ELEMENT(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_element(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val) \
    { if PREFIX##_elem_is_trivial(val) \
        { return pos; } \
      else \
        { PREFIX##_expand(M, pos); \
          M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = col, .val = val }; \
          return pos+1; \
        } \
    }

#define sparse_matrix_IMPLEMENT_ADD_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_row(MATRIX_TYPE *M, uint32_t pos, uint32_t row, ELEM_TYPE val[], uint32_t vals) \
    { if (row >= M->rows) { M->rows = row + 1; } \
      if (vals > M->cols) { M->cols = vals; } \
      int k; \
      for (k = 0; k < vals; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = k, .val = (*vk) }; \
              pos++; \
            } \
        } \
      return pos; \
    } \

#define sparse_matrix_IMPLEMENT_ADD_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
   uint32_t PREFIX##_add_column(MATRIX_TYPE *M, uint32_t pos, uint32_t col, ELEM_TYPE val[], uint32_t vals) \
    { if (vals > M->rows) { M->rows = vals; } \
      if (col >= M->cols) { M->cols = col + 1; } \
      int k; \
      for (k = 0; k < vals; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = k, .col = col, .val = (*vk) }; \
              pos++; \
            } \
        } \
      return pos; \
    }

#define sparse_matrix_IMPLEMENT_ADD_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_add_diagonal(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val[], uint32_t vals) \
    { row = row % M.rows; \
      col = col % M.cols; \
      int k; \
      for (k = 0; k < vals; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = k, .col = col, .val = (*vk) }; \
              pos++; \
            } \
          row++; if (row >= M.rows) { row = 0; } \
          col++; if (col >= M.cols) { col = 0; } \
        } \
      return pos; \
     } \

#define sparse_matrix_IMPLEMENT_FILL_DIAGONAL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_fill_diagonal(MATRIX_TYPE *M, uint32_t pos, uint32_t row, uint32_t col, ELEM_TYPE val, uint32_t vals) \
    { if (PREFIX##_is_trivial(val)) { return pos; }
      row = row % M.rows; \
      col = col % M.cols; \
      int k; \
      for (k = 0; k < vals; k++) \
        { PREFIX##_expand(M, pos); \
          M->e[pos] = (PREFIX##_entry_t){ .row = k, .col = col, .val = (*vk) }; \
          pos++; \
          row++; if (row >= M.rows) { row = 0; } \
          col++; if (col >= M.cols) { col = 0; } \
        } \
      return pos; \
    } \

#define sparse_matrix_IMPLEMENT_EXTRACT_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_extract_row(MATRIX_TYPE *M, uint32_t pos, uint32_t row, ELEM_TYPE val[], uint32_t vals) \
    { demand(vals == M.cols, "incompatible vector size"); \
      int k; \
      for (k = 0; k < vals; k++) { val[k] = PREFIX##_trivial_elem; } \
      while (pos < M->ents) \
        { PREFIX##_entry_type_t *eP = &(M->e[pos]); \
          if (eP->row != row) { break; } \
          assert(eP->col < M.cols); \
          val[eP->col] = eP->val; \
          pos++; \
        };
      return pos; \
    } \

#define sparse_matrix_IMPLEMENT_EXTRACT_COLUMN(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  uint32_t PREFIX##_EXTRACT_COLUMN(MATRIX_TYPE *M, uint32_t pos, uint32_t col, ELEM_TYPE val[], uint32_t vals) \
    { demand(vals == M.rows, "incompatible vector size"); \
      int k; \
      for (k = 0; k < vals; k++) { val[k] = PREFIX##_trivial_elem; } \
      while (pos < M->ents) \
        { PREFIX##_entry_type_t *eP = &(M->e[pos]); \
          if (eP->col != col) { break; } \
          assert(eP->row < M.rows); \
          val[eP->row] = eP->val; \
          pos++; \
        };
      return pos; \
    } \

#define sparse_matrix_IMPLEMENT_SORT_ENTRIES(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_sort_entries(MATRIX_TYPE *M, int orow, int ocol) \
    { if ((trow == 0) && (tcol == 0)) { return; } \
      auto int cmp(const void *avP, const void *bvP); \
      int cmp(const void *avP, const void *bvP) \
        { PREFIX##_entry_t *aP = avP; \
          PREFIX##_entry_t *bP = bvP; \
          return sparse_matrix_compare_indices(aP->row, aP->col, bP->row, bP->col, trow, tcol); \
        } \
      qsort(M->e, M->ents, sizeof(PREFIX##_entry_t), &cmp); \
    }

#define sparse_matrix_IMPLEMENT_TRANSPOSE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_transpose(MATRIX_TYPE *M) \
    { uint32_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_type_t *eP = &(M->e[pos]); \
          uint32_t t = eP->row; eP->row = eP->col; eP->col = t; \
        } \
    }

#endif
//@{sparse_matrix_io_def.h}
#ifndef sparse_matrix_io_def_H
#define sparse_matrix_io_def_H
/* The ugly entrails of {sparse_matrix_io.h}. */
/* Last edited on 2008-07-23 19:21:41 by stolfi */

/* DEFINITION MACROS */

#define sparse_matrix_IO_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_DECLARE_WRITE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_READ(MATRIX_TYPE,PREFIX,ELEM_TYPE)

#define sparse_matrix_DECLARE_ELEM_WRITE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_write(FILE *wr, ELEM_TYPE *valP)

#define sparse_matrix_DECLARE_ELEM_READ(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_read(FILE *rd, ELEM_TYPE *valP)

#define sparse_matrix_DECLARE_WRITE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_write(FILE *wr, MATRIX_TYPE *M)

#define sparse_matrix_DECLARE_READ(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_read(FILE *rd, MATRIX_TYPE *M)

/* IMPLEMENTATION MACROS */

#define sparse_matrix_IO_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_IMPLEMENT_WRITE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_READ(MATRIX_TYPE,PREFIX,ELEM_TYPE)

#define sparse_matrix_IMPLEMENT_WRITE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_write(FILE *wr, MATRIX_TYPE *M) \
    { sparse_matrix_write_header(wr, "##MATRIX_TYPE##", NULL, M->rows, M->cols, M->ents); \
      unit32_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_type_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          fprintf(wr, "%7d %7d ", eP->row, eP->col); \
          PREFIX##_write_elem(wr, &(eP->val)); \
          fprintf(wr, "\n"); \
        }\
      sparse_matrix_write_footer(wr, "##MATRIX_TYPE##"); \
    }

#define sparse_matrix_IMPLEMENT_READ(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_read(FILE *rd, MATRIX_TYPE *M) \
    { uint32_t ents; \
      sparse_matrix_read_header(wr, "##MATRIX_TYPE##", NULL, &(M->rows), &(M->cols), &ents); \
      PREFIX##_trim(M, ents); \
      unit32_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_type_t *eP = &(M->e[pos]); \
          eP->row = fget_uint32(rd, 10); demand(eP->row < M->rows, "invalid row index"); \
          eP->row = fget_uint32(rd, 10); demand(eP->col < M->cols, "invalid col index"); \
          fget_skip_spaces(rd); \
          PREFIX##_read_elem(rd, &(eP->val)); \
          fget_eol(rd); \
        }\
      sparse_matrix_read_footer(rd, "##MATRIX_TYPE##"); \
    }

#endif
//@{sparse_matrix_linalg_def.h}
#ifndef sparse_matrix_linalg_def_H
#define sparse_matrix_linalg_def_H
/* The ugly entrails of {sparse_matrix_linalg.h}. */
/* Last edited on 2008-07-23 19:21:41 by stolfi */

#define sparse_matrix_LINALG_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_DECLARE_IDENTITY(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_MIX(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_MUL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_MAP_COL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_DECLARE_MAP_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE)
  
#define sparse_matrix_DECLARE_IDENTITY(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_identity(MATRIX_TYPE *M, uint32_t size)

#define sparse_matrix_DECLARE_MIX(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mix(ELEM_TYPE a, MATRIX_TYPE *A, ELEM_TYPE b, MATRIX_TYPE *B, MATRIX_TYPE *C)

#define sparse_matrix_DECLARE_MUL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mul(MATRIX_TYPE *A, MATRIX_TYPE *B, MATRIX_TYPE *C)

#define sparse_matrix_DECLARE_MAP_COL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_col(MATRIX_TYPE *M, ELEM_TYPE a[], uint32_t na, ELEM_TYPE b[], uint32_t nb)

#define sparse_matrix_DECLARE_MAP_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_row(ELEM_TYPE a[], uint32_t na, MATRIX_TYPE *M, ELEM_TYPE b[], uint32_t nb)

/* IMPLEMENTATION MACROS */

#define sparse_matrix_LINALG_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  sparse_matrix_IMPLEMENT_IDENTITY(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_MIX(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_MUL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_MAP_COL(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  sparse_matrix_IMPLEMENT_MAP_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE)

#define sparse_matrix_IMPLEMENT_IDENTITY(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_identity(MATRIX_TYPE *M, uint32_t size) \
    { M->rows = M->cols = size; \
      PREFIX##_trim(M, size); \
      uint32_t pos = PREFIX##_fill_diagonal(M, 0, 0,0, PREFIX##_elem_one, size); \
      assert(pos == size); \
      assert(M->ents == size);  \
    }

#define sparse_matrix_IMPLEMENT_MIX(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mix(ELEM_TYPE a, MATRIX_TYPE *A, ELEM_TYPE b, MATRIX_TYPE *B, MATRIX_TYPE *C) \
    { demand(A->rows == B->rows, "incompatible row counts"); C->rows = A->rows; \
      demand(A->cols == B->cols, "incompatible column counts"); C->cols = A->cols; } \
      uint32_t posA = 0, posB = 0, posC = 0; \
      PREFIX##_entry_t *eCprev = NULL; /* Last entry stored into {C}. */ \
      while ((posA < A->ents) || (posB < B->ents)) \
        { PREFIX##_entry_t *aP = (posA < A->ents ? &(A->e[posA]) : NULL); \
          PREFIX##_entry_t *bP = (posB < B->ents ? &(B->e[posB]) : NULL); \
          int cmp; /* Which entry comes first? */ \
          if (aP == NULL) \
            { cmp = +1; } \
          else if (bP == NULL) \
            { cmp = -1; } \
          else \
            { cmp = sparse_matrix_compare_indices(aP->row, aP->col, bP->row, bP->col, +2, +1); } \
          uint32_t row, col; \
          ELEM_TYPE val; \
          if (cmp < 0) \
            { row = aP->row; col = aP->col; val = PREFIX##_elem_mul(a, aP->val); posA++; } \
          else if (cmp > 0) \
            { row = bP->row; col = bP->col; val = PREFIX##_elem_mul(b, bP->val); posB++; } \
          else \
            { assert(aP->row == bP->row); \
              assert(aP->col == bP->col); \
              ELEM_TYPE va = PREFIX##_elem_mul(a, aP->val); \
              ELEM_TYPE vb = PREFIX##_elem_mul(b, bP->val); \
              row = aP->row; col = aP->col; val = PREFIX##_elem_add(va, vb); posA++; posB++; \
            } \
          if (eCprev != NULL) \
            { demand((row != eCprev->row) || (col != eCprev->col), "duplicated entries"); \
              demand((row > eCprev->row) || ((row == eCprev->row) && (col > eCprev->col)), "entries out of order"); \
            } \
          if (! PREFIX##_elem_is_trivial(val)) \
            { PREFIX##_expand(C, posC); \
              PREFIX##_entry_t *eC = &(C->e[posC]); \
              eC->row = row; eC->col = col; eC->val = val; \
              eCprev = eC; \
              posC++; \
            } \
        } \
      PREFIX##_trim(C, posC); \
    }

#define sparse_matrix_IMPLEMENT_MUL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mul(MATRIX_TYPE *A, MATRIX_TYPE *B, MATRIX_TYPE *C) \
    { demand(A->cols == B->rows, "incompatible row/column counts"); \
      C->rows = A->rows; C->cols = B->cols; \
      ELEM_TYPE *av = malloc(A->cols * sizeof(ELEM_TYPE)); /* A row of {A} */ \
      ELEM_TYPE *cv = malloc(B->cols * sizeof(ELEM_TYPE)); /* A row of {B} */ \
      uint32_t posA = 0; /* Start of next row of {A}. */ \
      uint32_t posC = 0; /* Next free eentry in {C}. */ \
      while (posA < A->ents) \
        { uint32_t row = A->e[posA].row; \
          posA = PREFIX##_extract_row(A, posA, row, av, A->cols); \
          PREFIX##_map_row(av, A->cols, B, cv, B->cols); \
          posC = PREFIX##_add_row(C, posC, row, cv, B->cols); \
        } \
      PREFIX##_trim(C, posC); \
      free(av); free(cv); \
    }

#define sparse_matrix_IMPLEMENT_MAP_COL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_col(MATRIX_TYPE *M, ELEM_TYPE a[], uint32_t na, ELEM_TYPE b[], uint32_t nb) \
    { demand(M->cols == na, "incompatible column counts"); \
      demand(M->rows == nb, "incompatible row counts"); \
      uint32_t i; \
      for (i = 0; i < nb; i++) { b[i] = PREFIX##_elem_zero; } \
      uint32_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          ELEM_TYPE *bP = &(b[eP->row]); \
          ELEM_TYPE p = PREFIX##_elem_mul(eP->val, a[eP->col]); \
          ELEM_TYPE s = PREFIX##_elem_add(p, *bP); \
          (*bP) = s; \
        } \
    }

#define sparse_matrix_IMPLEMENT_MAP_ROW(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_row(ELEM_TYPE a[], uint32_t na, MATRIX_TYPE *M, ELEM_TYPE b[], uint32_t nb) \
    { demand(M->rows == na, "incompatible row counts"); \
      demand(M->cols == nb, "incompatible column counts"); \
      uint32_t j; \
      for (j = 0; j < nb; j++) { b[j] = PREFIX##_elem_zero; } \
      uint32_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          ELEM_TYPE *bP = &(b[eP->col]); \
          ELEM_TYPE p = PREFIX##_elem_mul(a[eP->row], eP->val); \
          ELEM_TYPE s = PREFIX##_elem_add(p, *bP); \
          (*bP) = s; \
        } \
    }

/* GENERIC SPARSE MATRIX DESCRIPTORS */

typedef struct sparse_matrix_t /* General sparse matrix descriptor. */
  { uint32_t rows,   /* Number of rows. */
    uint32_t cols,   /* Number of columns. */
    uint32_t ents;   /* Number of explicit entries. */
    void *e;      /* Pointer to the first explicit entry. */
  } sparse_matrix_t;
  /* An {sparse_matrix_t} {M} describes a uni-dimensional array of {M.ents} elements
    of the same size and type, appendd in consecutive memory locations,
    starting at address {M.e}. */

/* GENERIC SPARSE MATRIX FUNCTIONS

  The functions in this section may be called directly by clients,
  but their main purpose is to implement the
  functions {my_sparse_matrix_new}, {my_sparse_matrix_expand}, and {my_sparse_matrix_trim}
  defined by {sparse_matrix_typedef} and {sparse_matrix_impl}. */

sparse_matrix_t sparse_matrix_alloc(uint32_t rows, uint32_t cols, uint32_t ents, size_t elsz);
  /* Allocates a new sparse matrix with {ents} elements of size {elsz}.
    Bombs out if there is no space for the request.
    If {ents == 0}, the result has {e == NULL}. */

void sparse_matrix_expand(sparse_matrix_t *M, uint32_t pos, size_t elsz);
  /* Makes sure that element {M->e[pos]} exists, reallocating and
    copying the array {*M->e} if {pos >= M->ents}. If that happens,
    the new {M->ents} will be about twice as big as the old one. */

void sparse_matrix_trim(sparse_matrix_t *M, uint32_t ents, size_t elsz);
  /* Makes sure that {M.ents == ents}, reallocating and copying
    the  array {*M->e} if {M.ents != ents}. If {ents == 0}, the
    result will have {e == NULL}. */

int sparse_matrix_compare_indices(uint32_t arow, uint32_t acol, uint32_t brow, uint32_t bcol, int orow, int ocol);
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

#endif

//@{sparse_matrix.c}
/* See sparse_matrix.h */
/* Last edited on 2008-03-29 14:33:12 by stolfi */

#include <sparse_matrix.h>
#include <stdlib.h>
#include <stdint.h>
#include <affirm.h>

void *sparse_matrix_alloc(uint32_t ents, size_t esz)
  { void *e = (ents == 0 ? NULL : malloc(ents*esz));
    affirm((ents == 0) || (e != NULL), "out of mem");
    return e;
  }

void sparse_matrix_expand(uint32_t *entsP, void **eP, uint32_t index, size_t esz)
  { if (index >= (*entsP))
      { int ents = (*entsP) + index + 1;
        if ((*entsP) == 0) { affirm((*eP) == NULL, "bad elem pointer"); }
        (*eP) = realloc((*eP), ents*esz);
        affirm((*eP) != NULL, "out of mem");
        (*entsP) = ents;
      }
  }

void sparse_matrix_trim(uint32_t *entsP, void **eP, uint32_t ents, size_t esz)
  { if (ents != (*entsP))
      { if (ents == 0)
          { free((*eP)); (*eP) = NULL; }
        else
          { (*eP) = realloc((*eP), ents*esz);
            affirm((*eP) != NULL, "out of mem");
          }
        (*entsP) = ents;
      }
  }

int sparse_matrix_compare_indices(uint32_t arow, uint32_t acol, uint32_t brow, uint32_t bcol, int orow, int ocol)
  { unsigned int zrow = (orow < 0 ? -orow : orow);
    unsigned int zcol = (ocol < 0 ? -ocol : ocol);
    demand(zrow != zcol, "ambiguous sorting criterion");
    if (zrow > zcol)
      { /* Row index is more important: */
        if (aP->row < bP->row)
          { return -trow; }
        else if (aP->row > bP->row)
          { return +trow; }
        if (zcol > 0)
          { /* Break ties by column index: */
            if (aP->col < bP->col)
              { return -tcol; }
            else if (aP->col > bP->col)
              { return +tcol; }
            else
              { return 0; }
          }
        else
          { return 0; }
      }
    else
      { /* Column index is more important: */
        if (aP->col < bP->col)
          { return -tcol; }
        else if (aP->col > bP->col)
          { return +tcol; }
        if (zrow > 0)
          { /* Break ties by row index: */
            if (aP->row < bP->row)
              { return -trow; }
            else if (aP->row > bP->row)
              { return +trow; }
            else
              { return 0; }
          }
        else
          { return 0; }
      }
  }

//@{sparse_matrix_io.c}
/* See sparse_matrix_io.h */
/* Last edited on 2008-03-29 14:33:12 by stolfi */

#include <sparse_matrix.h>
#include <sparse_matrix_io.h>
#include <stdlib.h>
#include <stdint.h>
#include <affirm.h>

#define sparse_matrix_FILE_VERSION "2008-07-24"

void sparse_matrix_read_header(FILE *rd, char *type, char **cmtP, uint32_t *rowsP, uint32_t *colsP, uint32_t *entsP)
  {
    /* Read header line: */
    filefmt_read_header(rd, type, sparse_matrix_FILE_VERSION);

    /* Read the comment lines, if any: */
    char *cmt = filefmt_read_comment(rd, '#');
    if (cmtP == NULL)
      { free(cmt); }
    else
      { (*cmtP) = cmt; }

    /* Read the matrix dimensions: */
    (*rowsP) = nget_uint32(rd, "rows", 10); fget_eol(rd);
    (*colsP) = nget_uint32(rd, "columns", 10); fget_eol(rd);
    (*entsP) =  nget_uint(rd, "entries", 10); fget_eol(rd);
  }

void sparse_matrix_read_footer(FILE *rd, char *type)
  {
    filefmt_read_footer(rd, type);
  }

void sparse_matrix_write_header(FILE *wr, char *type, char *cmt, uint32_t rows, uint32_t cols, uint32_t ents)
  {
    /* Write header line: */
    filefmt_write_header(wr, type, sparse_matrix_FILE_VERSION);

    /* Write comment lines, if any: */
    if (cmt != NULL) { filefmt_write_comments(wr, '#', cmt); }

    /* Write the matrix dimensions: */
    fprintf(wr, "rows = %u\n", rows);
    fprintf(wr, "columns = %u\n", cols);
    fprintf(wr, "entries = %u\n", ents);
  }

void sparse_matrix_write_header(FILE *wr, char *type)
  {
    filefmt_write_footer(wr, type);
    fflush(wr);
  }

//@{sparse_matrix_linalg.c}
/* See sparse_matrix_linalg.h */
/* Last edited on 2008-03-29 14:33:12 by stolfi */

#include <sparse_matrix.h>
#include <sparse_matrix_linalg.h>
#include <stdlib.h>
#include <stdint.h>
#include <affirm.h>

//@{sparse_matrix_junk.c}
/* TRIVIAL VALUE */

#define sparse_matrix_DECLARE_TRIVIAL_ELEM(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  ELEM_TYPE PREFIX##_trivial_elem(void)
/*
  This macro declares the function {{PREFIX}_trivial_elem}.

  The call {val = {PREFIX}_trivial_elem()} (with no arguments)
  will append into the variable {val} (of type {ELEM_TYPE})
  the value . */

#define sparse_matrix_DECLARE_ELEM_IS_TRIVIAL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
/*
  This macro declares the function {{PREFIX}_elem_is_trivial}.

  The call {b = {PREFIX}_elem_is_trivial(val)}
  will return the boolean value {TRUE} iff {val} is the
  element value which is not stored explicitly in the matrix. */

