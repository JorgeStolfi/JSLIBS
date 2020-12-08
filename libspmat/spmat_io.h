#ifndef spmat_io_H
#define spmat_io_H
/* Reading and writing generic sparse matrices */

#define spmat_io_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 20:19:40 by stolfi */

/*
  This interface provides basic I/O operations on generic sparse
  matrices.

  For example, if {dspmat_t} was defined as a sparse matrix of {double}
  values, one can write

    ------------------------------------------------------------
    //Read a sparse matrix {M} from {stdin}:
    dspmat_t M = dspmat_new(0,0,0);
    dspmat_read(stdin, &M);

    //Write the result to {stdout}:
    dspmat_write(stdout, &M);
    ------------------------------------------------------------ */
  /* The I/O functions are made available by the macros
  {spmat_io_def} and {spmat_io_impl} (see below).
  For example,

    ------------------------------------------------------------
    / * Declarations * /
    #include <spmat_io.h>

    spmat_typedef(dspmat_t, dspmat, double);

    spmat_io_def(dspmat_t, dspmat, double);
    ------------------------------------------------------------ */
  /*   ------------------------------------------------------------
    / * Implementations  * /
    #include <spmat_io.h>

    spmat_typedef(dspmat_t, dspmat, double);

    void dspmat_elem_read(FILE *rd, double *valP)
      { fscanf(rd, "%f", valP); }

    void dspmat_elem_write(FILE *wr, double *valP)
      { fprintf(wr, "%6.3f", *valP); }

    spmat_io_impl(dspmat_t, dspmat, double);
    ------------------------------------------------------------ */

#include <stdlib.h>
#include <stdint.h>

#include <spmat.h>

#define spmat_io_def(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_io_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the prototype declarations of
  four procedures, explained later:

     ------------------------------------------------------------
     PREFIX##_write
     PREFIX##_read
     PREFIX##_elem_write
     PREFIX##_elem_read
     ------------------------------------------------------------ */

#define spmat_io_impl(MATRIX_TYPE,PREFIX,ELEMT_TYPE) \
  spmat_io_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  These macros expand into implementations of the procedures

     ------------------------------------------------------------
     PREFIX##_write
     PREFIX##_read
     ------------------------------------------------------------
  
  BEFORE calling this macro, the client must provide adequate
  implementations for the functions {PREFIX##_elem_write} and
  {PREFIX##_elem_read}, as explained below. */

/* ======================================================================

  The remainder of this interface describes the functions provided
  by the macros {spmat_io_def} and {spmat_io_impl}.  */

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

  The macro {spmat_io_def} also declares the following
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
    any output that may be produced by {PREFIX##_elem_write}. */

#include <spmat_io_def.h>

#endif
