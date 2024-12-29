/* imgsys - linear system solving for image-related problems. */

#ifndef pst_imgsys_H
#define pst_imgsys_H

/* Created on 2005-12-04 by Jorge Stolfi, unicamp, <stolfi@ic.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2024-12-24 18:57:06 by stolfi */
/* See the copyright and authorship notice at the end of this file. */

#include <float_image.h>
#include <pst_height_map.h>

#define MAXCOEFS 5

typedef struct pst_imgsys_equation_t
  { uint32_t nt;             /* {nt} is the number of terms in the equation. */
    uint32_t ix[MAXCOEFS];   /* {ix[k]} is the index of some unknown in the equation. */
    double cf[MAXCOEFS];     /* Coefficient of that unknown. */
    double rhs;              /* Right-hand side of equation. */
    double wtot;             /* Weight of the equation and of its main variable. */
  } pst_imgsys_equation_t;
  /* An {pst_imgsys_equation_t} record {eq} represents a linear equation 
    with at most {MAXCOEFS} nonzero terms, namely
      {SUM {eq.cf[j]*h[eq.ix[j]] : j = 0..eq.nt-1} == rhs},
    where {h[j]} is the unknown with index {j}.
    The variable with index {eq.ix[0]} is the main variable of the equation. */

typedef struct pst_imgsys_t
  { uint32_t N;                 /* Number of unknowns and equations. */
    pst_imgsys_equation_t* eq;  /* The equations are {eq[0..N-1]}. */
    /* Debugging info: */
    uint32_t *col;              /* Maps index of an unknown/equation to a {x} index. */
    uint32_t *row;              /* Maps index of an unknown/equation to a {y} index. */
    int32_t *ix;                /* Maps pixel index {x + y*NX_Z} to index of unknown/equation, or -1 if none. */
  } pst_imgsys_t;
  /* An {pst_imgsys_t} represents {N} linear equations {eq[0..N-1]} on {N} 
    unknowns {h[0..N-1}. In a proper system, the main variable
    of equation {k} has index {k} --- that is, {eq[k].ix[0] == k} ---
    for all {k} in {0..N-1}.  The tables {col,row,ix} define a mapping
    from equation/unknown indices to points of a 2D integer grid;
    they are used only when printing the system. */

pst_imgsys_t *pst_imgsys_new(int32_t NX, int32_t NY, uint32_t N, int32_t *ix, uint32_t *col, uint32_t *row);
  /* Creates a new linear system {S} with {N} equations on {N} unknowns
    for an image with {NX} columns and {NY} rows, using the given
    index mapping tables. Allocates the equations {S.eq} but does not
    initialize them. */

pst_imgsys_t *pst_imgsys_from_eqs(int32_t NX, int32_t NY, uint32_t N, pst_imgsys_equation_t *eq, int32_t *ix, uint32_t *col, uint32_t *row);
  /* Like {pst_imgsys_new} but uses the given {eq} vector instead of allocating a new one. */

void pst_imgsys_free(pst_imgsys_t *S);
  /* Deallocates all storage used by {S}, including the index tables. */

typedef void pst_imgsys_solution_report_proc_t(uint32_t iter, double change, bool_t final, uint32_t N, double Z[]);
  /* Type of a procedure that is used to report the progress of the solution. */

uint32_t* pst_imgsys_sort_equations(pst_imgsys_t *S);
  /* Returns and array {ord[0..S->N-1]} with the indices of the equations 
     of {S} in order of increasing {wtot} field. */ 

void pst_imgsys_solve
  ( pst_imgsys_t *S, 
    double Z[],
    uint32_t ord[],
    uint32_t maxIter, 
    double convTol,
    bool_t para, 
    bool_t szero,
    bool_t verbose,
    uint32_t indent,
    pst_imgsys_solution_report_proc_t *reportSol
  );
  /* Solves system {S}, and stores the solution into the vector {Z}.  
     
     Uses an iterative method, and therefore assumes that the unknown
     {Z[i]} appears in some term of equation {S->eq[i]}, where {i =
     0..N-1} and {N = S->N}, and that its coefficient is ``large enough''. Upon
     entry, the {OZ} image must contain the starting guess. Executes
     at most {max_iter} iterations, but stops whenever two consecutive
     iterations do not change any variable by more than the tolerance
     {tol}.

     If {para} is true, uses the ``parallel'' variant of the method
     (Jacobi), otherwise uses the sequential variant (Gauss-Seidel).
     If {ord} is not NULL, solves the equations in the order {ord[0..N-1]}
     If {szero} is true, adjusts the solution so that it adds to zero,
     after each iteration.
     
     If {verbose} is true, prints information about the iterations to {stderr}.
     All messages will be indented by {indent} spaces.
     
     if {reportSol} is not null, will call
     {reportSol(iter,change,FALSE,N,Z)} before every iteration,
     {reportSol(iter,change,TRUE,N,Z)} after the last one; where
     {iter} is the number of iterations done so far, and {change} is
     the maximum absolute change in any {Z} value from the previous
     iteration. When {iter} is zero, {change} is irrelevant.
     If {maxIter} is zero, makes only one call with {iter=0,final=TRUE}.
     When {maxIter>0}, the first call has {iter=0,final=FALSE}
     and the last call has {iter>0,final=TRUE}*/
          
/* DEBUGGING */
    
typedef void pst_imgsys_report_proc_t(uint32_t level, pst_imgsys_t *S); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the system used at each scale.
    Uses {col} and {row} to map indices of unknowns to pixel indices. */   

/* I/O */
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S);
  /* Writes the system {S} to stream {wr}.  The indices of unknowns and 
    equations are mapped to column and row indices
    with the tables {col[0..S.N-1]} and {row[0..S.N-1]}.  */

void pst_imgsys_write_report(pst_imgsys_t *S, char *filePrefix, int32_t level, char *tag, uint32_t indent);
  /* Writes the system {S} to a file called
    "{filePrefix}-{level}-{tag}.sys". If {tag} is null or empty the
    "-{tag}" is omitted. Uses {pst_imgsys_write}. Diagnostic messages
    are indented by {indent} spaces. */

#endif

/*
**
** Copyright (C) Jorge Stolfi, Unicamp.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty. Neither the author nor its employers are liable to
** any damages which may result from its use.
*/
