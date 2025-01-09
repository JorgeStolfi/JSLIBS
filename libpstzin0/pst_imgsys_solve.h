/* tools for solving of a {pst_imgsys_t} linear equation system. */

#ifndef pst_imgsys_solve_H
#define pst_imgsys_solve_H

/* Created on 2005-12-04 by Jorge Stolfi, unicamp, <stolfi@ic.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-08 00:21:23 by stolfi */
/* See the copyright and authorship notice at the end of this file. */

#include <stdint.h>

#include <bool.h>

#include <pst_imgsys.h>

uint32_t* pst_imgsys_sort_equations(pst_imgsys_t *S);
  /* Returns and array {ord[0..S->N-1]} with the indices of the equations 
     of {S} in order of increasing {wtot} field. */ 

typedef void pst_imgsys_solve_report_sol_proc_t(uint32_t indent, uint32_t iter, double change, bool_t final, uint32_t N, double h[]);
  /* Type of a procedure that is used to report the progress of the solution. */

void pst_imgsys_solve_iterative
  ( pst_imgsys_t *S, 
    double h[],
    uint32_t ord[],
    uint32_t maxIter, 
    double convTol,
    bool_t para, 
    bool_t szero,
    bool_t verbose,
    uint32_t indent,
    uint32_t reportStep,
    pst_imgsys_solve_report_sol_proc_t *reportSol
  );
  /* Solves system {S}, and stores the solution into the vector {h}.  

    Uses an iterative method, and therefore assumes that the unknown
    {h[i]} is the first term in equation {S->eq[i]}, for {i =
    0..N-1} and {N = S->N}; and that its coefficient is ``large enough''. Upon
    entry, the {h} vector must contain the starting guess. Executes
    at most {max_iter} iterations, but stops whenever two consecutive
    iterations do not change any variable by more than the tolerance
    {tol}.

    If {para} is true, uses the ``parallel'' (Jacobi) variant of the method,
    otherwise uses the sequential (Gauss-Seidel) variant.
    If {ord} is not NULL, solves the equations in the order {ord[0..N-1]}
    If {szero} is true, adjusts the solution so that it adds to zero,
    after each iteration.

    If {verbose} is true, prints information about the iterations to {stderr}.
    All messages will be indented by {indent} spaces.

    if {reportSol} is not null, the procedure calls
    {reportSol(indent,iter,change,final,N,h) one or more times during the
    solution. Here {iter} is the number of complete Gauss-Seidel or
    Gauss-Jacobi iterations performed before the call; and {change} is
    the max absolute change in any {h} element since the previous
    iteration (meaningless when {iter=0}). The procedure {reportSol} is
    always called once with {final=TRUE} after the last iteration, and
    once with {iter=0} and {final=FALSE} before the first iteration. If
    {reportStep} is not zero, {reportSol} is also called with {final=FALSE}
    before each iteration whose index {iter} is a positive multiple of
    {reportStep}. */

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
