/* Last edited on 2025-04-27 02:57:20 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef psgr_integrate_iterative_H
#define psgr_integrate_iterative_H

/* Single-scale graph based integration of slope maps. */ 

#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_integrate.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_integrate.h>

void psgr_integrate_iterative
  ( psgr_t *gr,
    double Z[],
    bool_t sortSys,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    int32_t level,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Computes a list of height values {Z[0..N-1]} for the vertices of
    the graph {gr}, based on the deltas and weights of its edges; where
    {N=gr->NV}.
    
    The vector {Z} is computed with the iterative Gauss-Seidel or
    Gauss-Jacobi method as implemented by {pst_imgsys_solve_iterative}
    (q.v.) with parameters {maxIter,convTol,para,szero}. as well as
    {verbose,level,reportStep}, and {reportHeights}. If {sortSys} is
    true, the equations will be sorted with {pst_imgsys_sort_equations}.
    On input, {Z} should contain an initial guess for the solution. 
    
    Entries of {Z} that correspond to vertices of {gr} that 
    have degree zero are set to {NAN}.
    
    If {reportSys} is not {NULL}, it is called once,
    after the equation system {S} has been built from the graph {gr}. */

#endif
