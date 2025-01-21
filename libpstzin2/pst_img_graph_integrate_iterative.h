/* Last edited on 2025-01-14 17:26:23 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_integrate_iterative_H
#define pst_img_graph_integrate_iterative_H

/* Single-scale graph based integration of slope maps. */ 

#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_integrate.h>
#include <pst_img_graph.h>
#include <pst_img_graph_integrate.h>

void pst_img_graph_integrate_iterative
  ( pst_img_graph_t *g,
    double Z[],
    bool_t topoSort,
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
    the graph {g}, based on the deltas and weights of its edges; where
    {N=g->NV}.
    
    The vector {Z} is computed with the iterative Gauss-Seidel or
    Gauss-Jacobi method as implemented by {pst_imgsys_solve_iterative}
    (q.v.) with parameters {maxIter,convTol,para,szero}. as well as
    {verbose,level,reportStep}, and {reportHeights}. If {topoSort} is
    true, the equations will be sorted with {pst_imgsys_sort_equations}.
    On input, {Z} should contain an initial guess for the solution. 
    
    Entries of {Z} that correspond to vertices of {g} that 
    have degree zero are set to {NAN}.
    
    If {reportSys} is not {NULL}, it is called once,
    after the equation system {S} has been built from the graph {g}. */

#endif
