/* Last edited on 2025-03-11 06:46:31 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_gr_integrate_recursive_H
#define pst_gr_integrate_recursive_H

/* Graph reduction by decimation of low-degree vertices. */ 

#include <stdint.h>
#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_integrate.h>
#include <pst_gr.h>
#include <pst_gr_integrate.h>

void pst_gr_integration_recursive
  ( pst_gr_t* gr,
    double Z[],
    int32_t level,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    pst_gr_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Computes a list of height values {Z[0..N-1]} for the vertices of
    the graph {gr}, based on the deltas and weights of its edges; where
    {N=gr->NV}.
    
    The heights computed by a multiscale version of the iterative
    Gauss-Seidel or Gauss-Jacobi method. The graph is recursively
    reduced by removing some vertices and rearranging the edges around
    them, as per {pst_gr_shrink}, until it is only two vertices.
    
    Then, at each scale, starting from the deepest recursion level and
    coming back up, the heights are computed with
    {pst_imgsys_solve_iterative} (q.v.) with parameters
    {para,verbose,level,reportStep}, and {reportHeights}. If {topoSort}
    is true, the equations will be sorted with
    {pst_imgsys_sort_equations}.
    
    At each scale, the initial guess for the Gauss-Seidel iteration is
    the solution obtained at the next lower (smaller) scale, with the
    heights of deleted edges estimated by interpolation of the
    neighbors.  The given values of {maxIter} and {convTol} are used at the
    highest (original) scale;they are ajdusted
    appropriately at each lower scale.
    
    The nput values of {Z[0..gr.NV-1]} are currently ignored.  On output,
    entries of {Z} that correspond to vertices of {gr} that 
    have degree zero are set to {NAN}. 
    
    The parameters {NX} and {NY} define the size of an image that 
    contains the pixel indices {gr.vdata[k].x} and {gr.vdata[k].y} 
    for all vertices.  */

#endif
