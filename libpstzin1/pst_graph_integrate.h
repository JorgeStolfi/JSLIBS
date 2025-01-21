#ifndef pst_graph_integrate_H
#define pst_graph_integrate_H

/* Last edited on 2025-01-11 14:57:22 by stolfi */
/* Created by Rafael F. V. Saracchini */
  
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <pst_graph.h>
#include <pst_imgsys.h>

pst_imgsys_t* pst_graph_build_system
  ( pst_graph_t* g,
    int32_t NX_Z,
    int32_t NY_Z
  );
  /* Builds the sparse equation system for the graph {g}. */

void pst_graph_solve_system
  ( pst_graph_t* g,
    pst_imgsys_t* S,
    float_image_t *OZ, 
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose 
  );
  /* Builds the sparse equation system for the graph {g},
    then solves it by weighted least squares. */

#endif
