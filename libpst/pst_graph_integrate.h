#ifndef pst_graph_integrate_H
#define pst_graph_integrate_H

/* Last edited on 2024-12-22 19:17:52 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <pst_graph.h>
#include <pst_imgsys.h>

pst_imgsys_t* pst_graph_build_integration_system(pst_graph_t* g, uint32_t NX_Z, uint32_t NY_Z);

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

#endif
