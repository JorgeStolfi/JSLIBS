/* Last edited on 2024-12-24 18:57:53 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_integrate_H
#define pst_img_graph_integrate_H

/* Improved version of {pst_graph_integrate.h} that uses {haf.h} for the topology info. */ 

#include <r2.h>
#include <float_image.h>

#include <pst_img_graph.h>
#include <pst_imgsys.h>

#include <pst_img_graph.h>

pst_imgsys_t *pst_img_graph_build_integration_system
  ( pst_img_graph_t* g,
    double *iW,
    int32_t NX_Z,
    int32_t NY_Z,
    int32_t **ref_tab
  );

void pst_img_graph_solve_system
  ( pst_img_graph_t *g,
    pst_imgsys_t *S,
    double *iZ,
    double *iW,
    int32_t *ref_tab,
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose 
  );

void pst_img_graph_integration
  ( pst_img_graph_t *g,
    double *iZ,
    double *iW,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose,
    uint32_t level,
    float_image_t *OZ, /*debug only*/
    float_image_t *RZ,
    char* out_prefix
  );

void pst_img_graph_put_solution_into_image
  ( pst_img_graph_t* g,
    double *iZ,
    float_image_t* OZ
  );

void pst_img_graph_put_error_into_image
  ( pst_img_graph_t *g,
    double *iZ,
    double *iW,
    float_image_t *RZ,
    float_image_t *OZ
  );

uint32_t *pst_img_graph_sort_equations
  ( pst_imgsys_t *S,
    int32_t *ref_tab,
    double *iW ,
    uint32_t NV 
  );

#endif
