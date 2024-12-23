/* Last edited on 2024-12-23 05:54:56 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_integrate_recursive_H
#define pst_img_graph_integrate_recursive_H

/* Improved version of {pst_graph_integrate_recursive.h} that uses {oct.h} for the topology info. */ 

#include <oct.h>
#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>
#include <pst_img_graph_integrate.h>

#define MARK_VERTEX_NONE 0
#define MARK_VERTEX_REMOVED 1
#define MARK_VERTEX_PRESERVED 2

#define DEFAULT_WMAG 1.0

void pst_img_graph_mark_vertex_removal(pst_img_graph_t* g, uint32_t degree[]);

void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t* g);

void pst_img_graph_shrink(pst_img_graph_t* g,double* w_i,double wmag);

void pst_img_graph_vertex_remove
  ( pst_img_graph_t* g,
    uint32_t vi,
    double* w_i,
    double wmag,
    bool_t verbose
  );
  
void pst_img_graph_remove_paralel_edges(pst_img_graph_t* g);

void pst_img_graph_copy_solution_from_shrunk(pst_img_graph_t* jg,double* jZ,pst_img_graph_t* ig, double *iZ);
  /* Also marks vertices of {ig} as MARK_VERTEX_PRESERVED or MARK_VERTEX_REMOVED*/

void pst_img_graph_estimate_from_shrunk(pst_img_graph_t* jg,double* jZ,pst_img_graph_t* g, double *iZ);

void pst_img_graph_integration_recursive
  ( pst_img_graph_t* g,
    double* iZ,
    double* iW,
    double wmag,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose,
    uint32_t level,
    float_image_t* OZ, /*debug only*/
    float_image_t* RZ,
    char* out_prefix,
    bool_t debug
  );

void pst_img_graph_vertex_remove_general
  ( pst_img_graph_t* g,
    uint32_t vi,
    uint32_t n_neighbours,
    double* w_i,
    double wmag,
    bool_t merge_diagonals,
    bool_t verbose
  );

#endif
