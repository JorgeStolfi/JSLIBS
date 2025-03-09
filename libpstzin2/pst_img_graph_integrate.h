/* Last edited on 2025-02-11 10:19:55 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_integrate_H
#define pst_img_graph_integrate_H

/* Improved version of {pst_graph_integrate.h} that uses {haf.h} for the topology info. */ 

#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>

typedef void pst_img_graph_integrate_report_data_proc_t(int32_t level, pst_img_graph_t *g); 
  /* Type of a client-given procedure that may be called by recursive integrators
    to report the input graph at each scale. */   

pst_imgsys_t *pst_img_graph_integrate_build_system
  ( pst_img_graph_t* g,
    int32_t NX_Z,
    int32_t NY_Z,
    int32_t kz_from_kv[],
    bool_t verbose
  );
  /* Builds an equation system {S} from the graph {g}. 
  
    Each vertex of {g} with index {kv} has a corresponding (unknown)
    height {Z[kz]} and a corresponding equation {S.eq[kz]}; except for
    vertices marked {DELETED} and vertices of degree zero. The indices
    {kz} are assigned sequentially in {0..S.N-1}, skipping the omitted
    vertices. The correspondence is returned in the table
    {kz_from_kv[0..g.NV-1]}, with {-1} marking the omitted vertices A
    vertex is considered deleted if and only if {g.vdata[kv].vmark} is
    {DELETED}.
    
    The fields {S.NX} and {S.NY} are set to {NX_Z} and {NY_Z}, respectively.
    The pixel indices {x,y} stored in each vertex data record of {g}
    must be either {-1} or in the ranges {0..NX-1} and {0..NY-1},
    respectively, and are stored in {S.col[kz]} and {S.row[kz]}, and
    used to build the inverse table {S.uid[0..S->N-1]}.
    
    The edge weights of {g} must be all finite and positive. */

#endif
