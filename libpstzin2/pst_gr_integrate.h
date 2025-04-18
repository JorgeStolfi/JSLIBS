/* Last edited on 2025-03-11 06:48:45 by stolfi */
/* Building the linear equation system from a graph for slope-to-height integration. */

#ifndef pst_gr_integrate_H
#define pst_gr_integrate_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <r2.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_gr.h>

typedef void pst_gr_integrate_report_data_proc_t(int32_t level, pst_gr_t *gr); 
  /* Type of a client-given procedure that may be called by recursive integrators
    to report the input graph at each scale. */   

pst_imgsys_t *pst_gr_integrate_build_system
  ( pst_gr_t* gr,
    int32_t NX_Z,
    int32_t NY_Z,
    int32_t kz_from_kv[],
    bool_t verbose
  );
  /* Builds an equation system {S} from the graph {gr}. 
  
    Each vertex of {gr} with index {kv} has a corresponding (unknown)
    height {Z[kz]} and a corresponding equation {S.eq[kz]}; except for
    vertices marked {DELETED} and vertices of degree zero. The indices
    {kz} are assigned sequentially in {0..S.N-1}, skipping the omitted
    vertices. The correspondence is returned in the table
    {kz_from_kv[0..gr.NV-1]}, with {-1} marking the omitted vertices A
    vertex is considered deleted if and only if {gr.vdata[kv].vmark} is
    {DELETED}.
    
    The fields {S.NX} and {S.NY} are set to {NX_Z} and {NY_Z}, respectively.
    The pixel indices {x,y} stored in each vertex data record of {gr}
    must be either {-1} or in the ranges {0..NX-1} and {0..NY-1},
    respectively, and are stored in {S.col[kz]} and {S.row[kz]}, and
    used to build the inverse table {S.uid[0..S->N-1]}.
    
    The edge weights of {gr} must be all finite and positive. */

#endif
