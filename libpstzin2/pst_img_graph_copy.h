/* Last edited on 2025-01-11 17:53:50 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_copy_H
#define pst_img_graph_copy_H

/* Routine to make a deep copy of a {pst_img_graph_t}. */ 

#include <haf.h>
#include <r2.h>
#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_path.h>

#include <pst_img_graph.h>
  
pst_img_graph_t* pst_img_graph_copy(pst_img_graph_t* g, int32_t ovid[], int32_t nvid[]);
  /* Creates a copy of {g}, excluding all deleted vertices and edges.
    New copies are made of all records of valid vertices and edges,
    including labels and {haf_edge_t} records. In the new copy, vertices
    and edges are renumbered consecutively starting from 0. */

#endif


