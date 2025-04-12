/* Last edited on 2025-03-11 06:51:46 by stolfi */
/* Creating a deep copy of a {pst_gr_t} graph structure. */ 

#ifndef pst_gr_copy_H
#define pst_gr_copy_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <r2.h>
#include <pst_gr_path.h>
#include <pst_gr.h>

#include <pst_gr_copy.h>
  
pst_gr_t* pst_gr_copy(pst_gr_t* gr, int32_t ovid[], int32_t nvid[]);
  /* Creates a copy of {gr}, excluding all deleted vertices and edges.
    New copies are made of all records of valid vertices and edges,
    including labels and {haf_edge_t} records. In the new copy, vertices
    and edges are renumbered consecutively starting from 0. */

#endif


