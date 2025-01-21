/* Writing {pst_img_graph_t} structures. */ 
/* Last edited on 2025-01-09 17:09:49 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_write_H
#define pst_img_graph_write_H

#include <stdio.h>
#include <stdint.h>

#include <pst_img_graph.h>

void pst_img_graph_write(FILE* wr,pst_img_graph_t* g);
  /* Writes the graph {g} to {wr}, in a format compatible 
     with {pst_img_graph_read}. */

#endif


