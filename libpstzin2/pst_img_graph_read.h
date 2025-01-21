/* Reading {pst_img_graph_t} structures. */ 
/* Last edited on 2025-01-09 17:09:43 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_read_H
#define pst_img_graph_read_H

#include <stdio.h>
#include <stdint.h>

#include <pst_img_graph.h>

pst_img_graph_t* pst_img_graph_read(FILE* rd);
  /* Reads from {rd} a descripton of a graph {g}. */

#endif


