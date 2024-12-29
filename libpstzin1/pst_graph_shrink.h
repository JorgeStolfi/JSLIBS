#ifndef pst_graph_shrink_H
#define pst_graph_shrink_H

/* Last edited on 2024-12-23 09:48:48 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>
#include <pst_graph.h>

pst_graph_t* pst_graph_shrink(pst_graph_t* g);
  /* Createsa graph {h} that is a shrunk version of {g}.  Roughly every
    two vertices of {g} become one vertex of {h}. */

#endif
