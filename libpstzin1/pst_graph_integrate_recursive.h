#ifndef pst_graph_integrate_recursive_H
#define pst_graph_integrate_recursive_H

/* Last edited on 2025-01-14 17:26:56 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

#include <pst_graph.h>

void pst_graph_integration_recursive( 
  pst_graph_t* g,
  float_image_t* OZ,
  uint32_t maxIter,
  double convTol, 
  bool_t para, 
  bool_t szero, 
  bool_t verbose,
  int32_t level
);


#endif
