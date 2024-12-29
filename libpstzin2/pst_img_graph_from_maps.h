/* Last edited on 2024-12-25 08:14:43 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_from_maps_H
#define pst_img_graph_from_maps_H

/* Improved version of {pst_graph_from_maps.h} that uses {haf.h} for the topology info. */ 

#include <r2.h>
#include <bool.h>
#include <float_image.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>

pst_img_graph_t* pst_img_graph_create_from_gradient_and_weight_maps(float_image_t* IG, float_image_t* IW, bool_t add_diags);

#endif
