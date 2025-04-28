/* See {psgr_from_slope_map.h} */
/* Last edited on 2025-04-27 03:08:58 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <affirm.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>

#include <psgr_curl_image.h>
  
#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE
  /* Shorter names. */  

void psgr_curl_image_fill(psgr_t *gr, float_image_t *U)
  { i2_t sz = psgr_image_size(gr); int32_t NX_gr = sz.c[0], NY_gr = sz.c[1];
    demand((NX_gr > 0) && (NY_gr > 0), "graph has no associated map");
    int32_t NC_U, NX_U, NY_U;
    float_image_get_size(U, &NC_U, &NX_U, &NY_U);
    demand(NC_U == 1, "invalid channels in curl image");
    demand((NX_U == NX_gr-1) && (NY_U == NY_gr-1), "wrong {U} image size");
    
    float_image_fill_channel(U, 0, 0.0);
    for (int32_t y_U = 0; y_U < NY_U; y_U++)
      { for (int32_t x_U = 0; x_U < NX_U; x_U++)
          { /* Map the center of the pixel {[x_U,y_U]} of {U} to the {gr} domain: */
            r2_t p = (r2_t){{ x_U + 0.5, y_U + 0.5 }};
            psgr_arc_t ai = psgr_find_enclosing_face(gr, &p);
            if (ai.ix  != ANONE.ix)
              { double curl = psgr_arc_left_face_curl(gr, ai);
                float_image_set_sample(U, 0,x_U,y_U, (float)curl);
              }
          }
      }
  }
              
