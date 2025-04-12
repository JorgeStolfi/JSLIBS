/* See {pst_gr_from_slope_map.h} */
/* Last edited on 2025-03-14 19:18:07 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <affirm.h>

#include <pst_gr.h>

#include <pst_gr_curl_image.h>
  
void pst_gr_curl_image_fill(pst_gr_t *gr, float_image_t *U)
  { demand((gr->NX > 0) && (gr->NY > 0), "graph has no associated map");
    int32_t NC_U, NX_U, NY_U;
    float_image_get_size(U, &NC_U, &NX_U, &NY_U);
    demand(NC_U == 1, "invalid channels in curl image");
    demand((NX_U == gr->NX-1) && (NY_U == gr->NY-1), "wrong {U} image size");
    
    float_image_fill_channel(U, 0, 0.0);
    for (int32_t y = 0; y < NY_U; y++)
      { for (int32_t x = 0; x < NX_U; x++)
          { pst_gr_arc_t ai = pst_gr_find_enclosing_face(gr, x, y);
            if (ai != pst_gr_NONE)
              { double curl = pst_gr_compute_left_face_curl(gr, ai);
                float_image_set_sample(U, 0,x,y, (float)curl);
              }
          }
      }
  }
              
