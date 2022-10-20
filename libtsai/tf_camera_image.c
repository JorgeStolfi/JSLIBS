/* See {tf_camera_image.h}.  */
/* Last edited on 2022-10-20 05:52:58 by stolfi */

#define _GNU_SOURCE

#include <float_image.h>
#include <stdint.h>
#include <float_image_interpolate.h>
#include <affirm.h>

#include <tf_camera.h>
#include <tf_camera_image.h>

float_image_t *image_apply_pincushion (tf_camera_params_t *cpar, float_image_t *img)
  {
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    
    float_image_t *omg = float_image_new(NC, NX, NY);
    double vv[NC];
    float v[NC];

    int32_t iy, ix, ic;
    for (iy = 0; iy < NY; iy++) {
      for (ix = 0; ix < NX; ix++) {

        r2_t pi;
        pi.c[0] = iy;
        pi.c[1] = ix;
        r2_t pd = tf_image_coords_to_sensor_coords (cpar, pi);
        r2_t pu = tf_dis_sensor_coords_to_und_sensor_coords (cpar, pd);
        r2_t npi = tf_sensor_coords_to_image_coords (cpar, pu);
 
        float_image_interpolate_pixel(img, npi.c[0], npi.c[1], 1, ix_reduction_EXTEND, vv);
        for (ic = 0; ic < NC; ic++) { v[ic] = (float)vv[ic]; }
        float_image_set_pixel(omg, ix, iy, v); 
      }
    }
    return omg;
  }
