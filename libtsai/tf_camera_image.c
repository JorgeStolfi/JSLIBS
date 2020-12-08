/* See {tf_camera_image.h}.  */
/* Last edited on 2015-10-04 20:50:29 by stolfilocal */

#define _GNU_SOURCE

#include <float_image.h>
#include <float_image_interpolate.h>
#include <affirm.h>

#include <tf_camera.h>
#include <tf_camera_image.h>

float_image_t *image_apply_pincushion (tf_camera_params_t *cpar, float_image_t *img)
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    float_image_t *omg = float_image_new(NC, NX, NY);
    double vv[NC];
    float v[NC];

    int iy, ix, ic;
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
