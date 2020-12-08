/* Operations on images for Tsai's routines. */
/* Last edited on 2011-05-15 00:41:28 by stolfi */

#ifndef tf_camera_image_H
#define tf_camera_image_H

#include <float_image.h>
#include <tf_camera.h>

float_image_t *image_apply_pincushion (tf_camera_params_t *cpar, float_image_t *img);

#endif
