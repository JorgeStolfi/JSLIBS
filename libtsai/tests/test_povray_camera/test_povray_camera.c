/* Last edited on 2023-04-23 11:09:06 by stolfi */

#define PROG_NAME "test_povray_camera"
#define PROG_DESC "tests the conversion from Tsai camera matrix to POV-Ray camera spec"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <r2x2.h> 
#include <r2.h>
#include <jsfile.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_write_pnm.h>

#include <tf_calib.h> 
#include <tf_camera.h> 

#define tpc_PNM_GAMMA 1.000
#define tpc_PNM_BIAS  0.000

void generate_artificial_cpar (tf_camera_params_t *cpar, double clock);

void set_generic_cpar (tf_camera_params_t *cpar);

void set_special_cpar (tf_camera_params_t *cpar);

int32_t main (int32_t argc, char *argv[])
{

  /* Get the fixed camera parameters: */
  tf_camera_specs_t *cspec = tf_camera_specs_for_povray_svga ();
  tf_camera_specs_print(stderr, cspec);

  tf_camera_params_t *cpar = tf_camera_specs_get_new_mean_params(cspec);
  float_image_t *img = float_image_new(1, 640, 480);
  float_image_fill_channel(img, 0, 1.0);

  //set_generic_cpar (cpar);
  set_special_cpar(cpar);
  tf_camera_params_print(cpar, stderr);

  /*Write the transformation matrix to a pov ray file*/
  FILE *inc = fopen("camera_dat.inc", "w");
  tf_camera_params_write_povray (inc, cpar, FALSE);
  fclose(inc);

  /*Write 4 points in a image file*/

  auto void paint_cross(double xw, double yw, double zw);

  void paint_cross(double xw, double yw, double zw) {
     r3_t pw = (r3_t){{ xw, yw, zw }};
     r2_t pi = tf_world_coords_to_image_coords (cpar, pw);
     fprintf(stderr, "cross at ( %10.6f %10.6f %10.6f ) world = ( %10.6f %10.6f ) image\n",
	     xw, yw, zw, 
             pi.c[0], pi.c[1]
	     );
     bool_t empty = FALSE;
     bool_t diagonal = TRUE;
     (void)float_image_paint_cross (img, 0, pi.c[0], pi.c[1], 5.0, empty, 1.0, diagonal, 0.0, 3);	     
  }
  
  paint_cross(0.0, 0.0, 0.0);
  paint_cross(500.0, 0.0, 0.0);
  paint_cross(0.0, 1200.0, 0.0);
  paint_cross(0.0, 0.0, 900.0);

  float_image_write_pnm_named ("out/crosses.pgm", img, FALSE, tpc_PNM_GAMMA, tpc_PNM_BIAS, FALSE, TRUE, FALSE);

  return 0;
}

void set_generic_cpar (tf_camera_params_t *cpar)
{
  cpar->S.c[0][0] = 1.0; 
  cpar->S.c[0][1] = 0.0; 
  cpar->S.c[0][2] = 0.0; 
  cpar->S.c[0][3] = 0.0;
  cpar->S.c[1][1] = -0.18013389693448786; 
  cpar->S.c[1][2] = 0.98223960595054982; 
  cpar->S.c[1][3] = -0.05250843434443486;
  cpar->S.c[2][1] = 0.23375941413578807; 
  cpar->S.c[2][2] = -0.00910513833225995; 
  cpar->S.c[2][3] = -0.9722518360789269;
  cpar->S.c[3][1] = -0.95546235691317971; 
  cpar->S.c[3][2] = -0.18740985288415049;  
  cpar->S.c[3][3] = -0.22796761077805028;
  cpar->S.c[1][0] = 39.32695350295859527; 
  cpar->S.c[2][0] = 287.21780834354387935;
  cpar->S.c[3][0] = 1588.9007105221079828;
  cpar->f = 41.28486083472594004;
  cpar->kappa = 0.00025220772348791;
  cpar->sx = 1.00499999999999989;
}

void set_special_cpar (tf_camera_params_t *cpar)
{
  cpar->S.c[0][0] = 1.0; 
  cpar->S.c[0][1] = 0.0; 
  cpar->S.c[0][2] = 0.0; 
  cpar->S.c[0][3] = 0.0;

  cpar->S.c[1][0] = 0.0; 
  cpar->S.c[1][1] = 0.0; 
  cpar->S.c[1][2] = 1.0; 
  cpar->S.c[1][3] = 0.0;

  cpar->S.c[2][0] = 0.0;
  cpar->S.c[2][1] = 0.0; 
  cpar->S.c[2][2] = 0.0; 
  cpar->S.c[2][3] = -1.0;

  cpar->S.c[3][0] = 4000.0;
  cpar->S.c[3][1] = -1.0; 
  cpar->S.c[3][2] = 0.0;  
  cpar->S.c[3][3] = 0.0;

  cpar->f = 40.0;
  cpar->kappa = 0.0;
  cpar->sx = 1.0;
}
