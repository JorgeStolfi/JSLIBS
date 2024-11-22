/* Last edited on 2024-11-20 06:10:01 by stolfi */

#define PROG_NAME "test_povray_camera"
#define PROG_DESC "tests the conversion from Tsai camera matrix to POV-Ray camera spec"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <affirm.h> 
#include <r2x2.h> 
#include <r2.h>
#include <jsfile.h>
#include <fget.h>
#include <float_image.h>
#include <float_image_write_pnm.h>
#include <float_image_paint.h>

#include <tf_calib.h> 
#include <tf_camera.h> 
#include <tf_camera_specs.h> 
#include <tf_targets.h> 

#define ttm_PNM_GAMMA 1.000
#define ttm_PNM_BIAS  0.000

int32_t main (int32_t argc, char *argv[])
  {
    /* Arguments: */
    char *data_dir = argv[1];      /* Directory for input file names. */
    char *out_dir = argv[2];       /* Directory for output file names */
    char *calib_tag = argv[3];     /* Tag of subject calibration, or "NONE". */
    char *true_tag = argv[3];      /* Tag of reference calibration, or "NONE". */
    int32_t frame_ini = atoi(argv[4]); /* Number of initial frame. */
    int32_t frame_fin = atoi(argv[5]); /* Number of final frame. */

    /* Get the fixed camera parameters: */
    char *cspec_fname = jsprintf("%s/cspec.txt", data_dir);
    FILE *cspec_file = open_read(cspec_fname, TRUE);
    char *cname = fget_string(cspec_file);
    fclose(cspec_file);

    tf_camera_specs_t *cspec = NULL;
    if (strcmp(cname,"optura") == 0)
      { cspec = tf_camera_specs_for_canon_optura (); }
    else if (strcmp(cname,"svga") == 0)
      { cspec = tf_camera_specs_for_povray_svga (); }
    else if (strcmp(cname,"hvga") == 0)
      { cspec = tf_camera_specs_for_povray_hvga (); }
    else
      { fprintf(stderr, "error: unknown camera type\n");
        exit(1);
      }  

    /* Allocate and initialize the camera parameters: */
    tf_camera_params_t *cal_cpar = tf_camera_specs_get_new_mean_params(cspec);
    tf_camera_params_t *ref_cpar = tf_camera_specs_get_new_mean_params(cspec);

    /* Get the world coordinates of fiducial marks: */
    tf_calib_data_t *cdat = tf_calib_data_read(data_dir);

    /* Loop on frames: */
    int32_t frame_num;
    for (frame_num = frame_ini; frame_num <= frame_fin; frame_num++)
    {
      fprintf(stderr, "--- begin frame %05d --------------------------------------------------\n", frame_num);
      
      /* Get frame-specific input and output directories: */
      char *frame_data_dir = jsprintf("%s/%05d", data_dir, frame_num);
      char *frame_out_dir = jsprintf("%s/%05d", out_dir, frame_num);
      assert(mkdir(frame_out_dir, 0777) == 0);
      
      /* Read the input frame image: */
      float_image_t *img = read_frame_image(frame_data_dir);

      /* Get the computed and reference calibrations of this frame: */
      if (strcmp(cal_tag,"NONE") != 0)
        { read_camera_parameters("calibrated", frame_data_dir, cal_tag, cal_cpar); }
      if (strcmp(ref_tag,"NONE") != 0)
        { read_camera_parameters("reference",  frame_data_dir, ref_tag, ref_cpar); }

      char *out_f_pi_name = jsprintf("%s/p_i.txt", frame_out_dir);
      FILE *f_pi = fopen(out_f_pi_name, "w");
      free(out_f_pi_name); 

      char *out_f_pw_name = jsprintf("%s/p_w.txt", frame_out_dir);
      FILE *f_pw = fopen(out_f_pw_name, "w");
      free(out_f_pw_name); 

      char *out_f_p_wgt_name = jsprintf("%s/p_wgt.txt", frame_out_dir);
      FILE *f_p_wgt = fopen(out_f_p_wgt_name, "w");
      free(out_f_p_wgt_name); 

      fprintf(f_pi, "%d\n", n_p_w);
      fprintf(f_pw, "%d\n", n_p_w);
      fprintf(f_p_wgt, "%d\n", n_p_w);

      /* Process marks: */
      for (int32_t i = 0; i < n_p_w; i++) {
          fprintf(stdout, "before tf_world_coords_to_image_coords\n"); 
          r2_t p_i = tf_world_coords_to_image_coords (cal_cpar, cdat->world[i]);
          fprintf(stdout, "after tf_world_coords_to_image_coords\n"); 
          
          double mark_radius = 5;
          bool_t empty = FALSE;
          bool_t diagonal = TRUE;
          (void)float_image_paint_cross(img, 0, p_i.c[0], p_i.c[1], mark_radius, empty, 1.0, diagonal, 1.0, 3); 
          (void)float_image_paint_cross(img, 0, p_i.c[0], p_i.c[1], mark_radius, empty, 0.5, diagonal, 0.0, 3); 

          fprintf(f_pi, "%f %f\n", p_i.c[0], p_i.c[1]);
          fprintf(f_pw, "%f %f %f\n", cdat->world[i].c[0], cdat->world[i].c[1], cdat->world[i].c[2]);
          fprintf(f_p_wgt, "%f\n", 1.00000);

      }

      fclose(f_pi);
      fclose(f_pw);
      fclose(f_p_wgt);

      fprintf(stdout, "loop finished\n"); 

      char *cross_fname = jsprintf("out/cross_%05d.pgm", frame_frame_num);
      float_image_write_pnm_named(cross_fname, img, FALSE, ttm_PNM_GAMMA, ttm_PNM_BIAS, FALSE, TRUE, FALSE);
      free(cross_fname);

      free(frame_fname);
      float_image_free(img);
      fprintf(stderr, "--- end frame %05d --------------------------------------------------\n", frame_num);
    }

    if (strcmp(f_true_fname, "none") != 0) {
       fclose(f_true);
    }

    return 0;
  }

        
      char *f_cpar_fname = jsprintf("%s/calib-%s.cpar", frame_data_dir, calib_tag);
      FILE *f_cpar_in = open_read(f_cpar_fname, TRUE);
      tf_camera_params_read(f_cpar_in, cpar);

          tf_camera_params_read_mutable (f_true, ctrue);


      fclose(f_cpar);
      free(f_cpar_in_fname);
      /* Show the calibrated camera parameters: */
       fprintf(stderr, "Given camera parameters:\n");
      tf_camera_params_print (cpar, stderr);


      char *frame_in_fname = jsprintf("%s/frame.pgm", frame_data_dir);
      float_image_t *img = float_image_read_pnm_named
        ( frame_in_fname, FALSE, ttm_PNM_GAMMA, ttm_PNM_BIAS, FALSE, TRUE, FALSE );
      free(frame_in_fname);

    char *wpos_fname = jsprintf("%s/world-coords.txt", data_dir);
    int32_t n_p_w;
    tf_calib_data_t * cdat = calibration_data_new();
    tf_calib_data_read_world_points (wpos_fname, &n_p_w, &(cdat->world));
