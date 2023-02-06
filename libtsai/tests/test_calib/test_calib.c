/* Last edited on 2023-02-04 09:04:43 by stolfi */

#define PROG_NAME "test_calib"
#define PROG_DESC "tests the Tsai camera calibration algorithm"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <r2x2.h>
#include <r2.h>
#include <assert.h>
#include <jsfile.h>
#include <affirm.h>

#include <tf_camera.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>
#include <tf_camera_plot.h>
#include <tf_calib_guess1.h>
#include <tf_calib_guess2.h>
#include <tf_calib_refine2.h>

tf_camera_params_t *read_camera_parameters(tf_camera_specs_t *cspec, char *in_dir, char *name);
tf_calib_data_t *read_calibration_data(char *in_dir, char *tag);
void plot_cpars
  ( tf_camera_specs_t *cspec,
    tf_calib_data_t *cdat,
    tf_camera_params_t *cal_cpar,
    tf_camera_params_t *ref_cpar,
    char *out_dir
  );
  
int32_t main (int32_t argc, char *argv[])
{
  /* Arguments: */
  char *camera_name = argv[1];
  char *in_dir = argv[2];
  char *out_dir = argv[3];

  int32_t which_calib = atoi(argv[4]);             /* Version of {tf_calib} to use. */
  int32_t nl_optimization = (atoi(argv[5]) != 0);  /* 0 to skip the NL optimization. */
  int32_t has_ref_calib = atoi(argv[6]);           /* 1 if has true calibration. */
  
  fprintf(stderr, "Starting...\n");
  
  /* Get the fixed camera parameters: */
  tf_camera_specs_t *cspec;
  if (strcmp(camera_name,"optura") == 0) {
    cspec = tf_camera_specs_for_canon_optura ();
  }
  else if (strcmp(camera_name,"svga") == 0) {
    cspec = tf_camera_specs_for_povray_svga ();
  }
  else if (strcmp(camera_name,"hvga") == 0) {
    cspec = tf_camera_specs_for_povray_hvga ();
  }
  else {
    fprintf(stderr, "error: define the camera specifications\n");
    exit(1);
  }
  fprintf(stderr, "Generic camera parameters:\n");
  tf_camera_specs_print(stderr, cspec);

  /* A priori uncertainty on mark image positions (for error computation): */
  r2_t q_dev = (r2_t){{ 5.0, 5.0 }};  /* Should be a parameter... */

  /* Shall we include the {cpar} Kalman error terms in the optimization? */
  bool_t use_cpar_dev_terms = FALSE;

  /*Output data directory*/
  mkdir(out_dir, 0777);

  /* Read the calibration data: */
  tf_calib_data_t *cdat = read_calibration_data(in_dir, "p");
  fprintf(stderr, "calibration data:\n");
  tf_calib_data_print(stderr, cdat); 
  fprintf(stderr, "-----------------------------\n");

  /* Read the true {cpar}: */
  fprintf(stderr, "True parameters:\n");
  tf_camera_params_t *ref_cpar = NULL;
  if (has_ref_calib > 0)
    { ref_cpar = read_camera_parameters(cspec, in_dir, "true");
      tf_write_cpar_and_errors("true", cdat, cspec, q_dev, FALSE, ref_cpar, out_dir); 
      tf_camera_params_print (ref_cpar, stderr);
    }

  /* Compute the calibrated {cpar}: */
  fprintf(stderr, "Calibrated parameters:\n");
  tf_camera_params_t *cal_cpar = tf_camera_specs_get_new_mean_params(cspec); /* Calibrated. */
  if (which_calib == 1) 
    { demand(FALSE, "calibration 1 not available");
      /* tf_calib1_tsai (cdat, cspec, cpar, TRUE, nl_optimization); */ 
    }
  else if (which_calib == 2)
    { tf_calib_guess2_initial_camera_parameters
        (cdat, cspec, cal_cpar, FALSE, TRUE);
      if (nl_optimization) {
        tf_calib_generic
          ( cdat, cspec, q_dev, use_cpar_dev_terms, cal_cpar,
            tf_calib_refine2_gather_optimization_params,
            tf_calib_refine2_scatter_optimization_params,
            TRUE);
      }
    }
  else
    { demand(FALSE, "invalid calibration"); }
  tf_camera_params_print (cal_cpar, stderr);
  tf_write_cpar_and_errors("calib", cdat, cspec, q_dev, FALSE, cal_cpar, out_dir); 

  plot_cpars(cspec, cdat, cal_cpar, ref_cpar, out_dir);
  
  tf_camera_params_free (cal_cpar);
  if (ref_cpar != NULL) 
    { tf_camera_params_free (ref_cpar); }

  fprintf(stderr, "Done.\n");
  return 0;
}

void plot_cpars
  ( tf_camera_specs_t *cspec,
    tf_calib_data_t *cdat,
    tf_camera_params_t *cal_cpar,
    tf_camera_params_t *ref_cpar,
    char *out_dir
  )
  {
    /* Camera image dimensions and aspect ratio: */
    double Nx = cal_cpar->Npx;
    double Ny = cal_cpar->Npy;
    double aspect = Nx/Ny;

    /* Postscript figure sizeand margin, in pt: */
    double hSize = 400.0*sqrt(aspect);
    double vSize = 400.0/sqrt(aspect);
    double hMrg = 4;
    double vMrg = 3;

    /* Plot the figures: */
    int32_t np = cdat->np;
    r2_t *q = cdat->image;
    r3_t *p = cdat->world;

    tf_camera_params_t *fig1[2] = { ref_cpar, cal_cpar };

    tf_plot_cameras(out_dir, "true_calib", np, p, 2, fig1, q, hSize, vSize, hMrg, vMrg, Nx, Ny, 1.0, 0);
  }
 
tf_camera_params_t *read_camera_parameters(tf_camera_specs_t *cspec, char *in_dir, char *name)
  { char *cpar_fname = NULL;
    asprintf(&cpar_fname, "%s/%s.cpar", in_dir, name);
    FILE *cpar_file = open_read(cpar_fname, TRUE);
    tf_camera_params_t *cpar = tf_camera_specs_get_new_mean_params(cspec);
    tf_camera_params_read(cpar_file, cpar);
    if (cpar_file != stdin) { fclose(cpar_file); }
    free(cpar_fname);
    return cpar;
  }

tf_calib_data_t *read_calibration_data(char *in_dir, char *tag)
  {  char *world_coords_fname = NULL;
    asprintf(&world_coords_fname, "%s/%s_w.txt", in_dir, tag);
    char *image_coords_fname = NULL;
    asprintf(&image_coords_fname, "%s/%s_i.txt", in_dir, tag);
    char *point_weights_fname = NULL;
    asprintf(&point_weights_fname, "%s/%s_wgt.txt", in_dir, tag);
    tf_calib_data_t *cdat = tf_calib_data_read
      ( world_coords_fname, image_coords_fname, point_weights_fname );
    free(image_coords_fname);
    free(world_coords_fname);
    free(point_weights_fname);
    return cdat;
  }
 
