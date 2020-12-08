
#include <tf_calib.h>
#include <tf_camera.h>
#include <tf_errors.h>
#include <tf_math.h>

#define DEBUG_CALIB_ERRORS 1

/* Code originaly written by Reg Willson (18-May-95). */

void tf_calib_summarize_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr)
{
  tf_camera_params_print (cpar, ferr);

  /* calculate the distorted image plane error statistics for the data set */
  tf_calib_summarize_image_errors (cdat, cpar, ferr);

  /* calculate the undistorted image plane error statistics for the data set */
  tf_calib_summarize_world_errors (cdat, cpar, ferr);
}

void tf_calib_summarize_image_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr)
{
  r2_t avg, rms, max;
  if (cdat->np < 1) {
    avg = rms = max = (r2_t) {{ 0, 0 }};
  }
  else {

  int i;
  r2_t sum_e = (r2_t) {{ 0, 0 }};
  r2_t sum_e2 = (r2_t) {{ 0, 0 }};
  double max_t2 = 0.0;
  max = (r2_t) {{ 0, 0 }};

  for (i = 0; i < cdat->np; i++) {
    /* Compute the error vector {e} in image coordinates: */
    r3_t p_w = cdat->world[i];
    r2_t p_i = cdat->image[i];
    r2_t e = tf_camera_compute_image_error(cpar, p_w, p_i);
    if (DEBUG_CALIB_ERRORS) {
      fprintf(stderr, "  mark %3d:", i);
      fprintf(stderr, "  p_i = ( %8.3lf %8.3lf )", p_i.c[0], p_i.c[1]);
      fprintf(stderr, "  error = ( %8.3lf %8.3lf )\n", e.c[0], e.c[1]);
    }
    r2_add(&e, &sum_e, &sum_e);
    r2_t e2 = (r2_t) {{ SQR (e.c[0]), SQR (e.c[1]) }};
    r2_add(&e2, &sum_e2, &sum_e2);
    double t2 = e2.c[0] + e2.c[1];
    if (t2 > max_t2) { max_t2 = t2; max = e; }
  }

  avg = (r2_t){{ sum_e.c[0]/cdat->np,  sum_e.c[1]/cdat->np }};
  rms = (r2_t){{ sqrt(sum_e2.c[0]/cdat->np),  sqrt(sum_e2.c[1]/cdat->np) }};

  }
  fprintf (ferr, "## Image position errors [pix]:");
  fprintf (ferr, "  avg = ( %8.3lf %8.3lf )", avg.c[0], avg.c[1]);
  fprintf (ferr, "  rms = ( %8.3lf %8.3lf )", rms.c[0], rms.c[1]);
  fprintf (ferr, "  max = ( %8.3lf %8.3lf ) = %8.3lf\n", max.c[0], max.c[1], r2_norm(&max));
}

void tf_calib_summarize_world_errors (tf_calib_data_t * cdat, tf_camera_params_t *cpar, FILE *ferr)
{
  double rms, max;
  if (cdat->np < 1) {
    rms = max = 0;
  }
  else {
    double sum_e2 = 0.0, max_e2 = 0.0;

    int i;
    for (i = 0; i < cdat->np; i++) {
      r3_t p_w = cdat->world[i];
      r2_t p_i = cdat->image[i];
      double e2 = tf_camera_compute_world_error_sqr(cpar, p_w, p_i);
      sum_e2 += e2;
      max_e2 = MAX(max_e2, e2);
    }

    rms = sqrt(sum_e2 / cdat->np);
    max = sqrt(max_e2);
  }
  fprintf (ferr, "## Object space error:");
  fprintf (ferr, "  rms = %10.3lf  max = %10.3lf [mm]\n", rms, max);
}
