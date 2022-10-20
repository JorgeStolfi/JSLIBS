/* Last edited on 2022-10-20 05:52:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <affirm.h>
#include <pswr.h>
#include <jsfile.h>

#include <tf_camera.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>
#include <tf_camera_plot.h>

void tf_plot_cameras
  ( PSStream *ps, 
    char *name, 
    int32_t np,
    r3_t p[], 
    int32_t ncpars,
    tf_camera_params_t *cpar[],
    r2_t q[], 
    double hSize,
    double vSize,
    double hMrg,
    double vMrg,
    double Nx,
    double Ny,
    double mRadius,
    int32_t mStyle
  ) 
  {
    pswr_new_canvas(ps, name);
    /* Negate all Y coordinates to get the proper orientation for the Y axis: */

    pswr_set_window(ps, 0.0, Nx,  -Ny, 0.0, hMrg, hSize+hMrg, vMrg, vSize+vMrg);
    /* Draw a frame around the image: */
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0,0);
    pswr_frame(ps);

    if (p != NULL) 
      { /* Plot computed image positions: */
        pswr_set_fill_color(ps, 1.000, 1.000, 0.400);
        pswr_set_pen(ps, 0.0,0.0,0.0, 0.20, 0,0);
        r2_t *q_aux = (r2_t *)notnull(malloc(np*sizeof(r2_t)), "no mem");
        int32_t k;
        for (k = 0; k < ncpars; k++) 
          { if (cpar[k] != NULL)
              { int32_t i;
                for (i = 0; i < np; i++)
                  { r3_t pi = p[i];
                    fprintf(stderr, "p[%d] = ( %f %f %f )\n", i, pi.c[0], pi.c[1], pi.c[2]);
                    q_aux[i] = tf_world_coords_to_image_coords(cpar[k], pi);
                    r2_t qi = q_aux[i];
                    fprintf(stderr, "q[%d] = ( %f %f )\n", i, qi.c[0], qi.c[1]);
                  }
                tf_plot_marks(ps, np, q_aux, mRadius, mStyle + k);
              }
          }
      }
    if (q != NULL) 
      { /* Plot the given image coordinates: */
        pswr_set_fill_color(ps, 0.700, 0.900, 1.000);
        pswr_set_pen(ps, 0.0,3.0,1.0, 0.20, 0,0);
        tf_plot_marks(ps, np, q, mRadius, mStyle + ncpars);
      }
  }

void tf_plot_marks
  ( PSStream *ps, 
    int32_t np,
    r2_t q[], 
    double mRadius,
    int32_t mStyle
  )
  {    
    int32_t i;
    double oRadius; /* Optical radius of a mark with unit radius: */
    for (i = 0; i < np; i++) {
      r2_t qi = q[i];
      switch (mStyle) {
      case 0:
        oRadius = 1.0000000;
        pswr_dot(ps, qi.c[0], -qi.c[1], mRadius/oRadius, TRUE, TRUE);
	break;
      case 1:
        oRadius = 0.8408964 * 0.75;
        pswr_cross(ps, qi.c[0], -qi.c[1], mRadius/oRadius, FALSE, TRUE);
	break;
      case 2:
        oRadius = 0.8408964 * 0.75;
        pswr_cross(ps, qi.c[0], -qi.c[1], mRadius/oRadius, TRUE, TRUE);
	break;
      case 3:
        oRadius = 0.9659258 * 0.83;
        pswr_asterisk(ps, qi.c[0], -qi.c[1], mRadius/oRadius, TRUE);
	break;
      case 4:
        oRadius = 0.8408964;
        pswr_square(ps, qi.c[0], -qi.c[1], mRadius/oRadius, TRUE, TRUE);
	break;
      default:
        oRadius = 0.8408964;
        pswr_diamond(ps, qi.c[0], -qi.c[1], mRadius/oRadius, mRadius/oRadius, TRUE, TRUE);
      }
    }
  }
