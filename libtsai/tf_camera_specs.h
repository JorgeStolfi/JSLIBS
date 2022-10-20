/* Specifying ranges for camera parameters. */
/* Last edited on 2022-10-20 05:52:01 by stolfi */

#ifndef tf_camera_specs_H
#define tf_camera_specs_H

#include <r2.h>
#include <stdint.h>
#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>
#include <interval.h>
#include <values.h>

#include <tf_camera.h>
#include <tf_camera_specs.h>

/* CAMERA PARAMETER RANGES */

typedef struct tf_camera_specs_t 
  { /* !!! The parameters Npx and Npy should be plain integers, not intervals. !!! */
    interval_t    Npx;        /* [pix]     Number of effective sensors (pixels) in camera's x direction */
    interval_t    Npy;        /* [pix]     Number of effective sensors (pixels) in camera's y direction */
    interval_t    dpx;        /* [mm/pix]  Effective X dimension of pixel in camera coords   */
    interval_t    dpy;        /* [mm/pix]  Effective Y dimension of pixel in camera coords   */
    interval_t    Cx;         /* [pix]     Z axis intercept of camera coordinate system      */
    interval_t    Cy;         /* [pix]     Z axis intercept of camera coordinate system      */
    interval_t    sx;         /* []        Scale factor to compensate for any error in dpx   */
    interval_t    f;          /* [mm]      Focal length of camera    */
    interval_t    kappa;      /* [1/mm^2]  Radial distortion coefficient    */
    /* The *camera's* position in world coords. */
    interval_t    v_w[3];     /* [mm]      Camera X,Y,Z coordinate of world's origin.  */
    /* The camera's orientation with respect to the world's system: */
    interval_t    R[3];       /* [radians] Euler angles.  */
  } tf_camera_specs_t;
/* A record that describes the possible {camera_parameter_t} settings
  for a particular camera.  If an interval is trivial, it means that
  the corresponding parameter has that fixed value.

  If the interval is not trivial, some functions assume that the
  parameter has an /a priori/ propability distribution that is uniform
  within the interval, and zero outside.  Other functions assume that
  is has a Gaussian distribution such that the given interval is {avg
  +/- 3*dev} where {avg} is the parameter's expected value and {dev}
  is its standard deviation.  However, for the {f} parameter, the
  Gaussian distribution is over {log(f)} rather than {f}.

  The {Rx,Ry,Rz} intervals should be either all trivial (for a camera
  with fixed orientation) or all non-trivial. Also, either all three
  are {[-INF,+INF]}, or all three have finite endpoints.
  
  The {Vx,Vy,Vz} intervals should be either all trivial (for a camera
  with fixed position) or all non-trivial. Also, either all three
  are {[-INF,+INF]}, or all three have finite endpoints.
*/

tf_camera_specs_t *tf_camera_specs_new (void);
/* Allocates a new camera specs record, uninitialized. */

tf_camera_specs_t *tf_camera_specs_copy (tf_camera_specs_t *cspec);
/* Creates a fresh copy of the camera specs {cspec}. */

void tf_camera_specs_print (FILE *wr, tf_camera_specs_t *cspec);
/* This routine prints all camera specs in {cspec} to file {wr},
  in human-readable format, with "//" in front of every line. */

void tf_camera_specs_write (FILE *wr, tf_camera_specs_t *cspec);
/* Writes {cspec} to file {wr}, in a format that can be read back by 
  {tf_camera_specs_read}. */

tf_camera_specs_t *tf_camera_specs_read (FILE *rd);
/* Reads a {tf_camera_specs_t} record from file {rd}. */

interval_t tf_camera_specs_read_range (FILE *rd, char *name);
/* Reads from file {rd} the name (a string of non-blanks) and 
  value range (two numbers) for a camera parameter. The fields
  must be separated by spaces.  The name must match {name}.
  !!! Details? !!! */

interval_t tf_camera_specs_get_param_range (tf_camera_specs_t *cspec, int32_t iparam, double vref);
/* Obtains the range of parameter {iparam} from {cspec}.  The parameter
  numbering and conventions are the same as those of {tf_camera_params_get_value_from_index}.
  For the {R} parameters, the resulting interval is adjusted to contain {vref}
  (if it is not NAN). */

void tf_camera_specs_set_param_range (tf_camera_specs_t *cspec, int32_t iparam, interval_t *range);
/* Sets the range of parameter {iparam} into {cspec} as {range}.  The parameter
  numbering and conventions are the same as those of {tf_camera_params_get_value_from_index}. */

tf_camera_specs_t *tf_camera_specs_for_canon_optura (void);
/*creates a camera specs structure for the canon optura 40*/

tf_camera_specs_t *tf_camera_specs_for_sony_dv40 (void);
/*creates a camera specs structure for the sony dv40*/

tf_camera_specs_t *tf_camera_specs_for_povray_svga (void);
/*creates a camera specs structure for POV-Ray-generated videos with 640x480. */

tf_camera_specs_t *tf_camera_specs_for_povray_svga_webcam (void);
/*creates a camera specs structure for POV-Ray-generated videos with 640x480.
  This structure try to simulate a webcam with fixed focal distance */

tf_camera_specs_t *tf_camera_specs_for_povray_hvga (void);
/*creates a camera specs structure for POV-Ray-generated videos with 320x240. */

tf_camera_specs_t *tf_camera_specs_for_povray_hvga_distorted (void);
/*creates a camera specs structure for POV-Ray-generated videos with 320x240
  with radial distorion generate by pnmdistrad program. */

void tf_camera_specs_get_mean_params(tf_camera_specs_t *cspec, tf_camera_params_t *cpar);
/* Stores into {cpar} the `average' or `default' values of 
   the camera parameters, compatible with {cspec}. */

tf_camera_params_t *tf_camera_specs_get_new_mean_params(tf_camera_specs_t *cspec);
/* Creates a {tf_camera_params_t *} record and fills it with `average' or `default' values of 
   the camera parameters, compatible with {cspec}. 
   !!! Do we need this routine? !!! */

#endif
