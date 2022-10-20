/* Camera model and camera parameters for most uses. */
/* Last edited on 2022-10-20 05:53:10 by stolfi */

#ifndef tf_camera_params_H
#define tf_camera_params_H

#include <r2.h>
#include <stdint.h>
#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>

/* INSTANTANEOUS CAMERA PARAMETERS */

/*********************************************************
 *                     TRANSFORMATION MATRIX             *
 *                    0      1       2      3            *
 *                ------------------------------         *
 *           0   |  1.00   0.00   0.00   0.00   |        *
 *     S  =  1   |   Tx     r1     r2     r3    |        *
 *           2   |   Ty     r4     r5     r6    |        *
 *           3   |   Tz     r7     r8     r9    |        *
 *                ------------------------------         *
*********************************************************/

typedef struct tf_camera_params_t {
    double    Npx;        /* [pix]     Number of effective sensors (pixels) in camera's x direction */
    double    Npy;        /* [pix]     Number of effective sensors (pixels) in camera's y direction */
    double    dpx;        /* [mm/pix]  Effective X dimension of pixel in camera coords   */
    double    dpy;        /* [mm/pix]  Effective Y dimension of pixel in camera coords   */
    double    Cx;         /* [pix]     Z axis intercept of camera coordinate system      */
    double    Cy;         /* [pix]     Z axis intercept of camera coordinate system      */
    /* !!! We should eliminate the {sx} parameter - adjust {dpx} instead !!! */
    double    sx;         /* []        Scale factor to compensate for any error in dpx   */
    double    f;          /* [mm]      Focal length of camera    */
    double    kappa;      /* [1/mm^2]  Radial distortion coefficient    */
    r4x4_t    S;          /* world to camera transformation matrix [mm->mm]*/
} tf_camera_params_t;
/* A record that describes a mapping from world coordinates to image coordinates,
   at a given moment in time. */

tf_camera_params_t *tf_camera_params_new (void);
/* Allocates a new {tf_camera_params_t *} record, un-initialized. */

tf_camera_params_t *tf_camera_params_copy (tf_camera_params_t *cpar);
/* Allocates a new {tf_camera_params_t *} record and copies into it
  the contents of the given {cpar}. */

void tf_camera_params_print (tf_camera_params_t *cpar, FILE *ferr);
/* This routine prints all camera parameters in {cpar} to file {ferr},
  in human-readable format, with "//" in front of every line. */

void tf_camera_matrix_print(r4x4_t *S, FILE *ferr);
/* This routine prints the matrix {S} to {ferr},
  in human-readable format, with "//" in front of every line. */

void tf_camera_matrix_print_rotation (r4x4_t *S, FILE *ferr);
/* This routine prints the rotation submatrix of matrix {S} to {ferr},
  in human-readable format, with "//" in front of every line. */

void tf_camera_params_write (FILE *wr, int32_t index, tf_camera_params_t *cpar);
/* Writes the {index} and all camera parameters from {cpar}
  to file {wr}, in a single line.

   The fields are
     "{index} 
      {S11} {S12} {S13}  
      {S21} {S22} {S23}  
      {S31} {S32} {S33}  
      {S10} {S20} {S30}
      {f} {kappa} {sx}
      {Npx} {Npy} {dpx} {dpy} {Cx} {Cy}". */

int32_t tf_camera_params_read (FILE *rd, tf_camera_params_t *cpar);
/* Reads frame {index} and all camera parameters from file {rd},
   in the format used by {tf_camera_params_write},
   and stores them into the {cpar} structure. Returns the frame index. */

void tf_camera_params_write_mutable (FILE *wr, int32_t index, tf_camera_params_t *cpar);
/* Writes the frame {index} and the mutable camera parameters from {cpar}
  ({S,f,kappa, sx}) to file {wr}, in a single line.

  The fields are
     "{index} 
      {S11} {S12} {S13}  
      {S21} {S22} {S23}  
      {S31} {S32} {S33}  
      {S10} {S20} {S30}
      {f} {kappa} {sx}".
   Ignores the fields {Npx} {Npy} {dpx} {dpy} {Cx} {Cy} of {cpar}. */

int32_t tf_camera_params_read_mutable (FILE *rd, tf_camera_params_t *cpar);
/* Reads the frame {index} and the mutable camera parameters from file {rd},
   in the format used by {tf_camera_params_write_mutable},
   and stores the latter in the {cpar} structure. Returns the frame index.
   Does nor change the fields {Npx} {Npy} {dpx} {dpy} {Cx} {Cy} of {cpar}. */

void tf_camera_params_print_changes (r3_t *Do, r3x3_t *DR, double Df, double Dkappa, FILE *ferr);
/* This routine print changes in the camera position, camera orientation matrix, focal distance
  and radial distortion. */

void tf_camera_params_free (tf_camera_params_t *cpar);
/*deallocates a camera parameters structure*/

void tf_camera_params_write_povray (FILE *wr, tf_camera_params_t *cpar, bool_t flip);
/* Writes to file {wr} the POV-Ray camera definition that is equivalent to
  the parameters {cpar} (except for radial distortion and non-central
  {Cx,Cy}). If {flip} is true,
  the images will be flipped left-to-right. */

/* ACCESSING PARAMETERS BY NUMBER */

#define tf_camera_num_parameters (15)
/* Number of indexable parameters in a {tf_camera_params_t}. */

double tf_camera_params_get_value_from_index (tf_camera_params_t *cpar, int32_t iparam, double vref);
/* Returns the parameter number {iparam} from {cpar}.
  If the parameter is an angle, return {tf_camera_adjust_angle(v,vref)}
  instead of the stored value {v}. */

void tf_camera_params_set_from_vector (tf_camera_params_t *cpar, double params[]);
/* Sets all fields of {cpar} from {params[0..n-1]} where 
  {n == tf_camera_num_parameters}. */

char *tf_camera_params_get_name_from_index (int32_t iparam);
/* The name of parameter number {iparam} in a {tf_camera_params_t *}. */

/* COORDINATE CONVERSION */

r3_t tf_world_coords_to_camera_coords (tf_camera_params_t *cpar, r3_t pw);
/* This routine takes the position of a point in world coordinates [mm]  and determines its position 
in camera coordinates [mm].  */

r2_t tf_camera_coords_to_und_sensor_coords (tf_camera_params_t *cpar, r3_t pc);
/* convert a point {pc} from camera coordinates [mm] to undistorted sensor plane coordinates [mm] */

r2_t tf_und_sensor_coords_to_dis_sensor_coords (tf_camera_params_t *cpar, r2_t pu);
/* convert a point {pu} from undistorted sensor plane coordinates [mm] to distorted 
sensor plane coordinates [mm]*/

r2_t tf_sensor_coords_to_image_coords (tf_camera_params_t *cpar, r2_t pd);
/* convert a point {pd} from sensor plane coordinates [mm] to image coordinates [pixel] */

r2_t tf_image_coords_to_sensor_coords (tf_camera_params_t *cpar, r2_t pi);
/* convert a point {pi} from image coordinates [pixels] into sensor plane coordinates [mm]*/

r2_t tf_dis_sensor_coords_to_und_sensor_coords (tf_camera_params_t *cpar, r2_t pd);
/* convert a point {pd} from distorted sensor plane coordinates [mm] into undistorted sensor plane coordinates [mm]*/

r2_t tf_world_coords_to_image_coords (tf_camera_params_t *cpar, r3_t pw);
/* This routine takes a point {pw} in world coordinates [mm] and determines the image coordinates 
of its projection [pixels]. */

r2_t tf_camera_apply_kappa(r2_t *p, double kappa);

/* PARAMETER MANIPULATION */

double tf_camera_adjust_angle (double ang, double aref);
/* Adds to {ang} an integer multiple of {2*PI} so that it 
  lies in {[aref-PI _ aref+PI]}.  If {aref} is {NAN}, 
  however, returns {ang} unchanged. */

void tf_camera_get_mark_position_and_shape 
  ( tf_camera_params_t *cpar,
    r3_t p_w,
    r3_t u_w,
    r3_t v_w,
    r2_t *p_i,
    r2_t *u_i,
    r2_t *v_i
  );
/* Given the world point position {pw}, two r3 vectors {uw} and {vw} ({uw} and
 {vw} describes the orientation of the point in space), the camera constants
 {cpar} (in our case cannon optura) and {cpar} which contains S
 (transformation matrix), f (focal distance) and kappa (distortion), returns
 the shape of the mark and a r2 point which is used as initial guess to find
 this mark.*/

void tf_camera_matrix_split(r4x4_t *S, r3x3_t *R, r3_t *o);
/* Extracts from the world-to-camera matrix {S} (whose first row
  should be [1,0,0,0]) the world coordinates {o} of the camera, and
  the camera rotation matrix {R} (whose columns are the world
  coordinates of the {x,y,z} camera axis directions. */

void tf_camera_matrix_assemble(r3x3_t *R, r3_t *o, r4x4_t *S);
/* Combines the world coordinates {o} of the camera, and the camera
  rotation matrix {R} (as returned by {}) into the world-to-camera
  matrix {S} (whose first row will be should be [1,0,0,0]). */

bool_t tf_camera_point_is_inside_image(r2_t p_i, tf_camera_params_t *cpar);
/* TRUE iff the point {p_i} is inside the domain of the image defined by {cpar}. */

void tf_camera_extrapolate 
  ( r4x4_t *S_2, double f_2, double kappa_2, 
    r4x4_t *S_1, double f_1, double kappa_1, 
    r4x4_t *S_0, double *f_0, double *kappa_0
  );
/* Estimates the camera parameters {*S_0,*f_0,*kappa_0}
  of a frame, by extrapolating the camera parameters
  {S_2,f_2,kappa_2} and {S_1,f_1,kappa_1} of the
  two previous frames. 
  !!! Should take {cpar} records instead of matrices etc. !!! */

double tf_camera_maximum_safe_kappa ( tf_camera_params_t *cpar );
  /* Computes the maximum safe value of {kappa} from 
     the sensor size data. */

/* CAMERA MODEL ERROR METRICS */

r2_t tf_camera_compute_image_error(tf_camera_params_t *cpar, r3_t p_w, r2_t p_i);
/* Vector difference between the image coordinates {p_i_pre} predicted by the model {cpar}
   from the world coordinates {p_w}, and the observed image coordinates {p_i}. */

void tf_camera_compute_all_image_errors
  ( tf_camera_params_t *cpar,
    int32_t nmarks,
    r3_t p_w[],
    r2_t p_i[],
    r2_t e_i[] );
/* Sets {e_i[k]} to {tf_camera_compute_image_error(cpar,p_w[k],p_i[k])} 
  for all {k} in {0..nmarks-1}; */

double tf_camera_compute_world_error_sqr(tf_camera_params_t *cpar, r3_t p_w, r2_t p_i);
/* Computes the quared distance of closest approach (i.e. 3D error) between the point
  with world coordinates {p_w} and the line of sight formed by back projecting the 
   measured image point {p_i} out through the camera model. */

interval_t tf_camera_adjust_angle_range (interval_t ang, double aref);
/* Shifts the range {ang} by an integer multiple of {2*PI} so that its
  midpoint is as close as possible to {aref}.  If {aref} is {NAN} 
  or {ang} is infinite, however, returns {ang} unchanged. */

void tf_camera_matrix_to_euler_angles ( r4x4_t *S, r3_t *R );
/* This procedure extracts the Euler angles from the rotation submatrix
  of the world-to-camera projective transformation matrix {S}. */

void tf_camera_matrix_from_euler_angles ( r3_t *R, r4x4_t *S );
/* This procedure converts the Euler angles to a 3x3 orthonormal
  matrix and stores it into the rotation submatrix of the
  world-to-camera projective transformation matrix {S}. */

void tf_camera_matrix_inverse_from_v_w_and_R (r3_t v_w, r4x4_t *S);
/**/

r3_t tf_camera_matrix_to_v_w (r4x4_t *S);
/**/

/* PARAMETER REMAPPING TOOLS */

double tf_camera_stretch_param (double x, interval_t *xr);
/* Maps {x} from the range {*xr} to the range {-INF,+INF}. */

double tf_camera_squeeze_param (double e, interval_t *er);
/* Maps {e} from the range {-INF,+INF} to the range {*er}. */

double tf_camera_safe_log (double x);
/* Logarithm of {x} clipped to avoid infinities and NaNs. */

double tf_camera_safe_exp (double e);
/* Exponential of {e} clipped to avoid infinities and NaNs. */

interval_t tf_camera_interval_safe_log (interval_t *xr);
/* Interval version of {tf_camera_safe_log}. */

interval_t tf_camera_interval_safe_exp (interval_t *er);
/* Interval version of {tf_camera_safe_exp}. */

#endif
