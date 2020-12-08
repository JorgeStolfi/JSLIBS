/* Routines to read target data and show placement of targets on images */
/* Last edited on 2011-05-15 00:45:42 by stolfi */

#ifndef tf_targets_H
#define tf_targets_H

#include <r2.h>
#include <r3.h>
#include <tf_camera.h>
#include <tf_calib_data.h>

typedef struct _targets_data_t {
  int ntargets;

  int *type;       /*type of each target. e.g: 0 circle, 1 square, etc ...*/
  /* Target positions and shapes in the scene: */ 
  r3_t *p_w;       /* world coords of center of each target [mm] */
  r3_t *u_w;       /* world coords of U-axis unit vector of each target [mm] */
  r3_t *v_w;       /* world coords of V-axis unit vector of each target [mm] */
  /* Target positions, shapes, colors, and confidences in the current frame: */ 
  r2_t *p_i;       /* Image coords of center of each target [pixels] */
  r2_t *u_i;       /* Image coords of U-axis unit vector of each target [pixels] */
  r2_t *v_i;       /* Image coords of V-axis unit vector of each target [pixels] */
  double *alpha;   /* Contrast of each target in image. */
  double *beta;    /* Black level of each target in image. */
  double *conf;    /* Confidence of each target in image. */
} *targets_data_t;

targets_data_t tf_targets_data_create (int size);
  /* Allocates a {_targets_data_t} record and its {p_w}, {p_i}, and {conf} vectors.
    The remaining fields are set to NULL. */

void tf_targets_data_free (targets_data_t tdat);

void tf_targets_data_copy (targets_data_t source, targets_data_t destination);
  /* Copies the {p_w}, {p_i}, and {conf} of all targets from {source} to {destination}. */

void tf_targets_data_clear (targets_data_t tdat);
  /* Sets the {p_w}, {p_i}, and {conf} of all targets in {tdat} to zero. */

void tf_targets_write_info
  ( int target, 
    r2_t pos, 
    double alpha, 
    double beta, 
    double errorSqr, 
    double confidence, 
    FILE *f
  );
  /* Writes the data of target number {i} to file {f} in a single line. */

targets_data_t tf_targets_data_read 
  ( char *p_w_fname,
    char *uv_w_fname,
    char *p_i_fname,
    char *conf_fname );
  /* Reads a {targets_info_t} structure from files:

    {p_w_fname}   - world coordinate of target center {p_w[k]},
    {tuv_w_fname} - target type {type[k]} and world U,V unit vectors {u_w[k],v_w[k]},
    {p_i_fname}   - approx image coordinates {p_i[k]} of target center in initial frame.
    {conf_fname}  - confidence {conf[k]} of target center in initial frame.

   Also sets all image shape vectors {u_i[k],v_i[k]} to zero,
   all target colors {alpha[k],beta[k]} to {1.0,0.0}.
   If {conf_fname} is NULL, also sets all confidence weights {conf[k]} to 1.0. */

void tf_read_types_and_world_shapes
  ( char *fname, 
    int *ntargets, 
    int **type,  
    r3_t **u_world_coords, 
    r3_t **v_world_coords
  );
/* Reads a file containing a number of targets {*ntargets} followed by {*ntargets}
  target types (*type)[0..*ntargets-1]}, world U-vectors {(*u_world_coords)[0..*ntargets-1]},
  and world V-vectors {(*v_world_coords)[0..*ntargets-1]},  each vector
  being three real numbers.   Allocates the vectors {*type},
  {*u_world_coords}, {*v_world_coords}. */

#endif

