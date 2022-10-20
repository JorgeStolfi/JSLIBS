/* General-purpose data and functions for calibration routines. */
/* Last edited on 2022-10-20 05:54:53 by stolfi */

#ifndef tf_calib_data_H
#define tf_calib_data_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <r2.h>
#include <r3.h>

/* Max number of world-image data pairs: */
#define MAX_POINTS 500

/* ---------------------------------------------------------------------- */
/* CALIBRATION DATA STRUCTURE */

/* A data set for weighted camera calibration: */
typedef struct tf_calib_data_t 
  {
    int32_t    np;       /* Number of calibration points. */
    r3_t   *world;   /* World coordinates of each calibration point [mm]. */
    r2_t   *image;   /* Image coordinates of each calibration point [pixels]. */
    double *weight;  /* Confidence weight of (world,image) data pair. */
  } tf_calib_data_t;

/* ---------------------------------------------------------------------- */
/* CALIBRATION DATA I/O */

tf_calib_data_t *tf_calib_data_new (void);
/* Allocates a new {tf_calib_data_t} record. Initializes {.np} with zero
  and the data tables with NULL. */ 

tf_calib_data_t *tf_calib_data_read (char *world_fname, char *image_fname, char *weight_fname);
/* Reads the world coordinates, image coordinates, and confidence
  weights of the marks from files {world_fname}, {image_fname}, and
  {weight_name}, respectively.  If {weight_name} is NULL, sets all weights
  to 1. 

  In each file, the first line must be the number of marks {np}.
  The following lines must contain the mark data, one mark per line,
  consisting of three, two, and one real numbers, respectively.
  The line may end optionally in a '#' followed by an arbitrary comment. */

void tf_calib_data_free (tf_calib_data_t *cdat);
/* Reclaims the record {*cdat} and its internal vectors. */

void tf_calib_data_print (FILE *wr, tf_calib_data_t *cdat);
/* Prints all the parameters of the calibration data structure. */

/* ---------------------------------------------------------------------- */
/* I/O OF SEPARATE POINT LISTS */


void tf_calib_data_read_world_points (char *fname, int32_t *np, r3_t **pp);
/* Reads a file containing a number of marks {*np} followed by {*np}
  world coordinates {(*pp)[0..*np-1]}, each being three real numbers. 
  Allocates the vector {*pp}. */

void tf_calib_data_write_world_points (char *fname, int32_t np, r3_t p[]);
/* Writes {p[0..np-1]} to file {fname}, in a format that can
  be read back by {tf_calib_data_read_world_points}. */

void tf_calib_data_print_world_points (FILE *wr, int32_t np, r3_t p[]);
/* Prints {p[0..np-1]} to {wr} in a human-readable format. */

void tf_calib_data_read_image_points (char *fname, int32_t *np, r2_t **qq);
/* Reads a file containing a number of marks {*np} followed by {*np}
  image coordinates {(*qq)[0..*np-1]}, each being two real numbers. 
  Allocates the vector {*qq}. */

void tf_calib_data_write_image_points (char *fname, int32_t np, r2_t q[]);
/* Writes {q[0..np-1]} to file {fname}, in a format that can
  be read back by {tf_calib_data_read_image_points}. */

void tf_calib_data_print_image_points (FILE *wr, int32_t np, r2_t q[]);
/* Prints {q[0..np-1]} to {wr} in a human-readable format. */


void tf_calib_data_read_weights (char *fname, int32_t *np, double **ww);
/* Reads a file containing a number of marks {*np} followed by {*np}
  confidence weights {(*ww)[0..*np-1]}, each being one real number. 
  Allocates the vector {*ww}. */

void tf_calib_data_write_weights (char *fname, int32_t np, double w[]);
/* Writes {w[0..np-1]} to file {fname}, in a format that can
  be read back by {tf_calib_data_read_weights}. */

void tf_calib_data_print_weights (FILE *wr, int32_t np, double w[]);
/* Prints {w[0..np-1]} to {wr} in a human-readable format. */

/* ---------------------------------------------------------------------- */
/* WEIGHT MANIPULATION */

void tf_calib_data_set_weights (tf_calib_data_t *cdat, double w[]);
/* Sets {cdat->weight[k] = w[k]} for all {k} in {0..cdat->np-1}. */

void tf_calib_data_set_weights_uniform (tf_calib_data_t *cdat);
/* Sets {cdat->weight[k] = 1} for all {k}. */

void tf_calib_data_skip_comment(FILE *rd);
/* Skips spaces until the next EOF or non-space character;
   if that character is '#', skips it and everthing else
   on that line (except the final '\n'). */

#endif
