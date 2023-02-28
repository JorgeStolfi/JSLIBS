/* See {tf_calib_data.h}.  */
/* Last edited on 2023-02-25 16:11:55 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>

#include <jsfile.h>
#include <affirm.h>
#include <r2.h>
#include <r3.h>
#include <fget.h>
#include <rn.h>

#include <tf_calib_data.h>

tf_calib_data_t *tf_calib_data_new (void)
  {
    tf_calib_data_t *cdat = notnull(malloc(sizeof(struct tf_calib_data_t)), "no mem");
    return cdat;
  }
  
/* ---------------------------------------------------------------------- */
/* READING */  

void tf_calib_data_read_world_points (char *fname, int32_t *np, r3_t **pp)
  {
    FILE *rd = open_read(fname, TRUE);
    int32_t n = fget_int32(rd); fget_eol(rd);
    demand((n > 0) && (n <= MAX_POINTS), "invalid number of points in file");
    r3_t *p = notnull(malloc(n * sizeof(r3_t)),"no mem");
    int32_t i;
    for (i = 0; i < n; i++)
      { r3_t *pi = &(p[i]);
        pi->c[0] = fget_double(rd);
        pi->c[1] = fget_double(rd);
        pi->c[2] = fget_double(rd);
        fget_comment_or_eol(rd, '#');
      }
    fclose(rd);
    (*pp) = p; 
    (*np) = n;
  }

void tf_calib_data_read_image_points (char *fname, int32_t *np, r2_t **qq)
  {
    FILE *rd = open_read(fname, TRUE);
    int32_t n = fget_int32(rd); fget_eol(rd);
    demand((n > 0) && (n <= MAX_POINTS), "invalid number of points in file");
    r2_t *q = notnull(malloc(n * sizeof(r2_t)),"no mem");
    int32_t i;
    for (i = 0; i < n; i++)
      { r2_t *qi = &(q[i]);
        qi->c[0] = fget_double(rd);
        qi->c[1] = fget_double(rd);
        fget_comment_or_eol(rd, '#');
      }
    fclose(rd);
    (*qq) = q; 
    (*np) = n;
  }

void tf_calib_data_read_weights (char *fname, int32_t *np, double **ww)
  {
    FILE *rd = open_read(fname, TRUE);
    int32_t n = fget_int32(rd); fget_eol(rd);
    demand((n > 0) && (n <= MAX_POINTS), "invalid number of points in file");
    double *w = notnull(malloc(n * sizeof(double)),"no mem");
    int32_t i;
    for (i = 0; i < n; i++)
      { w[i] = fget_double(rd);
        fget_comment_or_eol(rd, '#');
      }
    fclose(rd);
    (*ww) = w; 
    (*np) = n;
  }

tf_calib_data_t *tf_calib_data_read (char *world_fname, char *image_fname, char *weight_fname)
  {
    int32_t n_world, n_image, n_weight;
    tf_calib_data_t *cdat = tf_calib_data_new();
    tf_calib_data_read_world_points (world_fname, &n_world, &(cdat->world));
    tf_calib_data_read_image_points (image_fname, &n_image, &(cdat->image));
    demand(n_image == n_world, "inconsistent world and image point counts");
    if (weight_fname != NULL)
      { tf_calib_data_read_weights (weight_fname, &n_weight, &(cdat->weight));
        demand(n_weight == n_world, "inconsistent world point and weight counts");
      }
    else
      { /* Assume unit weight for all points */
        n_weight = n_world;
        cdat->weight = rn_alloc(n_weight);
        tf_calib_data_set_weights_uniform(cdat);
      }
    cdat->np = n_world;
    return cdat; 
  }
  
/* ---------------------------------------------------------------------- */
/* WRITING TO FILES */  

void tf_calib_data_write_world_points (char *fname, int32_t np, r3_t p[])
  {
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "%d\n", np);
    int32_t i;
    for (i = 0; i < np; i++)
      { r3_t pi = p[i];
        fprintf(wr, " %14.6f %14.6f %14.6f  # %03d\n", pi.c[0], pi.c[1], pi.c[2], i);
      }
    fclose(wr);
  }

void tf_calib_data_write_image_points (char *fname, int32_t np, r2_t q[])
  {
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "%d\n", np);
    int32_t i;
    for (i = 0; i < np; i++)
      { r2_t qi = q[i];
        fprintf(wr, " %9.3f %9.3f  # %03d\n", qi.c[0], qi.c[1], i);
      }
    fclose(wr);
  }

void tf_calib_data_write_weights (char *fname, int32_t np, double w[])
  {
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "%d\n", np);
    int32_t i;
    for (i = 0; i < np; i++)
      { fprintf(wr, " %11.9f  # %03d\n", w[i], i); }
    fclose(wr);
  }
  
/* ---------------------------------------------------------------------- */
/* PRINTING */  

void tf_calib_data_print (FILE *wr, tf_calib_data_t *cdat)
  {
    fprintf(wr, "--------- calibration data ---------------------------------\n");
    int32_t i;
    for (i = 0; i < cdat->np; i++)
      { r3_t pi = cdat->world[i];
        r2_t qi = cdat->image[i];
        double wi = cdat->weight[i];
        fprintf(wr, "[%03d]", i);
        fprintf(wr, "  world = ( %14.6f %14.6f %14.6f )", pi.c[0], pi.c[1], pi.c[2]);
        fprintf(wr, "  image = ( %9.3f %9.3f )", qi.c[0], qi.c[1]);
        fprintf(wr, "  weight = %11.9f\n", wi);
      }
    fprintf(wr, "------------------------------------------------------------\n");
  }

void tf_calib_data_print_world_points (FILE *wr, int32_t np, r3_t p[])
  {
    fprintf(wr, "--------- world points -------------------------------------\n");
    int32_t i;
    for (i = 0; i < np; i++)
      { r3_t pi = p[i];
        fprintf(wr, "[%03d]", i);
        fprintf(wr, "  p = ( %14.6f %14.6f %14.6f )\n", pi.c[0], pi.c[1], pi.c[2]);
      }
    fprintf(wr, "------------------------------------------------------------\n");
  }

void tf_calib_data_print_image_points (FILE *wr, int32_t np, r2_t q[])
  {
    fprintf(wr, "--------- image points ------------------------------------\n");
    int32_t i;
    for (i = 0; i < np; i++)
      { r2_t qi = q[i];
        fprintf(wr, "[%03d]", i);
        fprintf(wr, "  q = ( %9.3f %9.3f )\n", qi.c[0], qi.c[1]);
      }
    fprintf(wr, "------------------------------------------------------------\n");
  }

void tf_calib_data_print_weights (FILE *wr, int32_t np, double w[])
  {
    fprintf(wr, "--------- point weights ------------------------------------\n");
    int32_t i;
    for (i = 0; i < np; i++)
      { fprintf(wr, "[%03d]", i);
        fprintf(wr, "  w = %11.9f\n", w[i]);
      }
    fprintf(wr, "------------------------------------------------------------\n");
  }

void tf_calib_data_set_weights(tf_calib_data_t *cdat, double w[])
  { int32_t k;
    for (k = 0; k < cdat->np; k++) { cdat->weight[k] = w[k]; }
  }

void tf_calib_data_set_weights_uniform(tf_calib_data_t * cdat)
  { int32_t k;
    for (k = 0; k < cdat->np; k++) { cdat->weight[k] = 1.0; }
  }

void tf_calib_data_free (tf_calib_data_t * cdat)
{
    free(cdat->world);
    free(cdat->image);
    free(cdat->weight);
    free(cdat);
}


