/* Last edited on 2023-02-25 16:14:05 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include <r3.h>
#include <r2.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>

#include <tf_camera.h>
#include <tf_calib.h>
#include <tf_math.h>

#include <tf_targets.h>

targets_data_t tf_targets_data_create (int32_t size)
{
    targets_data_t tdat = (targets_data_t)malloc(sizeof(struct _targets_data_t));
    tdat->ntargets = size;
    tdat->p_w = (r3_t *)malloc(size*sizeof(r3_t));
    tdat->u_w = NULL;
    tdat->v_w = NULL;
    tdat->p_i = (r2_t *)malloc(size*sizeof(r2_t));
    tdat->u_i = NULL;
    tdat->v_i = NULL;
    tdat->alpha = NULL;
    tdat->beta = NULL;
    tdat->conf = rn_alloc(size);
    return tdat;
}

void tf_targets_data_free (targets_data_t tdat)
{
  if (tdat->type != NULL) free(tdat->type);
  if (tdat->p_w != NULL) free(tdat->p_w);
  if (tdat->u_w != NULL) free(tdat->u_w);
  if (tdat->v_w != NULL) free(tdat->v_w);
  if (tdat->p_i != NULL) free(tdat->p_i);
  if (tdat->u_i != NULL) free(tdat->u_i);
  if (tdat->v_i != NULL) free(tdat->v_i);
  if (tdat->alpha != NULL) free(tdat->alpha);
  if (tdat->beta != NULL) free(tdat->beta);
  if (tdat->conf != NULL) free(tdat->conf);
  free(tdat);
}


void tf_targets_data_clear (targets_data_t tdat)
{ 
    int32_t k;
    for (k = 0; k < tdat->ntargets; k++) {
        tdat->p_i[k].c[0] = 0.0; 
        tdat->p_i[k].c[1] = 0.0; 
        tdat->p_w[k].c[0] = 0.0;
        tdat->p_w[k].c[1] = 0.0;
        tdat->p_w[k].c[2] = 0.0;
        tdat->conf[k] = 0.0;
    }
}

void tf_targets_data_copy (targets_data_t source, targets_data_t destination)
{ 
    int32_t k;
    destination->ntargets = source->ntargets;
    for (k = 0; k < source->ntargets; k++) {
        destination->p_i[k] = source->p_i[k]; 
        destination->p_w[k] = source->p_w[k];
        destination->conf[k] = source->conf[k];
    }
}

void tf_targets_write_info
  ( int32_t target, 
    r2_t pos, 
    double alpha, 
    double beta, 
    double errorSqr, 
    double confidence, 
    FILE *f
  )
{
  fprintf(f, "target %3d  pos (%7.2f, %7.2f)",  target, pos.c[0], pos.c[1]);
  fprintf(f, "  alpha %6.4f  beta %6.4f",  alpha, beta);
  fprintf(f, "  error %6.4f  conf %10.8f\n",  sqrt(errorSqr), confidence);
}

targets_data_t tf_targets_data_read 
  ( char *p_w_fname,
    char *uv_w_fname,
    char *p_i_fname,
    char *conf_fname )
{
  int32_t n_p_w, n_p_i, n_uv_w, n_conf;
  targets_data_t tdat = (targets_data_t)notnull(malloc(sizeof(struct _targets_data_t)), "no mem");
  tf_calib_data_read_world_points (p_w_fname, &n_p_w, &(tdat->p_w));
  tf_read_types_and_world_shapes (uv_w_fname, &n_uv_w, &(tdat->type), &(tdat->u_w), &(tdat->v_w));
  demand(n_p_w == n_uv_w, "inconsistent target counts (world_coords vs. types/shapes)");
  tf_calib_data_read_image_points (p_i_fname, &n_p_i, &(tdat->p_i));
  demand(n_p_i == n_p_w, "inconsistent mark counts (world_coords vs. image_coords)");
  if (conf_fname != NULL) {
    /* read confidence weights from file */
    tf_calib_data_read_weights (conf_fname, &n_conf, &(tdat->conf));
    demand(n_conf == n_p_w, "inconsistent mark counts (world_coords vs. weights)");
  } 
  else {
    /* assume unit weight for all marks */
    tdat->conf = rn_alloc(n_p_w );
    int32_t k;
    for (k = 0; k < n_p_w; k++) { tdat->conf[k] = 1.0; }
  }
  
  /* Allocate other fields and set them to default values: */
  tdat->u_i = (r2_t *)notnull(malloc(n_p_w*sizeof(r2_t)), "no mem");
  tdat->v_i = (r2_t *)notnull(malloc(n_p_w*sizeof(r2_t)), "no mem");
  tdat->alpha = rn_alloc(n_p_w);
  tdat->beta = rn_alloc(n_p_w);
  int32_t i;
  for (i = 0; i < n_p_w; i++) {
    tdat->u_i[i] = (r2_t){{ 0, 0 }};
    tdat->v_i[i] = (r2_t){{ 0, 0 }};
    tdat->alpha[i] = 1.0;
    tdat->beta[i] = 0.0;
  }
  tdat->ntargets = n_p_w;
  return tdat;
}

void tf_read_types_and_world_shapes
  ( char *fname, 
    int32_t *ntargets, 
    int32_t **type,  
    r3_t **u_world_coords, 
    r3_t **v_world_coords
  )
{
  int32_t i, n;
  FILE *f = open_read(fname, TRUE);
  if (f == NULL) {
    fprintf(stderr, "error: the file %s doesn't exist\n", fname);
    exit(1);
  }
  int32_t ok = fscanf(f, "%d ", &n);
  if ((ok != 1) || (n < 1)) {
    fprintf(stderr, "error: number of targets missing or invalid\n");
    exit (-1);
  } 
  *type = (int32_t *)notnull(malloc(n * sizeof(int32_t)), "no mem");
  *u_world_coords = (r3_t *)notnull(malloc(n * sizeof(r3_t)), "no mem");
  *v_world_coords = (r3_t *)notnull(malloc(n * sizeof(r3_t)), "no mem");

    for (i = 0; i < n; i++)
    {
         /*reading type of the target*/
         fscanf(f, "%d ", &((*type)[i])); 

         /*reading u vector coordinates */
         fscanf(f, "%lf ", &((*u_world_coords)[i].c[0]));
         fscanf(f, "%lf ", &((*u_world_coords)[i].c[1]));
         fscanf(f, "%lf ", &((*u_world_coords)[i].c[2]));
        
         /*reading v vector coordinates */
         fscanf(f, "%lf ", &((*v_world_coords)[i].c[0]));
         fscanf(f, "%lf ", &((*v_world_coords)[i].c[1]));
         fscanf(f, "%lf ", &((*v_world_coords)[i].c[2]));
    }
  fclose(f);
  (*ntargets) = n;
}
