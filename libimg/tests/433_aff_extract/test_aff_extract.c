#define PROG_NAME "test_aff_extract"
#define PROG_DESC "test of {float_image_aff_extract.h}"
#define PROG_VERS "1.0"

/* Last edited on 2020-11-06 00:18:29 by jstolfi */ 
/* Created on 2020-11-05 by J. Stolfi, UNICAMP */

#define taffe_COPYRIGHT \
  "Copyright © 2020  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <uint16_image.h>
#include <float_image.h>
#include <float_image_read_gen.h>
#include <float_image_write_gen.h>
#include <float_image_aff_sampling.h>
#include <float_image_aff_extract.h>
#include <indexing.h>
#include <r2.h>
#include <i2.h>
#include <r2x2.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

typedef struct taffe_options_t 
  { char *prefix;     /* Output name prefix. */
    char *imageName;  /* First image name (sans directory and extension). */
    r2_aff_map_t A;   /* Optimum map for first image. */
    double Amag;      /* Extra nag factor for {A}. */
  } taffe_options_t;

float_image_t *taffe_read_image(char *imageName);
  /* Reads "in/{imageName}.png" as a float image. */

void taffe_write_image(float_image_t *img, char *prefix, char *imageName);
  /* Writes float image {img} to "out/{prefix}_{imageName}.png". */

taffe_options_t *taffe_parse_options(int argc, char **argv);
  /* Parses the command line options. */
 
void taffe_parse_next_affine_map(argparser_t *pp, r2_aff_map_t *A);
  /* Parses an affine map from the command line, as 6 numbers 
    {m[0][0] m[0][1] m[1][0] m[1][1] d[0] d[1]},
    where {m = A->mat} and {d = A->disp}. */
    
int main(int argn, char **argv);
  
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Parse arguments: */
    taffe_options_t *o = taffe_parse_options(argc, argv);

    /* Input image and affine map: */
    float_image_t *img = taffe_read_image(o->imageName);
    r2_aff_map_t A = o->A;
    r2x2_scale(o->Amag, &(A.mat), &(A.mat));
    
    /* Choose sampling step and sampling grid size: */
    r2_t dp = float_image_aff_sampling_choose_step(&(A.mat));
    fprintf(stderr, "step = (%.6f %.6f)\n", dp.c[0], dp.c[1]);
    double R = 3.5; /* Should be enough. */
    i2_t size = float_image_aff_sampling_grid_size(dp, R);
    fprintf(stderr, "size = (%d %d)\n", size.c[0], size.c[1]);
   
    /* Extract the feature: */
    float_image_t *res = float_image_aff_extract(img, &A, dp, size);
    
    /* Write it out: */
    taffe_write_image(res, o->prefix, o->imageName);
    
    return 0;
  }
  
float_image_t *taffe_read_image(char *imageName)
  {
    char *fname = NULL; 
    asprintf(&fname, "in/%s.png", imageName);
    image_file_format_t ffmt = image_file_format_PNG;
    uint16_t *maxval; /* Max file sample value per channel. */
    double gammaDec, bias; /* Sample decoding parameters. */
    float_image_t *img = float_image_read_gen_named(fname, ffmt, 0.0, 1.0, &maxval, &gammaDec, &bias, FALSE);
    fprintf(stderr, "gammaDec = %.6f 1/gammaDec = %.6f  bias = %.6f\n", gammaDec, 1/gammaDec, bias);
    free(fname);
    return img;
  }
  
void taffe_write_image(float_image_t *img, char *prefix, char *imageName)
  {
    char *fname = NULL; 
    asprintf(&fname, "out/%s_%s.png", prefix, imageName);
    image_file_format_t ffmt = image_file_format_PNG;
    float_image_write_gen_named(fname, img, ffmt, 0.0, 1.0, NAN, NAN, FALSE);
    free(fname);
  }
  
taffe_options_t *taffe_parse_options(int argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    taffe_options_t *o = notnull(malloc(sizeof(taffe_options_t)), "no mem");

    o->prefix = argparser_get_next(pp);
    
    o->imageName = argparser_get_next(pp);
    taffe_parse_next_affine_map(pp, &(o->A)); 
    o->Amag = argparser_get_next_double(pp, 0.1, 100.0);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void taffe_parse_next_affine_map(argparser_t *pp, r2_aff_map_t *A)
  {
    for (int32_t i = 0; i < 2; i++)
      { for (int32_t j = 0; j < 2; j++)
          { A->mat.c[i][j] = argparser_get_next_double(pp, -100.0, +100.0); }
      }
    for (int32_t j = 0; j < 2; j++)
      { A->disp.c[j] = argparser_get_next_double(pp, -100.0, +100.0); }
  }

    
