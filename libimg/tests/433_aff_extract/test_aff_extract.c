#define PROG_NAME "test_aff_extract"
#define PROG_DESC "test of {float_image_aff_extract.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-30 04:57:57 by stolfi */ 
/* Created on 2020-11-05 by J. Stolfi, UNICAMP */

#define taffe_COPYRIGHT \
  "Copyright © 2020  by the State University of Campinas (UNICAMP)"

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
#include <ix.h>
#include <r2.h>
#include <hr2.h>
#include <i2.h>
#include <r2x2.h>
#include <argparser_geo.h>
#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <affirm.h>

typedef struct taffe_options_t 
  { char *prefix;     /* Output name prefix. */
    char *imageName;  /* Image name (sans directory and extension). */
    hr2_pmap_t A;     /* Affine map to use. */
  } taffe_options_t;

float_image_t *taffe_read_image(char *imageName);
  /* Reads "in/{imageName}.png" as a float image. */

void taffe_write_image(float_image_t *img, char *prefix, char *imageName);
  /* Writes float image {img} to "out/{prefix}_{imageName}.png". */

void taffe_show_map(char *name, char *qualif, hr2_pmap_t *A);
  /* Prints the map {A} on {stderr}, labeled with the given {name}
    and {qualif}. */

taffe_options_t *taffe_parse_options(int argc, char **argv);
  /* Parses the command line options. */
   
int main(int argn, char **argv);
  
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Parse arguments: */
    taffe_options_t *o = taffe_parse_options(argc, argv);

    /* Input image and affine map: */
    float_image_t *img = taffe_read_image(o->imageName);
    hr2_pmap_t A = o->A;
    demand(hr2_pmap_is_affine(&A, 1.0e-12), "map is not affine");
    taffe_show_map("A", "", &A);
    
    /* Choose sampling step and sampling grid size: */
    r2_t dp = float_image_aff_sampling_choose_step(&(A.dir));
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
    char *fname = jsprintf("in/%s.png", imageName);
    image_file_format_t ffmt = image_file_format_PNG;
    uint16_t *maxval; /* Max file sample value per channel. */
    double gammaDec, bias; /* Sample decoding parameters. */
    bool_t yUp = FALSE;
    float_image_t *img = float_image_read_gen_named(fname, ffmt, yUp, 0.0, 1.0, &maxval, &gammaDec, &bias, FALSE);
    fprintf(stderr, "gammaDec = %.6f 1/gammaDec = %.6f  bias = %.6f\n", gammaDec, 1/gammaDec, bias);
    free(fname);
    return img;
  }
  
void taffe_write_image(float_image_t *img, char *prefix, char *imageName)
  {
    char *fname = jsprintf("out/%s_%s.png", prefix, imageName);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE;
    float_image_write_gen_named(fname, img, ffmt, yUp, 0.0, 1.0, NAN, NAN, FALSE);
    free(fname);
  }
 
void taffe_show_map(char *name, char *phase, hr2_pmap_t *A)
  { double w = fabs(A->dir.c[0][0]);
    fprintf(stderr, "%-4s %s = [", name, phase);
    fprintf(stderr, " [ %+10.4f %+10.4f ]", A->dir.c[1][1]/w, A->dir.c[1][2]/w);
    fprintf(stderr, " [ %+10.4f %+10.4f ]", A->dir.c[2][1]/w, A->dir.c[2][2]/w);
    fprintf(stderr, " ]");
    fprintf(stderr, " [ %10.4f %10.4f ]", A->dir.c[0][1]/w, A->dir.c[0][2]/w);
    fprintf(stderr, "\n");
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

    o->A = argparser_get_next_feature_map(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
    
