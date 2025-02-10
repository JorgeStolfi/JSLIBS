#define PROG_NAME "test_aff_compare"
#define PROG_DESC "test of {float_image_aff_compare.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-30 04:58:13 by stolfi */ 
/* Created on 2020-09-26 by J. Stolfi, UNICAMP */

#define taffc_COPYRIGHT \
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
#include <float_image_aff_compare.h>
#include <ix.h>
#include <r2.h>
#include <i2.h>
#include <r2x2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_pmap_affine.h>
#include <hr2_pmap_translation.h>
#include <argparser_geo.h>
#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <affirm.h>

typedef double squared_mismatch_func_t(int i, int scale, double x, double y); 
  /* Evaluates image number {i} reduced to scale {scale} at point {x,y}. */

typedef struct taffc_options_t 
  { char *prefix;     /* Output name prefix. */
    char *deform;     /* Kind of deformation to apply. */
    int32_t NS;       /* Number of steps to take on each side of zero. */
    double uvMax;     /* Max abs value of {u,v} parameters. */
    char *img1;       /* First image name (sans directory and extension). */
    hr2_pmap_t A1;    /* Optimum map for first image. */
    char *img2;       /* Second image name (sans directory and extension). */
    hr2_pmap_t A2;    /* Optimum map for second image. */
  } taffc_options_t;
  
typedef hr2_pmap_t taffc_make_map_func_t (double u, double v);
  /* Type of a procedure that returns an affine map {M(u,v)} a as a function 
    of two parameters {u,v} varying over {[-1 _ +1]}. 
    It should return the identity for {0,0}. It should assume
    a comparison region of radius 1. */

void taffc_plot_square_mismatch
  ( char *prefix,
    char *deform,
    float_image_t *img1,
    hr2_pmap_t *A1,
    float_image_t *img2,
    hr2_pmap_t *A2,
    int32_t NS,
    double uvMax,
    taffc_make_map_func_t *M
  );
  /* Writes a  file "out/{prefix}_{deform}_plot.dat" with a 2D slice of the function
    {float_image_aff_compare} with parameters {img1,C1(u,v),img2,C2(u,v)}, 
    for variable affine maps {C1(u,v),C2(u,v)}.
    
    The affine map {C1(u,v)} is defined as {C1(u,v)(p) = A1(M(u,v)(p))}. 
    
    The map {C2(u,v)} is defined as {C2(u,v)(p) = A2(M(u,v)^{-1}(p))}. .
    
    A total of {2*NS-1} values of {u} and {v} are used, varying over {[-uvMax _ +uvMax]}. */

#define taffc_translation_REF (0.5)
hr2_pmap_t taffc_make_xy_translation_map(double u, double v);
  /* Generates an affine map that displaces the origin by {(D*u,D*v)},
    where {D = taffc_translation_REF}. */
  
#define taffc_rotation_REF (0.60*M_PI)
#define taffc_scale_REF (2.0)
hr2_pmap_t taffc_make_rotate_scale_map(double u, double v);
  /* Generates an affine map that rotates the plane by {T*u}
    radians and magnifies it by {S^v}, where {T = taffc_rotation_REF}
    and {S = taffc_scale_REF}. */
 
#define taffc_stretch_REF (1.2)
hr2_pmap_t taffc_make_xy_stretch_map(double u, double v);
  /* Generates an affine map that scales {x} and {y} by {S^u}
    and {S^v}, respectively, where {S = taffc_stretch_REF}. */

#define taffc_shear_REF (0.5)
hr2_pmap_t taffc_make_xy_shear_map(double u, double v);
  /* Generates an affine map that shears the plane by {S*u} in the {y}
    direction and by {S*v} in the {x} diiection, where 
    {S = taffc_shear_REF}. Namely, maps {(1,0)} to {(1,S*u)} and
    {(0,1)} to {(S*u,1)}. */

float_image_t *taffc_read_image(char *imgname);
  /* Reads "in/{imgname}.png" as a float image. */

taffc_options_t *taffc_parse_options(int argc, char **argv);
  /* Parses the command line options. */
    
void taffc_scale_pmap_aff(hr2_pmap_t *A, double scale);
  /* Modifies the map {A} by composing it with a uniform scalinh by {scale},
    in that order. */

void taffc_show_map(char *name, char *qualif, hr2_pmap_t *A);
  /* Prints the map {A} on {stderr}, labeled with the given {name}
    and {qualif}. */

void taffc_show_map_pair(char *phase, hr2_pmap_t *C1, hr2_pmap_t *C2);
  /* Prints to {stderr} the maps {C1} and {C2}, as well as the composition
    of the inverse of {C1} with {C2}, applied in that order.  The latter is the
    map that morphs image 1 to match image 2. */

int main(int argn, char **argv);

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    srandom(4615*417);
    
    /* Parse arguments: */
    taffc_options_t *o = taffc_parse_options(argc, argv);

    /* Input images and affine maps: */
    float_image_t *img1 = taffc_read_image(o->img1);
    hr2_pmap_t A1 = o->A1;
    assert(hr2_pmap_is_affine(&A1, 1.0e-12));
    
    float_image_t *img2 = taffc_read_image(o->img2);
    hr2_pmap_t A2 = o->A2;
    assert(hr2_pmap_is_affine(&A2, 1.0e-12));
    
    int32_t NS = o->NS; /* Samples of {u,v} on each side of zero. */
    
    if (strcmp(o->deform, "xytrans") == 0)
      { taffc_plot_square_mismatch
          ( o->prefix, "xytrans", img1, &A1, img2, &A2, NS, o->uvMax, &taffc_make_xy_translation_map );
      }
    else if (strcmp(o->deform, "rotmag") == 0)
      { taffc_plot_square_mismatch
          ( o->prefix, "rotmag", img1, &A1, img2, &A2, NS, o->uvMax, &taffc_make_rotate_scale_map );
      }
    else if (strcmp(o->deform, "xystretch") == 0)
      { taffc_plot_square_mismatch
          ( o->prefix, "xystretch", img1, &A1, img2, &A2, NS, o->uvMax, &taffc_make_xy_stretch_map );
      }
    else if (strcmp(o->deform, "xyshear") == 0)
      { taffc_plot_square_mismatch
          ( o->prefix, "xyshear", img1, &A1, img2, &A2, NS, o->uvMax, &taffc_make_xy_shear_map );
      }
    else 
      { demand(FALSE, "invalid {deform}"); }
    
    return 0;
  }
  
hr2_pmap_t taffc_make_xy_translation_map(double u, double v)
  {
    double D = taffc_translation_REF;
    r2_t disp = (r2_t){{ D*u, D*v }};
    hr2_pmap_t M = hr2_pmap_translation_from_disp(&disp);
    return M;
  }
  
hr2_pmap_t taffc_make_rotate_scale_map(double u, double v)
  {
    double R = taffc_rotation_REF;
    double S = taffc_scale_REF;
    double ang = R*u;
    double scale = pow(S, v);
    hr2_pmap_t M = hr2_pmap_rotation_and_scaling(ang, scale);
    return M;
  }
  
hr2_pmap_t taffc_make_xy_stretch_map(double u, double v)
  {
    double S = taffc_stretch_REF;
    r2_t scale = (r2_t){{ pow(S,u), pow(S,v) }};
    hr2_pmap_t M = hr2_pmap_scaling(&scale);
    return M;
  }
  
hr2_pmap_t taffc_make_xy_shear_map(double u, double v)
  {
    double S = taffc_shear_REF;
    r2x2_t L;
    r2x2_ident(&L);
    L.c[0][1] = S*u;
    L.c[1][0] = S*v;
    r2_t disp = (r2_t){{ 0.0, 0.0 }};
    hr2_pmap_t M = hr2_pmap_affine_from_mat_and_disp(&L, &disp);
    return M;
  }
  
void taffc_scale_pmap_aff(hr2_pmap_t *A, double scale)
  {
    for (int32_t i = 0;  i < 3; i++)
      { for (int32_t j = 1;  j < 3; j++) 
         { A->dir.c[i][j] *= scale; A->inv.c[j][i] /= scale; }
      }
  }

float_image_t *taffc_read_image(char *imgname)
  {
    char *fname = jsprintf("in/%s.png", imgname);
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE;
    uint16_t *maxval; /* Max file sample value per channel. */
    double gammaDec, bias; /* Sample decoding parameters. */
    float_image_t *img = float_image_read_gen_named
      ( fname, ffmt, yUp, 0.0, 1.0, &maxval, &gammaDec, &bias, FALSE);
    free(fname);
    return img;
  }
  
void taffc_plot_square_mismatch
  ( char *prefix,
    char *deform,
    float_image_t *img1,
    hr2_pmap_t *A1,
    float_image_t *img2,
    hr2_pmap_t *A2,
    int32_t NS,
    double uvMax,
    taffc_make_map_func_t *M
  )
  {
    bool_t debug_maps = FALSE;
    char *fname = jsprintf("out/%s_%s_plot.dat", prefix, deform);
    FILE *wr = open_write(fname, TRUE);
    
    taffc_show_map_pair("initial", A1, A2);
    
    double fopt = +INF;  /* Optimum square msismatch value {(u,v)}. */
    r2_t popt = (r2_t){{ NAN, NAN }};    /* Optimum parameter pair {(u,v)}. */
    i2_t iopt = (i2_t){{ 9999, 9999 }};    /* Optimum parameter sample index {(iu,iv)}. */
    hr2_pmap_t C1opt, C2opt;    /* Optimum affine maps. */
    r2_t dp; /* Sampling steps used by {float_image_aff_compare}. */
    i2_t size; /* Sampling grid size used by {float_image_aff_compare}. */

    for (int32_t iu = -NS; iu <= +NS; iu++)
      { double u = uvMax*((double)iu)/((double)NS);
        for (int32_t iv = -NS; iv <= +NS; iv++)
          { double v = uvMax*((double)iv)/((double)NS);
            
            /* Construct the modified maps {C1, C2}: */
            hr2_pmap_t Muv = M(u,v);
            hr2_pmap_t C1 = hr2_pmap_compose(&Muv, A1);
            if (debug_maps) { taffc_show_map("C1", "", &C1); }
            
            hr2_pmap_t Nuv = hr2_pmap_inv(&Muv);
            hr2_pmap_t C2 = hr2_pmap_compose(&Nuv, A2);
            if (debug_maps) { taffc_show_map("C2", "", &C2); }
            
            double f = float_image_aff_compare(img1, &C1, img2, &C2, &dp, &size);
            if (debug_maps) { fprintf(stderr, "point %4d %4d u = %10.7f v = %10.7f mismatch = %12.9f\n", iu, iv, u, v, f); }
            if (f < fopt)
              { /* Update optimum: */
                fopt = f;  popt = (r2_t){{ u, v }}; iopt = (i2_t){{ iu, iv }};
                C1opt = C1; C2opt = C2;
              }
            if (debug_maps) { fprintf(stderr, "\n"); }
              
            fprintf(wr, "0 %4d %4d  %10.7f %10.7f  %12.9f\n", iu, iv, u, v, f);
          }
        fprintf(wr, "\n"); /* For gnuplot. */
      }
    //  fprintf
    //    ( wr, "1 %4d, %4d %10.7f %10.7f %16.8e\n", 
    //      iopt.c[0], iopt.c[0], popt.c[0], popt.c[1], fopt
    //    );
    fclose(wr);
    fprintf
      ( stderr, "lowest mismatch %.8e at [%d, %d] = (%.7f %.7f)\n", 
        fopt, iopt.c[0], iopt.c[0], popt.c[0], popt.c[1]
      );
      
    taffc_show_map_pair("optimum", &C1opt, &C2opt);
  }
    
taffc_options_t *taffc_parse_options(int argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, "HELP!");
    argparser_set_info(pp, "KNOW THYSELF");
    argparser_process_help_info_options(pp);
    
    taffc_options_t *o = notnull(malloc(sizeof(taffc_options_t)), "no mem");

    o->prefix = argparser_get_next(pp);
    o->deform = argparser_get_next(pp);
    
    o->NS = (int32_t)argparser_get_next_uint(pp, 1, 500);
    
    o->uvMax = argparser_get_next_double(pp, 0.01, 100.0);
    
    o->img1 = argparser_get_next(pp);
    o->A1 = argparser_get_next_feature_map(pp);
    
    o->img2 = argparser_get_next(pp);
    o->A2 = argparser_get_next_feature_map(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }

void taffc_show_map(char *name, char *phase, hr2_pmap_t *A)
  { double w = fabs(A->dir.c[0][0]);
    fprintf(stderr, "%-4s %s = [", name, phase);
    fprintf(stderr, " [ %+10.4f %+10.4f ]", A->dir.c[1][1]/w, A->dir.c[1][2]/w);
    fprintf(stderr, " [ %+10.4f %+10.4f ]", A->dir.c[2][1]/w, A->dir.c[2][2]/w);
    fprintf(stderr, " ]");
    fprintf(stderr, " [ %10.4f %10.4f ]", A->dir.c[0][1]/w, A->dir.c[0][2]/w);
    fprintf(stderr, "\n");
  }
            
void taffc_show_map_pair(char *phase, hr2_pmap_t *C1, hr2_pmap_t *C2)
  { 
    /* Print the maps from {\RR^2} to {img1} and {img2}: */
    taffc_show_map("C1", phase, C1);
    taffc_show_map("C2", phase, C2);
   
    /* Show the implied map from {img1} to {img2}: */
    hr2_pmap_t C1inv = hr2_pmap_inv(C1);
    hr2_pmap_t C12 = hr2_pmap_compose(&C1inv, C2);
    taffc_show_map("C12", phase, &C12);

    /* Print the normalized map from {img1} to {img2} and the relative scale radius: */
    hr2_pmap_t A12 = C12;
    double R12 = hypot(hypot(A12.dir.c[1][1], A12.dir.c[1][2]), hypot(A12.dir.c[2][1], A12.dir.c[2][2]));
    taffc_scale_pmap_aff(&(A12), 1/R12);
    taffc_show_map("A12", phase, &A12);
    fprintf(stderr, "R12 %s = %10.4f\n", phase, R12);
  }
    
