/* See {float_image_transform.h}. */
/* Last edited on 2021-06-09 20:56:05 by jstolfi */

#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>
#include <affirm.h>
#include <indexing.h>
#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_average.h>

/* INTERNAL PROTOTYPES */

#define SQRT3 (1.73205080756887729352)

/* IMPLEMENTATIONS */

void float_image_transform_all
  ( float_image_t *iimg,    /* Input image. */
    ix_reduction_t red,     /* Index reduction method. */ 
    r2_map_jacobian_t *map, /* Output-to-input coordinate transformation. */
    float undef,            /* Sample value for undefined output pixels. */
    bool_t avg,             /* TRUE to average pixels, FALSE to add them. */
    int order,              /* Interpolation order. */
    r2_pred_t *debugp,      /* Tells whether pixel should be debugged. */
    float_image_t *oimg     /* Output image. */
  )
  { 
    int ocols = (int)oimg->sz[1];
    int orows = (int)oimg->sz[2];
    float_image_transform_sub(iimg, red, map, undef, avg, order, 0, 0, ocols, orows, debugp, oimg);
  }    

void float_image_transform_sub
  ( float_image_t *iimg,    /* Input image. */
    ix_reduction_t red,     /* Index reduction method. */ 
    r2_map_jacobian_t *map, /* Output-to-input coordinate transformation. */
    float undef,            /* Sample value for undefined output pixels. */
    bool_t avg,             /* TRUE to average pixels, FALSE to add them. */
    int order,              /* Interpolation order. */
    int x0,                 /* First output image column. */
    int y0,                 /* First output image row. */
    int NX,                 /* Number of output image columns. */
    int NY,                 /* Number of output image rows. */
    r2_pred_t *debugp,      /* Tells whether pixel should be debugged. */
    float_image_t *oimg     /* Output image. */
  )
  { demand(iimg->sz[0] == oimg->sz[0], "images must have the same channels");
    int chns = (int)oimg->sz[0];
    float fo[chns];
    /* Scan rows from top to bottom to make debugging easier: */
    int row, col;
    for (row = y0+NY-1; row >= y0; row--)
      { for (col = x0; col < x0+NX; col++)
          { r2_t p = (r2_t){{ col + 0.5, row + 0.5 }};
            bool_t debug = (debugp == NULL ? FALSE : debugp(&p));
            if (debug) { fprintf(stderr, "------------------------------------------------------------\n"); }
            /* bool_t debug = FALSE; */
            /* Get transformed pixel: */
            float_image_transform_get_pixel(iimg, red, col, row, map, undef, avg, order, fo, debug);
            if (debug) { float_image_debug_pixel("  fo", p.c[0], p.c[1], chns, fo, "\n"); }
            /* Store it in the image: */
            float_image_set_pixel(oimg, col, row, fo);
            if (debug) { fprintf(stderr, "------------------------------------------------------------\n"); }
          }
      }
  }    

#define MAX_TAYLOR_ERROR (+INF)  
  /* For now */

void float_image_transform_get_pixel
  ( float_image_t *img, 
    ix_reduction_t red, /* Index reduction method. */ 
    int col, 
    int row, 
    r2_map_jacobian_t *map, 
    float undef, 
    bool_t avg,
    int order, 
    float f[],
    bool_t debug        /* If TRUE, prints debugging info. */
  )
  { 
    int chns = (int)img->sz[0];

    /* To make the sampling integrals manageable, we replace the transform
      {map} is by its 1st degree Taylor expansion.  
      
      More precisely, for each output pixel center {op} we approximate
      {map(op+r)}, for small {r}, by {aff(r) = ip + r*J}, where {ip}
      is {map(op) and {J} is the Jacobian matrix of {map} evaluated at
      {op}. This approximation assumes that the higher-order terms of
      the Taylor expansion is negligible as long as {r} lies within
      the sampling kernel's support (typically, a few output pixels).
      The integral then reduces to a sum of input pixel values times
      affinely-transformed images of the reconstruction kernel {h}. */

    /* Get Cartesian coordinates {xc,yc} of destination pixel's center: */
    r2_t op = (r2_t) {{ (double)col + 0.5, (double)row + 0.5 }};
    
    /* Compute corresp point {pt=(xt,yt)} in source image, and Jacobian: */
    r2x2_t J; /* Jacobian {{dxt/dxc, dxt/dyc}, {dyt/dxc, dyt/dyc}} */
    r2_t ip = op;
    r2x2_ident(&J);
    map(&ip, &J);
    bool_t invalid = (isnan(ip.c[0]) || isnan(ip.c[1]));
    
    if (debug) { fprintf(stderr, "  computing output pixel with indices (%d,%d):\n", col,row); }
    if (debug) { r2_debug_point_jac("    po", &op, NULL, "\n"); }
    if (debug) { r2_debug_point_jac("    pi", &ip, &J, "\n"); }
    if (debug & (!invalid)) { r2_map_check_jacobian(&op, map, "map", 1.0e-5, FALSE); }

   /* Get input image value at point {pt}: */
    int ic;
    if (invalid)
      { /* There is no image on the negative side of the two-sided plane: */
        for (ic = 0; ic < chns; ic++) { f[ic] = undef; }
      }
    else
      { /* demand(err <= MAX_TAYLOR_ERROR, "excessive warp, should use subsampling"); */
        /* Sample points of input image: */
        float_image_average_parallelogram(img, red, &ip, &J, avg, order, f, debug);
        for (ic = 0; ic < chns; ic++) { if (isnan(f[ic])) { f[ic] = undef; } }
      }
  }

void float_image_transform_copy_persp_rectangle
  ( float_image_t *iimg,
    ix_reduction_t red, /* Index reduction method. */ 
    double xlo,         /* Min X in true coords. */
    double xhi,         /* Max X in true coords. */
    double ylo,         /* Min Y in true coords. */
    double yhi,         /* Max Y in true coords. */
    r3x3_t *T2I,        /* Projective map from true coords to image coords. */
    float undef,        /* Defaut for undefined pixels. */
    bool_t avg,         /* TRUE to compute average. */
    int order,          /* Interpolation order to use. */
    int x0,             /* First output image column. */
    int y0,             /* First output image row. */
    int NX,             /* Number of output image columns. */
    int NY,             /* Number of output image rows. */
    r2_pred_t *debugp,  /* Tells whether pixel should be debugged. */
    float_image_t *oimg /* Output image. */
  )
  { 
    if ((NX == 0) || (NY == 0)) { /* Nothing to do: */ return; }

    double xscale = (xhi - xlo)/NX;
    double yscale = (yhi - ylo)/NY;
    
    /* Assemble scaling matrix {S} from output image coords to true coords: */
    r3x3_t S = (r3x3_t){{{ 1.0, xlo - xscale*x0, ylo - yscale*y0}, { 0.0, xscale, 0.0 }, { 0.0, 0.0, yscale }}};
    
    /* Compute total projective map {M}: */
    r3x3_t M;
    r3x3_mul(&S, T2I, &M);
    
    auto void projmap(r2_t *p, r2x2_t *J);
      /* Applies the projective map to the point {p}. Also post-multiplies {J} by the local Jacobian, as
        expected from a {r2_map_jacobian_t}. */
    
    float_image_transform_sub(iimg, red, &projmap, undef, avg, order, x0, NX, y0, NY, debugp, oimg);
    
    return;
    
    void projmap(r2_t *p, r2x2_t *J)
      { r2_map_projective(p, &M, J); }
  }

