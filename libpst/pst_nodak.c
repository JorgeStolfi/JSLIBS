/* See pst_nodak.h */
/* Last edited on 2024-11-30 22:55:24 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_average.h>
#include <r2.h> 
#include <r3.h> 
#include <r3x3.h> 
#include <rn.h> 
#include <hr2.h>
#include <hr2_pmap.h>
#include <qmin_simplex.h> 
#include <affirm.h>
#include <argparser.h>
#include <ix_reduce.h>

#include <pst_basic.h>
#include <pst_nodak.h>

/* INTERNAL PROTOTYPES */

r3x3_t pst_nodak_get_matrix_from_4_points(r2_t p[]);
  /* Computes the matrix that maps the unit square corners {(0,0), (1,0), (1,1), (0,1)} 
    to the points {p[0..3]}. */

r3x3_t pst_nodak_get_matrix_4(r2_t geo_p[], r2_t img_p[]);
  /* Computes the matrix that maps the four chart points {geo_p[0..3]} 
    to the four image points {img_p[0..3]}. */

void debug_matrix(char *head, r3x3_t *M, char *foot, double R);
  /* Prints on {stderr} the matrix {M} and its effect on the corners 
    of a rectangle {[-R_+R]^2}. */

/* IMPLEMENTATIONS */

void debug_matrix(char *head, r3x3_t *M, char *foot, double R)
  { if (head != NULL) { fprintf(stderr, "  %s\n", head); }
    int x, y;
    for (y = -1; y <= 1; y++)
      for (x = -1; x <= 1; x++)
        { r3_t p = (r3_t){{ 1.0, x*R, y*R }};
          r3_t q;
          r3x3_map_row(&p, M, &q);
          r2_t t;
          t.c[0] = q.c[1]/q.c[0]; 
          t.c[1] = q.c[2]/q.c[0];
          rn_gen_print(stderr, 3, p.c, "%+7.2f", "  [ ", " ", " ]");
          rn_gen_print(stderr, 3, q.c, "%+7.2f", " = [ ", " ", " ]");
          rn_gen_print(stderr, 2, t.c, "%+7.2f", " = ( ", " ", " )");
          fprintf(stderr, "\n");
        }
     if (foot != NULL) { fprintf(stderr, "  %s\n", foot); }
  }

r3x3_t pst_nodak_get_matrix(r2_vec_t* geo_ctr, double geo_radius, int32_vec_t* img_num, r2_vec_t* img_ctr)
  {
    bool_t debug = TRUE;

    /* Compute a number of matrices from 4 point pairs. */
    int n = img_num->ne; /* Number of data points. */
    demand(img_ctr->ne == n,"Wrong list length");
    int m = (n+3)/4;    /* Number of matrices. */
    int step = n/m;     /* Displacement between data subsets. */
    r3x3_t sum_Q;       /* Sum of matrices computed from each data subset. */
    r3x3_zero(&sum_Q);
    int i;
    for(i = 0; i < m; i++)
      { /* Extract the data subset number {i}. */
        r2_t geo_p[4]; 
        r2_t img_p[4];
        int j;
        for(j = 0; j < 4; j++)
          { int ij = (i*step + j)%n;
            img_p[j] = img_ctr->e[ij];
            int k = img_num->e[ij];
            demand( (k >= 0) && (k < geo_ctr->ne), "Invalid spot index");
            geo_p[j] = geo_ctr->e[k];
          }
        /* Compute the matrix {Q} from the data subset. */
        r3x3_t Q = pst_nodak_get_matrix_4(geo_p,img_p);
        if (debug) { debug_matrix("Global Matrix:", &Q, "",geo_radius); }
        /* Normalize the matrix and add to {sum_Q}. */
        double Qmag = r3x3_norm(&Q);
        r3x3_mix(1.0, &sum_Q, 1.0/Qmag, &Q, &sum_Q);
      }
    return sum_Q;
  }

r3x3_t pst_nodak_get_matrix_from_4_points(r2_t p[]) 
  {
    hr2_point_t p0 = (hr2_point_t){{{1.0,p[0].c[0],p[0].c[1]}}};
    hr2_point_t p1 = (hr2_point_t){{{1.0,p[1].c[0],p[1].c[1]}}};
    hr2_point_t p2 = (hr2_point_t){{{1.0,p[2].c[0],p[2].c[1]}}};
    hr2_point_t p3 = (hr2_point_t){{{1.0,p[3].c[0],p[3].c[1]}}};
    hr2_pmap_t P =  hr2_pmap_from_four_points(&p0, &p1, &p2, &p3);
    return P.dir;
  }

r3x3_t pst_nodak_get_matrix_4(r2_t geo_p[], r2_t img_p[])
  {
    bool_t debug = TRUE;
   
    r3x3_t img_P = pst_nodak_get_matrix_from_4_points(img_p);
    r3x3_t geo_P = pst_nodak_get_matrix_from_4_points(geo_p);

    r3x3_t geo_IP;
    r3x3_inv(&geo_P,&geo_IP);
    
    if (debug) { debug_matrix("Image Matrix:", &img_P, "",1.0); }
    if (debug) { debug_matrix("Chart Matrix:", &geo_P, "",1.0); }

    r3x3_t M;
    r3x3_mul(&geo_IP,&img_P,&M);
    return M;
  }

float_image_t *pst_nodak_extract_chart
  ( float_image_t *img,   /* Photo of a scene that includes a N-Spot chart. */
    double rad,           /* Chart radius in chart coordinates. */
    r3x3_t *C2I,          /* Chart-to-image projective map matrix. */
    int OSZ               /* Width and height of output image (pixels) */
  )
  { /* Get input image channels: */
    int NC = (int)(img->sz[0]);
    
    /* Output image dimensions (pixels): */
    int ONX = OSZ;
    int ONY = OSZ;
    
    /* Edges of chart, in chart coordinates: */
    double xlo = -rad;
    double xhi = +rad;
    double ylo = -rad;
    double yhi = +rad;
    
    /* Create and fill the chart image: */
    float_image_t *omg = float_image_new(NC, ONX, ONY);

    ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
    float undef = 0.5;
    bool_t avg = TRUE;
    int order = 1;
    float_image_transform_copy_persp_rectangle
      (img, red, xlo, xhi, ylo, yhi, C2I, undef, avg, order, 0, 0, ONX, ONY, NULL, omg);

    /* Should clear the background around the chart. */
    return omg;
  }

float_image_t *pst_nodak_extract_gray_scale
  ( float_image_t *img,    /* Photo of a scene that includes a N-Spot chart. */
    r2_vec_t* geo_ctr,     /* Center of each spot in chart coordinates. */
    double_vec_t* geo_rad, /* Radius of each spot in chart coordinates. */
    r3x3_t *C2I,           /* Chart-to-image projective map matrix. */
    double mrg,            /* Safety margin width in pixels. */
    int NX,                /* Width of each patch in the output image (pixels) */
    int NY                 /* Height of each patch in the output image (pixels) */
  )
  { 
    /* Get number of channels:*/
    int NC = (int)(img->sz[0]);
    /* Create output image: */
    int NS = geo_ctr->ne;
    demand(geo_rad->ne == NS, "inconsistent spot {ctr}/{rad}");
    int ONX = NX*NS;
    int ONY = NY;
    float_image_t *omg = float_image_new(NC, ONX, ONY);
    /* Fill in the patches: */
    int i;
    for (i = 0; i < NS; i++)
      { /* Get spot center and radius in chart coordinates:*/
        r2_t *ctr = &(geo_ctr->e[i]);
        double rad = geo_rad->e[i];
        /* Compute spot average color: */
        float avg[NC];
        ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
        float_image_average_persp_disk(img, red, ctr, rad, C2I, mrg, avg);
        /* Fill patch in output image: */
        int x, y;
        for (y = 0; y < NY; y++)
          { for (x = 0; x < NX; x++)
              { float_image_set_pixel(omg, i*NX + x, y, avg); }
          }
      }
    return omg;
  }


