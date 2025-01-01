#define PROG_NAME "test_hxg"
#define PROG_DESC "tests the {libhxg} routines"
#define PROG_VERSION "1.0"

/* Last edited on 2025-01-01 02:39:55 by stolfi */ 

#define PROG_C_COPYRIGHT "  Run \"" PROG_NAME " -info\" for details"

#define PROG_HELP \
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -mapName {MAP_NAME} ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2024-12-31 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2024-12-31 Created based on {test_cpk.c}. J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a map description and draws it on an regular hexagonal grid."
  
#define PROG_INFO_OPTS \
  "  -mapName {MAP_NAME}\n" \
  "    This mandatory option specifies the name of the map, that will be read from \"in/{MAP_NAME}-urb.txt\"."


#define USAGE "test_hxg [-verbose] [-plot] [-mag4] INPUT_DIR OUTPUT_DIR"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <fget.h>
#include <r2.h>
#include <i2.h>
#include <interval.h>

#include <hxg_canvas.h>
#include <hxg_paint.h>
#include <hxg_eps.h>

/* INTERNAL PROTOTYPES */

void txhg_scale_point(hxg_canvas_t *cvs, r2_t *p);
void txhg_scale_length(hxg_canvas_t *cvs, double *r);
void thxg_paint_a_polygon(hxg_canvas_t *cvs, hxg_paint_op_t *op);
void thxg_paint_a_rect_stroke(hxg_canvas_t *cvs, hxg_paint_op_t *op);
void thxg_paint_a_circle(hxg_canvas_t *cvs, hxg_paint_op_t *op);
void thxg_paint_a_sausage(hxg_canvas_t *cvs, hxg_paint_op_t *op);

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {   
    uint32_t nx =  50;
    uint32_t ny =  70;
    hxg_canvas_t *cvs = hxg_canvas_new(nx, ny, 0.05, TRUE);
    interval_t BB[2];
    hxg_canvas_bbox(cvs, BB);
    
    double *pixv = talloc(nx*ny, double);
    for (int32_t k = 0; k < nx*ny; k++) { pixv[k] = 0.0; }

    /* Seed for random numbers: */
    uint32_t seed = 4635*417;
    srand(seed);
    srandom(seed);
    
    double val;
    auto void set_pixv(uint32_t ix, uint32_t iy, uint32_t k);
      /* A canvas-painting procedure that sets {pixv[k]}
        to {val}.  The indices {ix,iy} are ignored. */ 
    
    val = 1;
    thxg_paint_a_polygon(cvs, set_pixv);
    
    val = 2;
    thxg_paint_a_rect_stroke(cvs, set_pixv);
    
    val = 3;
    thxg_paint_a_circle(cvs, set_pixv);
    
    val = 4;
    thxg_paint_a_sausage(cvs, set_pixv);
    
    uint32_t maxv = 4;
    double *r =  talloc(maxv+1, double);
    frgb_t *fill_color = talloc(maxv+1, frgb_t);
    frgb_t *draw_color = talloc(maxv+1, frgb_t);
    
    for (uint32_t v = 0; v <= maxv; v++)
      { r[v] = (v == 0 ? 0 : 0.05 + 0.10*v); txhg_scale_length(cvs, &(r[v]));
        double f = ((double)v)/maxv;
        float red = (float)(0.500 + sqrt(sqrt(f))*0.500);
        float grn = (float)(0.500 + f*0.500);
        float blu = (float)0.500;
        fill_color[v] = (frgb_t){{ red, grn, blu }};
        for (int32_t c = 0; c < 3; c++)
          { float val = (float)(0.250 + 0.250*fill_color[v].c[c]);
            draw_color[v].c[c] = val;
          }
      }
    
    /* Create the EPS writer: */
    char *fname = "out/test";
    epswr_figure_t *eps = hxg_eps_new_figure(BB, fname);
    hxg_eps_plot_canvas(eps, cvs, pixv, maxv, r, fill_color, draw_color);
    
    epswr_set_pen_color(eps, 1.000, 0.000, 0.000); 
    epswr_rectangle(eps, BB[0].end[0], BB[0].end[1], BB[1].end[0], BB[1].end[1], FALSE, TRUE);
    epswr_end_figure(eps);
    
    return 0;
    
    void set_pixv(uint32_t ix, uint32_t iy, uint32_t k)
      { assert((ix >= 0) && (ix < nx));
        assert((iy >= 0) && (iy < ny));
        pixv[k] = val;
      }
  }

void txhg_scale_point(hxg_canvas_t *cvs, r2_t *p)
  { interval_t BB[2];
    hxg_canvas_bbox(cvs, BB);
    double xlo = BB[0].end[0];
    double xhi = BB[0].end[1];
    double ylo = BB[1].end[0];
    double yhi = BB[1].end[1];
    p->c[0] = xlo + p->c[0]*(xhi - xlo);
    p->c[1] = ylo + p->c[1]*(yhi - ylo);
  }
  
void txhg_scale_length(hxg_canvas_t *cvs, double *r)
  { interval_t BB[2];
    hxg_canvas_bbox(cvs, BB);
    double xlo = BB[0].end[0];
    double xhi = BB[0].end[1];
    double ylo = BB[1].end[0];
    double yhi = BB[1].end[1];
    double scale = sqrt((xhi - xlo)*(yhi-ylo));
    (*r) = (*r) * scale;
  }

void thxg_paint_a_polygon(hxg_canvas_t *cvs, hxg_paint_op_t *op)
  {
    uint32_t n = 11;
    r2_vec_t p = r2_vec_new(n);
    uint32_t k = 0;
    p.e[k] = (r2_t){{ 0.1, 0.1 }}; k++;
    p.e[k] = (r2_t){{ 0.3, 0.1 }}; k++;
    p.e[k] = (r2_t){{ 0.3, 0.5 }}; k++;
    p.e[k] = (r2_t){{ 0.6, 0.1 }}; k++;
    p.e[k] = (r2_t){{ 0.9, 0.1 }}; k++;
    p.e[k] = (r2_t){{ 0.6, 0.5 }}; k++;
    p.e[k] = (r2_t){{ 0.8, 0.6 }}; k++;
    p.e[k] = (r2_t){{ 0.8, 0.8 }}; k++;
    p.e[k] = (r2_t){{ 0.6, 0.9 }}; k++;
    p.e[k] = (r2_t){{ 0.1, 0.9 }}; k++;
    p.e[k] = (r2_t){{ 0.1, 0.1 }}; k++;
    assert(k == n);
    r2_vec_trim(&p, n);
    for (int32_t i = 0; i < n; i++) { txhg_scale_point(cvs, &(p.e[i])); }
    hxg_paint_polygon(cvs, &p, op);
    free(p.e);
  }

void thxg_paint_a_rect_stroke(hxg_canvas_t *cvs, hxg_paint_op_t *op)
  {
    r2_t a = (r2_t){{ 0.2, 0.8 }}; txhg_scale_point(cvs, &a);
    r2_t b = (r2_t){{ 0.3, 0.2 }}; txhg_scale_point(cvs, &b);
    double r = 0.05; txhg_scale_length(cvs, &(r));
    hxg_paint_rect_stroke(cvs, &a, &b, r, op);
  }

void thxg_paint_a_circle(hxg_canvas_t *cvs, hxg_paint_op_t *op)
  {
    r2_t ctr = (r2_t){{ 0.7, 0.7 }}; txhg_scale_point(cvs, &ctr);
    double r = 0.1; txhg_scale_length(cvs, &(r));
    hxg_paint_circle(cvs, &ctr, r, op);
  }

void thxg_paint_a_sausage(hxg_canvas_t *cvs, hxg_paint_op_t *op)
  {
    uint32_t n = 7;
    r2_vec_t p = r2_vec_new(n);
    uint32_t k = 0;
    p.e[k] = (r2_t){{ 0.70, 0.20 }}; k++;
    p.e[k] = (r2_t){{ 0.75, 0.25 }}; k++;
    p.e[k] = (r2_t){{ 0.60, 0.30 }}; k++;
    p.e[k] = (r2_t){{ 0.65, 0.35 }}; k++;
    p.e[k] = (r2_t){{ 0.50, 0.40 }}; k++;
    p.e[k] = (r2_t){{ 0.55, 0.45 }}; k++;
    p.e[k] = (r2_t){{ 0.40, 0.50 }}; k++;
    assert(k == n);
    r2_vec_trim(&p, n);
    for (int32_t i = 0; i < n; i++) { txhg_scale_point(cvs, &(p.e[i])); }
    double r = 0.035; txhg_scale_length(cvs, &(r));
    hxg_paint_sausage(cvs, &p, r, op);
  }
