/* See {delaunay_plot.h}. */
/* Last edited on 2016-04-01 01:04:46 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <quad.h>
#include <pswr.h>
#include <bool.h>
#include <i3.h>
#include <r2.h>

#include <delaunay.h>
#include <delaunay_plot.h>

PSStream *open_ps_stream(double r, char *prefix, bool_t eps);
void close_ps_stream(PSStream *fps, bool_t eps);
r2_t far_left(delaunay_site_t *ap, delaunay_site_t *bp);
r2_t left_voronoi_vertex(quad_arc_t e);
double max_site_radius(delaunay_site_t *st, int nsites);
r2_t circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp);

#define RBIG 10.0

void draw_voronoi_edges(PSStream *fps, quad_arc_t e);
void draw_delaunay_edges(PSStream *fps, quad_arc_t e);
void draw_sites(PSStream *fps, delaunay_site_t *st, int nsites);

void plot_delaunay (quad_arc_t e, delaunay_site_t *st, int nsites, char *prefix, bool_t eps)
  { 
    double r = max_site_radius(st, nsites);
  
    /* Create Postscript document or EPS figure stream. */
    PSStream *fps = open_ps_stream(r, prefix, eps);
    
    /* Start a new picture: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;
    pswr_new_picture(fps, xmin,xmax, ymin, ymax);
    
    /* Plot Voronoi edges, dashed: */    
    float penwd = (eps ? 0.20f : 0.10f);
    pswr_set_pen(fps, 0,0,0, penwd, 1.0, 1.0);
    draw_voronoi_edges(fps, e);
    
    /* Plot Delaunay edges, solid: */
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    draw_delaunay_edges(fps, e);
    
    /* Plot sites: */
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    pswr_set_fill_color(fps, 1.00,0.00,0.75);
    draw_sites(fps, st, nsites);

    if (! eps)
      { /* Add caption and frame: */
        pswr_set_pen(fps, 0,0,0, 0.10, 0.0, 0.0);
        pswr_add_caption(fps, "Voronoi/Delaunay diagram", 0.5);
        pswr_set_pen(fps, 0,0,0, 0.20, 0.0, 0.0);
        pswr_frame(fps);
      }
    /* We are done: */
    pswr_close_stream(fps);
  }

void draw_voronoi_edges(PSStream *fps, quad_arc_t e)
  {
    auto void draw_voronoi_edge(quad_arc_t e);
    
    quad_arc_vec_t root = quad_arc_vec_make_desc(&e, 1);
    quad_enum(&root, draw_voronoi_edge);

    void draw_voronoi_edge(quad_arc_t e)
      { r2_t lp = left_voronoi_vertex(e);
        r2_t rp = left_voronoi_vertex(quad_sym(e));
        pswr_segment(fps, lp.c[0], lp.c[1], rp.c[0], rp.c[1]);
      }
   }

void draw_delaunay_edges(PSStream *fps, quad_arc_t e)
  {
    auto void draw_delaunay_edge(quad_arc_t e);
    
    quad_arc_vec_t root = quad_arc_vec_make_desc(&e, 1);
    quad_enum(&root, draw_delaunay_edge);
    
    void draw_delaunay_edge(quad_arc_t e)
      { delaunay_site_t *op = ORG(e);
        delaunay_site_t *dp = DST(e);
        r2_t rop = delaunay_r2_from_hi2(&(op->pt));
        r2_t rdp = delaunay_r2_from_hi2(&(dp->pt));
        pswr_segment(fps, rop.c[0], rop.c[1], rdp.c[0], rdp.c[1]);
      }
 }

void draw_sites(PSStream *fps, delaunay_site_t st[], int nsites)
  {
    int i;
    for (i = 0; i < nsites; i++) 
      { delaunay_site_t *si = &(st[i]);
        r2_t pi = delaunay_r2_from_hi2(&(si->pt));
        pswr_dot(fps, pi.c[0], pi.c[1], 0.5,  TRUE, TRUE); 
      }
  }

  
r2_t left_voronoi_vertex(quad_arc_t e)
  { delaunay_site_t *ap =  ORG(e);
    delaunay_site_t *bp =  DST(e);
    delaunay_site_t *cp =  DST(quad_lnext(e));
    if (delaunay_orient(ap, bp, cp) > 0)
      { /* Internal face, compute circumcenter: */
        return circumcenter(ap, bp, cp);
      }
    else
      { /* External face, compute an "almost infinte" point on the left: */
        return far_left(ap, bp);
      }
  }

r2_t circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp)
  { 
    r2_t rap = delaunay_r2_from_hi2(&(ap->pt));
    r2_t rbp = delaunay_r2_from_hi2(&(bp->pt));
    r2_t rcp = delaunay_r2_from_hi2(&(cp->pt));
    return r2_circumcenter(&rap, &rbp, &rcp);
  }
  
double max_site_radius(delaunay_site_t *st, int nsites)
  { double r2max = 0.0;
    int i;
    for (i = 0; i < nsites; i++) 
      { delaunay_site_t *si = &(st[i]);
        r2_t pi = delaunay_r2_from_hi2(&(si->pt));
        double r2i = r2_norm_sqr(&pi);
        if (r2i > r2max) { r2max = r2i; }
      }
    return sqrt(r2max);
  }

r2_t far_left(delaunay_site_t *ap, delaunay_site_t *bp)
  { r2_t a = delaunay_r2_from_hi2(&(ap->pt));
    r2_t b = delaunay_r2_from_hi2(&(bp->pt));
    r2_t m, d, v; 
    r2_mix(0.5, &a, 0.5, &b, &m);
    r2_sub(&b, &a, &d);
    (void)r2_dir(&d, &d);
    v.c[0] = m.c[0] - RBIG * d.c[1];
    v.c[1] = m.c[1] + RBIG * d.c[0];
    return v;
  }
  
PSStream *open_ps_stream(double r, char *prefix, bool_t eps)
  { 
    double mm = (72.0/25.4); /* One mm in pt. */
    double xfigsz = 150.00*mm; /* Figure X size excluding margin (pt). */
    double yfigsz = 150.00*mm; /* Figure Y size excluding margin (pt). */
    double fmrg = 3.0; /* Figure margin width (pt). */
    double pmrg = 2.0; /* Picture margin width (pt). */

    /* Add caption only if there is a user caption, or it is not EPS. */
    /* Select a good figure size: */
    PSStream *fps = pswr_new_stream
      ( /* prefix */                txtcat(prefix, "-"),
        /* file */                  NULL,
        /* eps */                   eps,
        /* docName */               "doc",
        /* paperSize */             "letter",
        /* landscape */             FALSE,
        /* hPageSize, vPageSize */  xfigsz + 2*fmrg, yfigsz + 2*fmrg
      );
    pswr_set_canvas_layout
      ( fps,
        /* hPicSize, vPicSize */     xfigsz, yfigsz,
        /* adjustPicSize */          FALSE,
        /* hPicMargin,vPicMargin */  pmrg, pmrg,
        /* captionLines */           (eps ? 0 : 1),  
        /* vCount, hCount */         0, 0  /* Let {pswr} choose it. */
      ); 
    return fps;
  }  
