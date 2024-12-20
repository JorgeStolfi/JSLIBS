/* See {delaunay_plot.h}. */
/* Last edited on 2024-12-05 10:25:06 by stolfi */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <jsstring.h>
#include <quad.h>
#include <epswr.h>
#include <bool.h>
#include <i3.h>
#include <r2.h>

#include <delaunay.h>
#include <delaunay_plot.h>

r2_t delaunay_plot_far_left(delaunay_site_t *ap, delaunay_site_t *bp);
r2_t delaynay_plot_left_voronoi_vertex(quad_arc_t e);
double delaunay_plot_max_site_radius(delaunay_site_t *st, int32_t nsites);
r2_t delaunay_plot_circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp);

#define RBIG 10.0

void delaunay_plot_delaunay_plot_draw_voronoi_edges(epswr_figure_t *eps, quad_arc_t e);
void delaunay_plot_draw_edges(epswr_figure_t *eps, quad_arc_t e);
void delaunay_plot_draw_sites(epswr_figure_t *eps, delaunay_site_t *st, int32_t nsites);

void delaunay_plot_diagram(epswr_figure_t *eps, quad_arc_t e, delaunay_site_t *st, int32_t nsites)
  { 
    /* Plot Voronoi edges, dashed: */    
    epswr_set_pen(eps, 0,0,0, 0.20f, 1.0, 1.0);
    delaunay_plot_delaunay_plot_draw_voronoi_edges(eps, e);
    
    /* Plot Delaunay edges, solid: */
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    delaunay_plot_draw_edges(eps, e);
    
    /* Plot sites: */
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    epswr_set_fill_color(eps, 1.00,0.00,0.75);
    delaunay_plot_draw_sites(eps, st, nsites);

    /* Add frame: */
    epswr_set_pen(eps, 0,0,0, 0.10f, 0.0, 0.0);
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    epswr_frame(eps);
  }

void delaunay_plot_delaunay_plot_draw_voronoi_edges(epswr_figure_t *eps, quad_arc_t e)
  {
    auto void delaunay_plot_draw_voronoi_edge(quad_arc_t e);
    
    quad_arc_vec_t root = quad_arc_vec_make_desc(&e, 1);
    quad_enum(&root, delaunay_plot_draw_voronoi_edge);

    void delaunay_plot_draw_voronoi_edge(quad_arc_t e)
      { r2_t lp = delaynay_plot_left_voronoi_vertex(e);
        r2_t rp = delaynay_plot_left_voronoi_vertex(quad_sym(e));
        epswr_segment(eps, lp.c[0], lp.c[1], rp.c[0], rp.c[1]);
      }
   }

void delaunay_plot_draw_edges(epswr_figure_t *eps, quad_arc_t e)
  {
    auto void delaunay_plot_draw_edge(quad_arc_t e);
    
    quad_arc_vec_t root = quad_arc_vec_make_desc(&e, 1);
    quad_enum(&root, delaunay_plot_draw_edge);
    
    void delaunay_plot_draw_edge(quad_arc_t e)
      { delaunay_site_t *op = ORG(e);
        delaunay_site_t *dp = DST(e);
        r2_t rop = delaunay_r2_from_hi2(&(op->pt));
        r2_t rdp = delaunay_r2_from_hi2(&(dp->pt));
        epswr_segment(eps, rop.c[0], rop.c[1], rdp.c[0], rdp.c[1]);
      }
 }

void delaunay_plot_draw_sites(epswr_figure_t *eps, delaunay_site_t st[], int32_t nsites)
  {
    int32_t i;
    for (i = 0; i < nsites; i++) 
      { delaunay_site_t *si = &(st[i]);
        r2_t pi = delaunay_r2_from_hi2(&(si->pt));
        epswr_dot(eps, pi.c[0], pi.c[1], 0.5,  TRUE, TRUE); 
      }
  }
  
r2_t delaynay_plot_left_voronoi_vertex(quad_arc_t e)
  { delaunay_site_t *ap =  ORG(e);
    delaunay_site_t *bp =  DST(e);
    delaunay_site_t *cp =  DST(quad_lnext(e));
    if (delaunay_orient(ap, bp, cp) > 0)
      { /* Internal face, compute delaunay_plot_circumcenter: */
        return delaunay_plot_circumcenter(ap, bp, cp);
      }
    else
      { /* External face, compute an "almost infinte" point on the left: */
        return delaunay_plot_far_left(ap, bp);
      }
  }

r2_t delaunay_plot_circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp)
  { 
    r2_t rap = delaunay_r2_from_hi2(&(ap->pt));
    r2_t rbp = delaunay_r2_from_hi2(&(bp->pt));
    r2_t rcp = delaunay_r2_from_hi2(&(cp->pt));
    return r2_circumcenter(&rap, &rbp, &rcp);
  }
  
double delaunay_plot_max_site_radius(delaunay_site_t *st, int32_t nsites)
  { double r2max = 0.0;
    int32_t i;
    for (i = 0; i < nsites; i++) 
      { delaunay_site_t *si = &(st[i]);
        r2_t pi = delaunay_r2_from_hi2(&(si->pt));
        double r2i = r2_norm_sqr(&pi);
        if (r2i > r2max) { r2max = r2i; }
      }
    return sqrt(r2max);
  }

r2_t delaunay_plot_far_left(delaunay_site_t *ap, delaunay_site_t *bp)
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
  
epswr_figure_t *delaunay_plot_new_figure
  ( char *prefix, 
    char *tag,
    int32_t capLines, 
    delaunay_site_t st[], 
    int32_t nsites
  )
  { 
    double r = delaunay_plot_max_site_radius(st, nsites);
  
    /* Set the client window: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;

    double fontHeight = 10;
    double maxSize = 150*epswr_pt_per_mm;
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_captioned_figure
      ( NULL, prefix, NULL, nsites, tag,
        xmin,xmax, ymin,ymax, 
        maxSize, maxSize, capLines, fontHeight,
        eps_verbose
      );
    epswr_set_fill_color(eps, 0,0,0);
    return eps;
  }  
