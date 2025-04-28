/* See {psgr_plot.h} */
/* Last edited on 2025-04-27 02:59:29 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <bool.h>
#include <epswr.h>
#include <interval.h>
#include <box.h>
#include <i2.h>
#include <r2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_path.h>

#include <psgr_plot.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

void psgr_plot(epswr_figure_t *eps, psgr_t *gr, double pixelSize, double fontSize, double vertexRadius)
  {
    uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);
    
    auto void plot_vertex(psgr_vertex_t vi);

    auto void plot_edge(psgr_edge_t ei);

    if (fontSize > 0)
      { epswr_set_label_font(eps, "Courier", fontSize); }

    /* Plot the edges: */
    for (uint32_t e_ix = 0; e_ix < NE; e_ix++)
      { psgr_edge_t e = (psgr_edge_t){ e_ix };
        plot_edge(e);
      }

    /* Plot the vertices: */
    for (uint32_t v_ix = 0; v_ix < NV; v_ix++)
      { psgr_vertex_t v = (psgr_vertex_t){ v_ix };
        plot_vertex(v);
     }

    void plot_vertex(psgr_vertex_t v)
      { r2_t p = psgr_plot_coords_from_pixel(psgr_vertex_pixel(gr, v), pixelSize);
        epswr_set_pen(eps, 0,0,0, 0.1, 0, 0);
        bool_t cFill = TRUE, cDraw = TRUE;
        switch (psgr_vertex_mark(gr, v))
          { case CLEAR:   epswr_set_fill_color(eps, 1.000,1.000,1.000); break;
            case KEEP:    epswr_set_fill_color(eps, 0.000,0.800,0.200); break;
            case KILL:    epswr_set_fill_color(eps, 1.000,0.000,0.000); break;
            case DELETED: epswr_set_fill_color(eps, 0.400,0.300,0.300); break;
            default: demand(FALSE, "invalid vertex mark");
          }
        epswr_circle(eps, p.c[0],p.c[1], vertexRadius, cFill, cDraw);

        if (fontSize > 0)
          { char *label = jsprintf("%ld", v);
            double rot = 0.0, hAlign = -0.5, vAlign = 0.5;
            bool_t clipped = FALSE; /* Labels may extend outside plot window. */
            bool_t tFill = TRUE, tDraw = FALSE;
            epswr_set_fill_color(eps, 1,0,0);
            epswr_label
              ( eps, label, "0", p.c[0], p.c[1], 
                rot, clipped, hAlign, vAlign, tFill, tDraw
              );
            free(label);
          }
      }

    void plot_edge(psgr_edge_t e)
      {
        psgr_arc_t a = psgr_orient_edge(e, 0);
        
        psgr_vertex_t org = psgr_arc_org(gr, a);
        psgr_vertex_t dst = psgr_arc_dst(gr, a);
        
        i2_t ipo = psgr_vertex_pixel(gr, org);
        i2_t ipd = psgr_vertex_pixel(gr, dst);
         
        epswr_set_pen_width(eps, 0.1);
        epswr_set_pen_dashing(eps, 0,0);
        switch (psgr_edge_mark(gr, e))
          { case CLEAR:     epswr_set_pen_color(eps, 0.000,0.000,0.000); break;
            case DELETED:   epswr_set_pen_color(eps, 0.800,0.800,0.800); break;
            default: demand(FALSE, "invalid edge mark");
          }
        psgr_path_t P = psgr_arc_path(gr, a);
        psgr_plot_path(eps, pixelSize, &ipo, P, &ipd);

        if (fontSize > 0) 
          { char *label = jsprintf("%ld", e);
            r2_t ctr = psgr_path_center(&ipo, P, &ipd);
            r2_scale(pixelSize, &ctr, &ctr);
            r2_t dir = psgr_path_central_dir(&ipo, P, &ipd);
            double rot = 0.0;
            bool_t clipped = FALSE; /* Labels may extend outside plot window. */
            double hAlign = (dir.c[0] > 0 ? -0.2 : +1.2);
            double vAlign = (dir.c[1] > 0 ? +1.2 : -0.2);
            bool_t tFill = TRUE, tDraw = FALSE;
            epswr_set_fill_color(eps, 0,0,0);
            epswr_label
              ( eps, label, "0", ctr.c[0],ctr.c[1], 
                rot, clipped, hAlign, vAlign, tFill, tDraw
              );
            free(label);
          }
      }
  }
  
void psgr_plot_named(char *fname, psgr_t *gr, double pixelSize, double fontSize, double vertexRadius)
  {
    uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);

    /* Get bounding box (mm) of vertex centers: */
    interval_t bbox[2];
    box_empty(2, bbox);
    for (uint32_t v_ix = 0; v_ix < NV; v_ix++)
      { psgr_vertex_t v = (psgr_vertex_t){ v_ix };
        r2_t p = psgr_plot_coords_from_pixel(psgr_vertex_pixel(gr, v), pixelSize);
        box_include_point(2, bbox, p.c, bbox);
      }
    
    /* Enlarge bbox for vertex radius plusline width: */
    double mrg_vrad = vertexRadius + 1.0; /* Vertex radius and line width. */
    box_widen(2, bbox, mrg_vrad, bbox);
     
    /* Margins to account for vertex radii, line widths, and labels: */
    double fsmm = fontSize/epswr_pt_per_mm;
    
    double mrg_vlab_h = (NV == 0 ? 0 : fsmm*(1 + digits(NV-1))); /* Vertex labels (horz). */
    double mrg_elab_h = (NE == 0 ? 0 : fsmm*(1 + digits(NE-1))); /* Edge labels (horz). */
    double mrg_h = fmax(0, fmax(mrg_vlab_h, mrg_elab_h) - mrg_vrad);
    
    double mrg_vlab_v = (NV == 0 ? 0 : fsmm); /* Vertex labels (vert). */
    double mrg_elab_v = (NE == 0 ? 0 : fsmm); /* Edge labels (vert). */
    double mrg_v = fmax(0, fmax(mrg_vlab_v, mrg_elab_v) - mrg_vrad);
    
    double hPlotSize = epswr_pt_per_mm * (HI(bbox[0]) - LO(bbox[0]));
    double vPlotSize = epswr_pt_per_mm * (HI(bbox[1]) - LO(bbox[1]));
    double hMargin = epswr_pt_per_mm * mrg_h;
    double vMargin = epswr_pt_per_mm * mrg_v;

    FILE *wr = open_write(fname, TRUE);
    epswr_figure_t *eps = epswr_new_figure(wr, hPlotSize, vPlotSize, hMargin, hMargin, vMargin, vMargin, FALSE);
    epswr_set_client_window(eps, LO(bbox[0]), HI(bbox[0]), LO(bbox[1]), HI(bbox[1]));
    psgr_plot(eps, gr, pixelSize, fontSize, vertexRadius);
    epswr_end_figure(eps);
  }
  
void psgr_plot_path(epswr_figure_t *eps, double pixelSize, i2_t *ipo, psgr_path_t P, i2_t *ipd)
  { if (P.reverse) { i2_t *ipt = ipo; ipo = ipd; ipd = ipt; }
    
    r2_t p = psgr_plot_coords_from_pixel(*ipo, pixelSize);
    for (int32_t k = 0; k <= P.n; k++)
      { r2_t q = psgr_plot_coords_from_path(k, ipo, P, ipd, pixelSize);
        epswr_segment(eps, p.c[0],p.c[1], q.c[0],q.c[1]);
        p = q;
      }
    return;
  }

r2_t psgr_plot_coords_from_pixel(i2_t ip, double pixelSize)
  { demand((ip.c[0] >= 0) && (ip.c[1] >= 0), "negative pixel indices");
    r2_t p = (r2_t){{ pixelSize*ip.c[0], pixelSize*ip.c[1] }};
    return p;
  }

r2_t psgr_plot_coords_from_path(int32_t k, i2_t *ipo, psgr_path_t P, i2_t *ipd, double pixelSize)
  {
    r2_t p ;
    if (k == -1)
      { p = (r2_t){{ ipo->c[0], ipo->c[1] }}; }
    else if (k == P.n)
      { p = (r2_t){{ ipd->c[0], ipd->c[1] }}; }
    else
      { int32_t kr = (P.reverse ? (int32_t)P.n - k : k);
        assert((kr >= 0) && (kr <= P.n));
        p = P.v[kr];
      }
    r2_scale(pixelSize, &p, &p);
    return p;
  }
