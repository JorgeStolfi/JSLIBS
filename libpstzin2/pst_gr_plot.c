/* See {pst_gr_plot.h} */
/* Last edited on 2025-03-14 06:47:55 by stolfi */
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

#include <pst_gr.h>
#include <pst_gr_path.h>

#include <pst_gr_plot.h>

void pst_gr_plot(epswr_figure_t *eps, pst_gr_t* gr, double fontSize, double vertexRadius)
  {
    auto void plot_vertex(pst_gr_vertex_t vi, pst_gr_vertex_data_t *vd);

    auto void plot_edge(pst_gr_edge_t ei, pst_gr_edge_data_t *ed);

    if (fontSize > 0)
      { epswr_set_label_font(eps, "Courier", fontSize); }

    /* Plot the edges: */
    for (pst_gr_edge_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_edge_data_t *ed = &(gr->edata[ei]);
        plot_edge(ei, ed);
      }

    /* Plot the vertices: */
    for (pst_gr_vertex_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_vertex_data_t* vd = &(gr->vdata[vi]);
        plot_vertex(vi, vd);
     }

    void plot_vertex(pst_gr_vertex_t vi, pst_gr_vertex_data_t *vd)
      { r2_t *vp = &(vd->coords);
        
        epswr_set_pen(eps, 0,0,0, 0.1, 0, 0);
        bool_t cFill = TRUE, cDraw = TRUE;
        if (vd->vmark == 0)  {epswr_set_fill_color(eps, 1.000,1.000,1.000); }
        if (vd->vmark == 1)  {epswr_set_fill_color(eps, 1.000,0.000,0.000); }
        if (vd->vmark == 2)  {epswr_set_fill_color(eps, 1.000,0.800,0.000); }
        if (vd->vmark == 3)  {epswr_set_fill_color(eps, 0.000,0.900,0.000); }
        if (vd->vmark == 4)  {epswr_set_fill_color(eps, 0.000,0.400,1.000); }
        epswr_circle(eps, vp->c[0],vp->c[1], vertexRadius, cFill, cDraw);

        if (fontSize > 0)
          { char *label = jsprintf("%ld", vi);
            double rot = 0.0, hAlign = -0.5, vAlign = 0.5;
            bool_t clipped = FALSE; /* Labels may extend outside plot window. */
            bool_t tFill = TRUE, tDraw = FALSE;
            epswr_set_fill_color(eps, 1,0,0);
            epswr_label
              ( eps, label, "0", vp->c[0],vp->c[1], 
                rot, clipped, hAlign, vAlign, tFill, tDraw
              );
            free(label);
          }
      }

    void plot_edge(pst_gr_edge_t ei, pst_gr_edge_data_t *ed)
      {
    
        pst_gr_vertex_t org = ed->org[0]; r2_t *po = &(gr->vdata[org].coords);
        pst_gr_vertex_t dst = ed->org[1]; r2_t *pd = &(gr->vdata[dst].coords);
        
        epswr_set_pen(eps, 0,0,1, 0.1, 0,0);
        pst_gr_path_plot(eps, po, ed->path, pd);

        if (fontSize > 0) 
          { char *label = jsprintf("%ld", ei);
            r2_t ctr, dir;
            pst_gr_path_ctr_dir(po, ed->path, pd, &ctr, &dir);
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
  
void pst_gr_path_plot(epswr_figure_t *eps, r2_t *po, pst_gr_path_t P, r2_t *pd)
  { r2_t *p = po;
    for (int32_t i = 0; i <= P.n; i++)
      { r2_t *q = (i == P.n ? pd : &(P.v[i]));
        epswr_segment(eps, p->c[0],p->c[1], q->c[0],q->c[1]);
        p = q;
      }
  }
  
void pst_gr_plot_named(char *fname, pst_gr_t *gr, double fontSize, double vertexRadius)
  {
    /* Get bounding box (mm) of vertex centers: */
    interval_t bbox[2];
    box_empty(2, bbox);
    for (pst_gr_vertex_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_vertex_data_t* vd = &(gr->vdata[vi]);
        r2_t *vp = &(vd->coords);
        box_include_point(2, bbox, vp->c, bbox);
      }
    
    /* Enlarge bbox for vertex radius plusline width: */
    double mrg_vrad = vertexRadius + 1.0; /* Vertex radius and line width. */
    box_widen(2, bbox, mrg_vrad, bbox);
     
    /* Margins to account for vertex radii, line widths, and labels: */
    double fsmm = fontSize/epswr_pt_per_mm;
    
    double mrg_vlab_h = (gr->NV == 0 ? 0 : fsmm*(1 + digits(gr->NV-1))); /* Vertex labels (horz). */
    double mrg_elab_h = (gr->NE == 0 ? 0 : fsmm*(1 + digits(gr->NE-1))); /* Edge labels (horz). */
    double mrg_h = fmax(0, fmax(mrg_vlab_h, mrg_elab_h) - mrg_vrad);
    
    double mrg_vlab_v = (gr->NV == 0 ? 0 : fsmm); /* Vertex labels (vert). */
    double mrg_elab_v = (gr->NE == 0 ? 0 : fsmm); /* Edge labels (vert). */
    double mrg_v = fmax(0, fmax(mrg_vlab_v, mrg_elab_v) - mrg_vrad);
    
    double hPlotSize = epswr_pt_per_mm * (HI(bbox[0]) - LO(bbox[0]));
    double vPlotSize = epswr_pt_per_mm * (HI(bbox[1]) - LO(bbox[1]));
    double hMargin = epswr_pt_per_mm * mrg_h;
    double vMargin = epswr_pt_per_mm * mrg_v;

    FILE *wr = open_write(fname, TRUE);
    epswr_figure_t *eps = epswr_new_figure(wr, hPlotSize, vPlotSize, hMargin, hMargin, vMargin, vMargin, FALSE);
    epswr_set_client_window(eps, LO(bbox[0]), HI(bbox[0]), LO(bbox[1]), HI(bbox[1]));
    pst_gr_plot(eps, gr, fontSize, vertexRadius);
    epswr_end_figure(eps);
  }
