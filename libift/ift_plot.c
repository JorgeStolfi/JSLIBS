/* Implementation of ift_plot.h. */
/* Last edited on 2023-02-03 22:22:14 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <epswr.h>
#include <frgb.h>
#include <affirm.h>
#include <bool.h>

#include <ift.h>
#include <ift_plot.h>

void ift_plot_pixel
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col, 
    ift_pixel_index_t row,
    frgb_t *rgb,
    int outline
  )
  {
    double xclo = (double)col;
    double yclo = (double)(G->rows - 1 - row);
    double xchi = xclo + 1;
    double ychi = yclo + 1;
    
    epswr_set_fill_color(eps, rgb->c[0], rgb->c[1], rgb->c[2]);
    epswr_rectangle(eps, xclo,xchi, yclo,ychi, TRUE, outline);
  }

void ift_plot_node
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col, 
    ift_pixel_index_t row,
    double radius,
    frgb_t *rgb
  )
  { double xp = 0.5+(double)col;
    double yp = 0.5+(double)(G->rows - 1 - row);
    epswr_set_fill_color(eps, rgb->c[0], rgb->c[1], rgb->c[2]);
    epswr_dot(eps, xp,yp, radius, TRUE, TRUE);
  }

void ift_plot_arc
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col1, 
    ift_pixel_index_t row1,
    ift_pixel_index_t col2, 
    ift_pixel_index_t row2,
    int arrow
  );
  /* Draws an arc from pixel (col1, row1) to pixel (col2, row2).
    If {arrow} is TRUE, also draws the arrowhead. */

void ift_plot_pixel_values
  ( epswr_figure_t *eps, 
    ift_graph_t *G,
    frgb_t rgb[],
    double whiten
  )
  { int k;
    for (k = 0; k < G->nodes; k++)
      { ift_node_t *p = &(G->node[k]);
        frgb_t yy = rgb[k];
        if (whiten != 0.0)
          { yy.c[0] = (float)(whiten + (1-whiten)*yy.c[0]);
            yy.c[1] = (float)(whiten + (1-whiten)*yy.c[1]);
            yy.c[2] = (float)(whiten + (1-whiten)*yy.c[2]);
          }
        ift_plot_pixel(eps, G, p->col, p->row, &yy, 0);
      }
  }

void ift_plot_forest_edges
  ( epswr_figure_t *eps, 
    ift_graph_t *G
  )
  { int k;
    for (k = 0; k < G->nodes; k++)
      { ift_node_t *p = &(G->node[k]);
        ift_node_t *q = p->P;
        if ((q != NULL) && (q != p))
          { double xp = 0.5+(double)p->col;
            double yp = 0.5+(double)(G->rows - 1 - p->row);
            double xq = 0.5+(double)q->col;
            double yq = 0.5+(double)(G->rows - 1 - q->row);
            epswr_segment(eps, xp,yp, xq,yq);
          }
      }
  }

void ift_plot_forest_nodes
  ( epswr_figure_t *eps, 
    ift_graph_t *G,
    bool_t roots,
    double radius,
    frgb_t *rgb
  )
  { int k;
    for (k = 0; k < G->nodes; k++)
      { ift_node_t *p = &(G->node[k]);
        ift_node_t *q = p->P;
        bool_t isroot = ((q == NULL) | (q == p));
        if (isroot == roots)
          { ift_plot_node(eps, G, p->col, p->row, radius, rgb); }
      }
  }

