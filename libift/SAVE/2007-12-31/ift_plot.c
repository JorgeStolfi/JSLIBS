/* Implementation of ift_plot.h. */
/* Last edited on 2007-12-26 16:23:10 by stolfi */

#define _GNU_SOURCE
#include <pswr.h>
#include "ift.h"
#include "ift_plot.h"

void ift_plot_pixel(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col, 
    PixelIndex row,
    double r, double g, double b,
    int outline
  )
  {
    double xclo = (double)col;
    double yclo = (double)(G->rows - 1 - row);
    double xchi = xclo + 1;
    double ychi = yclo + 1;
    
    pswr_set_fill_color(ps, r, g, b);
    pswr_rectangle(ps, xclo,xchi, yclo,ychi, TRUE, outline);
  }

void ift_plot_node(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col, 
    PixelIndex row,
    double radius,
    double r, double g, double b
  )
  { double xp = 0.5+(double)col;
    double yp = 0.5+(double)(G->rows - 1 - row);
    pswr_set_fill_color(ps, r, g, b);
    pswr_dot(ps, xp,yp, radius, TRUE, TRUE);
  }

void ift_plot_arc(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col1, 
    PixelIndex row1,
    PixelIndex col2, 
    PixelIndex row2,
    int arrow
  );
  /* Draws an arc from pixel (col1, row1) to pixel (col2, row2).
    If `arrow' is TRUE, also draws the arrowhead. */

void ift_plot_pixel_values(
    PSStream *ps, 
    ImageGraph *G,
    double whiten
  )
  { int k;
    for (k = 0; k < G->nodes; k++)
      { PixelNode *p = &(G->node[k]);
        double r = 0, g = 0, b = 0;
        if (G->channels == 1)
          { double yy = p->y.c[0];
            r = yy; g = yy; b = yy;
          }
        else if (G->channels == 1)
          { r = p->y.c[0];
            g = p->y.c[1];
            b = p->y.c[2];
          }
        else
          { IFT_ERROR("bad number of channels"); }
        if (whiten != 0.0)
        { r = whiten + (1-whiten)*r;
          g = whiten + (1-whiten)*g;
          b = whiten + (1-whiten)*b;
        }
        ift_plot_pixel(ps, G, p->col, p->row, r, g, b, 0);
      }
  }

void ift_plot_forest_edges(
    PSStream *ps, 
    ImageGraph *G
  )
  { int k;
    for (k = 0; k < G->nodes; k++)
      { PixelNode *p = &(G->node[k]);
        PixelNode *q = p->P;
        if ((q != NULL) && (q != p))
          { double xp = 0.5+(double)p->col;
            double yp = 0.5+(double)(G->rows - 1 - p->row);
            double xq = 0.5+(double)q->col;
            double yq = 0.5+(double)(G->rows - 1 - q->row);
            pswr_segment(ps, xp,yp, xq,yq);
          }
      }
  }

void ift_plot_forest_nodes(
    PSStream *ps, 
    ImageGraph *G,
    int roots,
    double radius,
    double r, double g, double b
  )
  { int isroot, k;
    for (k = 0; k < G->nodes; k++)
      { PixelNode *p = &(G->node[k]);
        PixelNode *q = p->P;
        isroot = ((q == NULL) | (q == p));
        if (isroot == roots)
          { ift_plot_node(ps, G, p->col, p->row, radius, r, g, b); }
      }
  }

