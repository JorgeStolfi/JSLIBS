/* See {sheet_cut_plot.h} */
/* Last edited on 2023-02-04 09:05:35 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <r2.h>
#include <bool.h>
#include <epswr.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsfile.h>

#include <sheet_cut.h>
#include <sheet_cut_plot.h>

/* IMPLEMENTATIONS */  

epswr_figure_t *sheet_cut_plot_new_figure
  ( char *outPrefix, 
    char *sheet_tag,
    char *sheet_mat, 
    double sheet_thk,
    r2_t sheet_dim,
    double scrap_mrg,
    double cur_px,
    double cur_py,
    double cur_ht
  )
  { 
    double mrg = sheet_cut_plot_FIG_MARGIN_PT;
    /* Plot window size in device coordinates (pt): */
    double hPlotSize = sheet_cut_plot_SCALE*sheet_dim.c[0];
    double vPlotSize = sheet_cut_plot_SCALE*sheet_dim.c[1];
    
    double topMrg = mrg + sheet_cut_plot_FONT_SIZE_PT + mrg;
    bool_t verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, outPrefix, sheet_tag, -1, NULL,
        hPlotSize, vPlotSize, mrg, mrg, mrg, topMrg, 
        verbose
      );
    
    /* Compute the size of the sheet's usable area (plot window) in mm: */
    r2_t scrapVec = (r2_t){{ 2*scrap_mrg, 2*scrap_mrg }};
    r2_t usableSize; r2_sub(&sheet_dim, &scrapVec, &usableSize);

    /* Set the font for all labels: */
    epswr_set_label_font(eps, sheet_cut_plot_FONT_NAME, sheet_cut_plot_FONT_SIZE_PT);

    /* Write the sheet info text at top: */
    char *label = NULL;
    asprintf
      ( &label, 
        "sheet %s: %s %.1f x %.1f x %.1f mm (usable %.1f x %.1f mm)", 
        sheet_tag, sheet_mat, sheet_thk, 
        sheet_dim.c[0], sheet_dim.c[1], 
        usableSize.c[0], usableSize.c[1]
      );
    epswr_label(eps, label, "X", mrg, mrg + vPlotSize + mrg, 0.0, FALSE, 0.0,0.0, TRUE, FALSE);
    free(label);
    
    /* Set the plot window (in mm) to the physical sheet: */
    double xMin = 0.0, xMax = sheet_dim.c[0];
    double yMin = 0.0, yMax = sheet_dim.c[1];
    double hMin = mrg, hMax = mrg + hPlotSize;
    double vMin = mrg, vMax = mrg + vPlotSize;
    epswr_set_window(eps, hMin, hMax, vMin, vMax, FALSE, xMin, xMax, yMin, yMax);
    
    /* Draw the outline of the physical sheet in black: */
    epswr_set_pen(eps,  0.000,0.000,0.000,  0.50, 0,0);
    epswr_frame(eps);
    
    /* Reset the plot window (in mm) to the usable area: */
    r2_t size = (r2_t){{ sheet_dim.c[0] - 2*scrap_mrg, sheet_dim.c[1] - 2*scrap_mrg }};
    xMin = 0.0, xMax = size.c[0];
    yMin = 0.0, yMax = size.c[1];
    double ptScrap = sheet_cut_plot_SCALE*scrap_mrg;
    fprintf(stderr, "scrap margin = %.2f pt\n", ptScrap);
    hMin = mrg + ptScrap, hMax = mrg + hPlotSize - ptScrap;
    vMin = mrg + ptScrap, vMax = mrg + vPlotSize - ptScrap;
    epswr_set_window(eps, hMin, hMax, vMin, vMax, FALSE, xMin, xMax, yMin, yMax);

    if (! (isnan(cur_px) || isnan(cur_py) || isnan(cur_ht)))
      { /* Draw the staircase outline: */
        epswr_set_pen(eps,  0.000,0.000,1.000,  0.30, 0,0);
        epswr_segment(eps, xMin,yMin, xMin,cur_ht);
        epswr_segment(eps, xMin,cur_ht, cur_px,cur_ht);
        epswr_segment(eps, cur_px,cur_ht, cur_px,cur_py);
        epswr_segment(eps, cur_px,cur_py, xMax,cur_py);
        epswr_segment(eps, xMax,cur_py, xMax,yMin);
        epswr_segment(eps, xMax,yMin, xMin,yMin);
      }
      
    return eps;
  }
  
void sheet_cut_plot_all_nodes(epswr_figure_t *eps, sheet_cut_node_t *pc, r2_t org_pc)
  {
    auto void visit(sheet_cut_node_t *pr, r2_t org_pr, bool_t post);
      /* Visiting procedure for {sheet_cut_enum}. 
      Will be called for each plate or block {pr} reachable from {pc}.
      
      If {pr} is a block, on the first visit ({post=FALSE})
      draws the bounding box in red, and does nothing in the second visit.
      
      If it is a plate, just draws it in black, labeled with 
      its tag and other info. */
      
    sheet_cut_enum(pc, org_pc, visit);
    return;
    
    /* Inpternal implementations: */
      
    void visit(sheet_cut_node_t *pr, r2_t org_pr, bool_t post)
      { r2_t low_pr; r2_add(&org_pr, &(pr->pos), &low_pr); /* Plot coords of {pr}'s low corner. */
        if (pr->nsub == 0)
          { /* Plate node; just draw it: */
            sheet_cut_plot_plate(eps, low_pr, pr->size, pr->tag);
          }
        else if (! post)
          { /* First visit to block. Just draw the bounding box : */
            sheet_cut_plot_block_bbox(eps, low_pr, pr->size);
          }
      }
  }

void sheet_cut_plot_block_bbox(epswr_figure_t *eps, r2_t pos, r2_t size)
  {
    /* Get the bounding box position and size: */
    double xPos = pos.c[0], yPos = pos.c[1];
    double xDim = size.c[0], yDim = size.c[1];

    /* Draw it: */
    epswr_set_pen(eps,  8.000,0.600,0.200,  0.25, 0,0);
    epswr_rectangle(eps, xPos, xPos + xDim, yPos, yPos + yDim, FALSE, TRUE);
  }

void sheet_cut_plot_plate(epswr_figure_t *eps, r2_t pos, r2_t size, char *tag)
  { 
    /* Draw the rectangle: */
    double xPos = pos.c[0], yPos = pos.c[1];
    double xDim = size.c[0], yDim = size.c[1];
    epswr_set_pen(eps,  0.000,0.000,0.000,  0.25, 0,0);
    epswr_set_fill_color(eps, 1.000, 1.000, 0.850);
    epswr_rectangle(eps, xPos, xPos + xDim, yPos, yPos + yDim, TRUE, TRUE);

    /* Write the text with dimensions and tag: */
    sheet_cut_plot_labels(eps, pos, size, tag);
  }

void sheet_cut_plot_labels(epswr_figure_t *eps, r2_t pos, r2_t size, char *tag)
  {
    double mrg = sheet_cut_plot_FIG_MARGIN_PT/sheet_cut_plot_SCALE; /* Margin space width in mm */
    double vfs = sheet_cut_plot_FONT_SIZE_PT/sheet_cut_plot_SCALE;  /* Vertical font size in mm */
    
    bool_t rotate = (size.c[0] < size.c[1]); /* Should labels be rotated? */
    double rot = (rotate ? -90.0 : 0.0); /* CCW rotation angle. */
    
    epswr_set_fill_color(eps,  0.000,0.000,0.000);
    epswr_set_pen(eps,  0.000,0.000,0.000,  0.25, 0,0);

    /* Plot coordinates at each corner: */
    for (int32_t px = 0; px < 2; px++)
      { for (int32_t py = 0; py < 2; py++)
          { /* Coordinates of corner: */
            r2_t cpos = (r2_t){{ pos.c[0] + px*size.c[0], pos.c[1] + py*size.c[1] }}; 
            /* Coordinates or label ref point and text rotation angle: */
            r2_t tpos = (r2_t){{ cpos.c[0] + mrg*(1-2*px), cpos.c[1] + mrg*(1-2*py) }};
            /* Label alignment (0.0 = bottom/left, 1.0 = top/right). */
            double xAlign, yAlign;
            if (rotate)
              { xAlign = (double)(1-py); yAlign = (double)px; }
            else
              { xAlign = (double)px; yAlign = (double)py; }
            /* Generate the label text: */
            char *txt = NULL;
            asprintf(&txt, "(%.1f, %.1f)", cpos.c[0], cpos.c[1]);
            epswr_label(eps, txt, "0", tpos.c[0], tpos.c[1], rot, FALSE, xAlign, yAlign, TRUE, FALSE);
            free(txt);
          }
      }
     
    /* Plot plate tag and size at center: */
    r2_t ctr = (r2_t){{ pos.c[0] + 0.5*size.c[0], pos.c[1] + 0.5*size.c[1] }};
    r2_t vst; /* Displacement vector for each additional line of text. */
    if (rotate)
      { vst = (r2_t){{ mrg + vfs, 0 }}; }
    else
      { vst = (r2_t){{ 0, -(mrg + vfs) }}; }
    
    /* Plate tag is on first line: */
    r2_t tpos = (r2_t){{ ctr.c[0] - 0.5*vst.c[0], ctr.c[1] - 0.5*vst.c[1] }};
    epswr_label(eps, tag, "X", tpos.c[0], tpos.c[1], rot, FALSE, 0.5, 0.5, TRUE, FALSE);

    /* Plate dimensions on second line: */
    char *txt = NULL;
    asprintf(&txt, "%.1f x %.1f", size.c[0], size.c[1]);
    r2_add(&tpos, &vst, &tpos);
    epswr_label(eps, txt, "0", tpos.c[0], tpos.c[1], rot, FALSE, 0.5, 0.5, TRUE, FALSE);
    free(txt);
  }

void sheet_cut_plot_end_figure(epswr_figure_t *eps)
  { 
    epswr_end_figure(eps);
  }

