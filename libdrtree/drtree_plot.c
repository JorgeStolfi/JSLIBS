/* See {grd_plot.h} */
/* Last edited on 2024-11-04 07:31:31 by stolfi */

#define drtree_plot_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsfile.h>
#include <vec.h> 
#include <affirm.h> 
#include <frgb.h> 
#include <frgb_ops.h> 

#include <epswr.h> 

#include <drtree.h>

#include <drtree_plot.h>

double drtree_plot_cell_center(int32_t i, int32_t n, double step);
  /* Return the coordinate of the center of column (or row) {i}
    given the number {n} of columns (or rows) and the width {step} of each 
    column (or row).  The index {i} is complemented if {step} is negative. */

/* IMPLEMENTATIONS */

void drtree_plot_individuals
  ( epswr_figure_t *eps, 
    int32_t tMin, 
    int32_t ncols, 
    int32_t nrows, 
    double Xstep, 
    double Ystep, 
    int32_t ni, 
    drtree_node_t dt[], 
    int32_t rdr[], 
    bool_t *fill,
    int32_t *chf
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "    > %s\n", __FUNCTION__); }
    
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");
    demand(ncols > 0, "invalid {ncols}");
    
    demand((Xstep != 0) && (Ystep != 0), "invalid{Xstep,Ystep}");

    int32_t tMax = tMin + ncols - 1;
    
    double lifw = 0.50*Ystep; /* Width of node traces. */

    double dbrr = 0.40*lifw;  /* Radius of birth dot. */

    double arlw = 0.33*lifw;  /* Arrow line width. */
    double alen = 4.00*arlw;  /* Arrowhead length. */
    double awid = 1.50*arlw;  /* Arrowhead width. */
    double halw = 0.30*arlw;  /* Width of arrows white halo. */
    double datr = 0.50*arlw;  /* Radius of arrow tail dot. */
    
    /* Colors for each individual: */
    frgb_t *color = NULL;
    frgb_t defColor; /* Default color for 'no family' */
    double phi = (sqrt(5) - 1)/2; /* Golden ratio. */
    
    if (chf == NULL)
      { defColor = (frgb_t){{ 0.900f, 0.600f, 0.200f }}; }
    else
      { defColor = (frgb_t){{ 0.750f, 0.720f, 0.700f }};
        color = (frgb_t*)notnull(malloc(ni*sizeof(frgb_t)), "no mem");
        double nextHue = 0.000; /* Next hue to assign. */
        for (uint32_t iq = 0;  iq < ni; iq++)
          { int32_t ic = chf[iq]; /* Chief of {q}'s family. */
            if (ic == -1)
              { /* Default color: */
                color[iq] = defColor;
              }
            else if (ic == iq)
              { /* {q} is a family chielf; assign a new color: */
                frgb_t colorq = (frgb_t){{ (float)nextHue, 1.000f, 0.500f }};
                frgb_from_HTY(&colorq);
                color[iq] = colorq;
                nextHue += phi; 
                while (nextHue > 1.0) { nextHue -= 1.0; }
              }
            else
              { /* Non-chief family member: */
                assert((ic >= 0) && (ic < iq));
                color[iq] = color[ic];
              }
          }
      }
      
    auto void draw_indiv_pass_0(int32_t iq);
    auto void draw_indiv_pass_1(int32_t iq);
    auto void draw_indiv_pass_2(int32_t iq);
    auto void draw_indiv_pass_3(int32_t iq);
      /* Draws parts of an individual {iq} that belong to layer 0 and layer 1,
        respectively.  */
    
    auto void draw_life_trace(int32_t col0, int32_t col1, int32_t row, frgb_t *rgb);
      /* Draws the trace of an individual on the given {row},  
        spanning columns {col0..col1} (corresponding to its 
        times of birth and death). */
        
    auto void draw_arrow(int32_t col, int32_t row0, int32_t row1, bool_t shadow);
      /* Draws a vertical arrow from row {row0}
        to row {row1} on column {col}.  Meant to be 
        from some point on an individual's trace to the 
        start of the trace of one of its children.  If {shadow},
        draws the white shadow of the arrow, else writes the 
        arrow itself. */
      
    auto void draw_tail_dot(int32_t col, int32_t row);
      /* Draws a dot at the given row and column, meant to be
        the point on the trace of some individial when it had a child. */
      
    auto void draw_birth_dot(int32_t col, int32_t row, bool_t fill);
      /* Draws a dot at the given column and row, meant to indicate 
        the birth of an individual. */
   
    /* Plot the individuals: */
    for (uint32_t pass = 0;  pass < 4; pass++) 
      { /* Pass 0 draws lives, pass 1 draws arrows. */
        for (uint32_t iq = 0;  iq < ni; iq++)
          { if (rdr[iq] != -1)
              { if (pass == 0)
                  { draw_indiv_pass_0(iq); }
                else if (pass == 1)
                  { draw_indiv_pass_1(iq); }
                else if (pass == 2)
                  { draw_indiv_pass_2(iq); }
                else if (pass == 3)
                  { draw_indiv_pass_3(iq); }
              }
          }
      }

    if (color != NULL) { free(color); }

    if (debug) { fprintf(stderr, "    < %s\n", __FUNCTION__); }
    return;

    void draw_indiv_pass_0(int32_t iq)
      { /* Life traces */
        int32_t tbrq = dt[iq].tbr;
        int32_t tdtq = dt[iq].tdt;
        /* Draw fat life trace segment: */
        assert(tdtq >= tbrq);
        int32_t col0 = tbrq - tMin;
        int32_t col1 = tdtq - tMin;
        int32_t row = rdr[iq];
        if (debug) 
          { fprintf(stderr, "      life span of iq = %d times {%d .. %d}", iq, tbrq, tdtq);
            fprintf(stderr, " cols = {%d .. %d} row = %d\n", col0, col1, row);
          }
        frgb_t *rgb = (color == NULL ? &(defColor) : &(color[iq]));
        draw_life_trace(col0, col1, row, rgb);
      }

    void draw_indiv_pass_1(int32_t iq)
      { /* Shadows of arrows: */
        int32_t tbrq = dt[iq].tbr; /* Time of birth. */
        int32_t ip = dt[iq].par; /* Relevant parent index. */
        if ((tbrq >= tMin) && (tbrq <= tMax) && (ip != -1))
          { /* Birth of {q} is visible and has a relevant parent: */
            assert((ip >= 0) && (ip < ni));
            int32_t col = tbrq - tMin;
            int32_t rowq = rdr[iq];
            int32_t rowp = rdr[ip];
            draw_arrow(col, rowp, rowq, TRUE);
          }
      }

    void draw_indiv_pass_2(int32_t iq)
      { /* Arrows proper and tail dots: */
        int32_t tbrq = dt[iq].tbr; /* Time of birth. */
        int32_t ip = dt[iq].par; /* Relevant parent index. */
        if ((tbrq >= tMin) && (tbrq <= tMax) && (ip != -1))
          { /* Birth of {q} is visible and has a relevant parent: */
            assert((ip >= 0) && (ip < ni));
            int32_t col = tbrq - tMin;
            int32_t rowq = rdr[iq];
            int32_t rowp = rdr[ip];
            draw_arrow(col, rowp, rowq, FALSE);
            draw_tail_dot(col, rowp);
          }
      }

    void draw_indiv_pass_3(int32_t iq)
      { /* Birth dots: */
        int32_t tbrq = dt[iq].tbr; /* Time of birth. */
        if ((tbrq >= tMin) && (tbrq <= tMax))
          { /* Birth of {q} is visible: */
            int32_t col = tbrq - tMin;
            int32_t rowq = rdr[iq];
            bool_t fillq = (fill == NULL ? TRUE : fill[iq]);
            draw_birth_dot(col, rowq, fillq);
          }
      }

    void draw_life_trace(int32_t col0, int32_t col1, int32_t row, frgb_t *rgb)
      { if (debug) { fprintf(stderr, "      %s: columns {%d..%d}\n", __FUNCTION__, col0, col1); }
        if (col0 < 0) { col0 = 0; }
        if (col1 >= ncols) { col1 = ncols-1; }
        if (col0 > col1) { return; }
        demand((0 <= col0) && (col0 <= col1) && (col1 < ncols), "bad column range");
        demand((0 <= row) && (row < nrows), "bad row");
        
        double X0 = drtree_plot_cell_center(col0, ncols, Xstep);
        double Y0 = drtree_plot_cell_center(row, nrows, Ystep);
        
        double X1 = drtree_plot_cell_center(col1, ncols, Xstep);
        double Y1 = Y0;

        epswr_set_pen(eps, rgb->c[0], rgb->c[1], rgb->c[2], lifw, 0,0);
        epswr_segment(eps, X0, Y0, X1, Y1);
      }

    void draw_arrow(int32_t col, int32_t row0, int32_t row1, bool_t shadow)
      { if ((col < 0) || (col >= ncols)) { return; } 
        demand((0 <= row0) && (row0 < nrows), "bad row {row0}");
        demand((0 <= row1) && (row1 < nrows), "bad row {row1}");
        demand(row0 != row1, "row collision");

        double X0 = drtree_plot_cell_center(col, ncols, Xstep);
        double Y0 = drtree_plot_cell_center(row0, nrows, Ystep);
        
        double X1 = X0;
        double Y1 = drtree_plot_cell_center(row1, nrows, Ystep);

        /* Shorten arrow to end on edge of birth dot: */
        Y1 += (Y1 > Y0 ? -dbrr : + dbrr);
        
        if (shadow)
          { /* Draw the white shadow: */
            epswr_set_pen(eps, 1.000,1.000,1.000, arlw + 2*halw, 0,0);
            epswr_segment(eps, X0, Y0, X1, Y1);
            epswr_arrowhead(eps, X0, Y0, X1, Y1, awid/2, awid/2, alen, 1.0, TRUE, TRUE);
          }
        else
          { /* Draw the arrow proper: */
            epswr_set_pen(eps, 0.000,0.000,0.000, arlw, 0,0);
            epswr_segment(eps, X0, Y0, X1, Y1);
            epswr_arrowhead(eps, X0, Y0, X1, Y1, awid/2, awid/2, alen, 1.0, TRUE, TRUE);
          }
      }

    void draw_birth_dot(int32_t col, int32_t row, bool_t exp)
      { if ((col < 0) || (col >= ncols)) { return; } 
        demand((0 <= row) && (row < nrows), "bad row {row}");
        
        double X = drtree_plot_cell_center(col, ncols, Xstep);
        double Y = drtree_plot_cell_center(row, nrows, Ystep);

        epswr_set_pen(eps, 0.000,0.000,0.000, arlw, 0,0);
        epswr_dot(eps, X, Y, dbrr, exp, TRUE);
      }

    void draw_tail_dot(int32_t col, int32_t row)
      { if ((col < 0) || (col >= ncols)) { return; } 
        demand((0 <= row) && (row < nrows), "bad row {row}");
        
        double X = drtree_plot_cell_center(col, ncols, Xstep);
        double Y = drtree_plot_cell_center(row, nrows, Ystep);

        epswr_set_pen(eps, 0.000,0.000,0.000, arlw, 0,0);
        epswr_dot(eps, X, Y, datr, TRUE, TRUE);
      }

  }
  
double drtree_plot_cell_center(int32_t i, int32_t n, double step)
  { int32_t i_p = (step > 0 ? i : n - 1 - i);
    return fabs(step)*(i_p + 0.5);
  }

void drtree_plot_time_line
  ( epswr_figure_t *eps, 
    int32_t tMin,
    int32_t ncols, 
    int32_t nrows, 
    double Xstep, 
    double Ystep, 
    frgb_t *rgb,
    int32_t tRef
  )
  {
    int32_t col = tRef - tMin;
    if ((col < 0) || (col >= ncols)) { return; }

    double X = drtree_plot_cell_center(col, ncols, Xstep);
    double Y0 = drtree_plot_cell_center(0, nrows, Ystep);
    double Y1 = drtree_plot_cell_center(nrows - 1, nrows, Ystep);

    double linw = 0.75*Xstep;
    
    epswr_set_pen(eps, rgb->c[0],rgb->c[1],rgb->c[2], linw, 0,0);
    epswr_segment(eps, X, Y0, X, Y1);
  }

epswr_figure_t *drtree_plot_create_eps_figure
  ( char *name,
    int32_t ncols, 
    int32_t nrows,
    double Xstep,
    double Ystep
  )
  {
    Xstep = fabs(Xstep);
    Ystep = fabs(Ystep);
   
    double Xsize = Xstep*ncols;
    double Ysize = Ystep*nrows;
    
    /* Device plot window size and margins (pt): */
    double hPlotSize = Xsize*epswr_pt_per_mm;
    double vPlotSize = Ysize*epswr_pt_per_mm;
    double figMrg = 2.0*epswr_pt_per_mm;  /* Margin outside device plot window. */
    
    bool_t eps_verbose = FALSE;

    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, NULL, name, -1, NULL, 
        hPlotSize, vPlotSize, figMrg, figMrg, figMrg, figMrg, 
        eps_verbose
      );
    
    /* Set the client plot window: */
    double cliMrg = 2.0; /* Margin between device window and client window (mm): */
    double xMin = -cliMrg, xMax = Xsize + 2*cliMrg;
    double tMin = -cliMrg, tMax = Ysize + 2*cliMrg;
    epswr_set_client_window(eps, xMin, xMax, tMin, tMax);
    
    return eps;
  }

