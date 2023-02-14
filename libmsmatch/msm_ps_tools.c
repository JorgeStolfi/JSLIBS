/* See msm_ps_tools.h */
/* Last edited on 2023-02-14 17:45:58 by stolfi */

#define msm_ps_tools_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <stdint.h>
#include <msm_ps_tools.h>

#include <epswr.h>
#include <affirm.h>
#include <jsmath.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

/* INTERNAL REPRESENTATION */

struct msm_ps_tools_t
  { epswr_figure_t *eps;  /* The EPS figure-wrting object. */
    double size[2];       /* Usable figure dimensions (mm). */
    double fontSize;      /* Label font size (pt). */
    int32_t maxXLabChars;     /* Assumed width of X-scale labels (chars). */
    int32_t maxYLabChars;     /* Assumed width of Y-scale labels (chars). */
    /* Graph-to-epswr_figure_t coordinate conversion: */
    double epswrMin[2], epswrMax[2]; /* Min and max {epswr} client coords in each axis (mm). */
    double graphMin[2], graphMax[2];   /* Min and max Graph coords in each axis. */
  };
  /* The Postscript stream {eps} is set up so that the plot area is
    {size[0]} by {size[1]} mm, and the {epswr.h} client coordinates
    (here called /Epswr coordinates/) are measured in mm from the
    lower left corner.
    
    The {msm_ps_tools.h} Graph coordinates are such that the rectangle
    {[graphMin[0] _ graphMax[0]] × [graphMin[1] _ graphMax[1]]} in Graph coords gets
    plotted to the rectangle {[epswrMin[0]_epswrMax[0]] × [epswrMin[1]_epswrMax[1]]}
    in Epswr coords. */

void msm_ps_tools_debug_rectangle(char *title, double vMin[], double vMax[]);
  /* Prints the rectangle {[vMin[0] _ vmax[0]] x [vMin[1] x vMax[1]]} to {stderr},
    labeled with {title}. */

msm_ps_tools_t *msm_ps_tools_new
  ( FILE *wr,
    char *name,
    char *tag,
    double hSize, 
    double vSize,
    double fontSize,
    int32_t maxXLabChars,
    int32_t maxYLabChars,
    double mrg
  )
  { /* Join {name} and {tag} into a file name prefix: */
    if (wr == NULL)
      { /* Add the extension ".eps" and open the file: */
        wr = msm_open_write(name, tag, ".eps", TRUE);
      }
    /* Create an Encapsulated Postscript stream from {wr}: */
    double ptPmm = msm_PT_PER_MM;
    double hSize_pt = (hSize + 2*mrg) * ptPmm; /* (pt) */
    double vSize_pt = (vSize + 2*mrg) * ptPmm; /* (pt) */
    double mrg_pt = mrg*ptPmm;
    epswr_figure_t *eps = epswr_new_figure(wr, hSize_pt, vSize_pt, mrg_pt, mrg_pt, mrg_pt, mrg_pt, TRUE);
    /* Set {epswr} client coordinates to mm: */
    epswr_set_client_window(eps, 0, hSize, 0, vSize); 
    /* Create the {msm_ps_tools_t} object: */
    msm_ps_tools_t *mps = (msm_ps_tools_t *)notnull(malloc(sizeof(msm_ps_tools_t)), "no mem");
    mps->eps = eps;
    mps->size[0] = hSize; mps->size[1] = vSize;
    /* Set default Graph and Epswr reference window: */
    int32_t ax;
    for (ax = 0; ax < 2; ax++)
      { /* Epswr ref window is the whole plot area (canvas minus margins): */
        mps->epswrMin[ax] = 0; mps->epswrMax[ax] = mps->size[ax];
        /* Graph ref window is initially set so that cient and Epswr coords are the same: */
        mps->graphMin[ax] = mps->epswrMin[ax]; 
        mps->graphMax[ax] = mps->epswrMax[ax];
      }
    mps->fontSize = fontSize;
    mps->maxXLabChars = maxXLabChars;
    mps->maxYLabChars = maxYLabChars;
    return mps;
  }

msm_ps_tools_t *msm_ps_tools_new_graph
  ( FILE *wr,
    char *name,
    char *tag,
    double hGraphSize, 
    double vGraphSize,
    bool_t scaleL, bool_t titleL,
    bool_t scaleR, bool_t titleR,
    bool_t scaleB, bool_t titleB,
    bool_t scaleT, bool_t titleT,
    double fontSize,
    int32_t maxXLabChars,
    int32_t maxYLabChars,    
    double mrg
  )
  {
    double ptPmm = msm_PT_PER_MM;
    
    /* Estimate the character dimensions in mm: */
    double vCharSize = fontSize/ptPmm; /* Estimated character height (mm), with skosh. */
    double hCharSize = 0.8*vCharSize;  /* Estimated character width (mm), with skosh; guess. */
    
    /* Compute the max label dimensions in mm: */
    double hLabelSize = maxXLabChars*hCharSize;
    double vLabelSize = vCharSize;
    
    /* Initialize title/scale margin widths in mm: */
    double mrgL = 0.0; /* Left. */
    double mrgR = 0.0; /* Right. */
    double mrgB = 0.0; /* Bottom. */
    double mrgT = 0.0; /* Top. */
    
    /* Add space for X and Y titles. Assume that Y titles are rotated. */
    if (titleL) { mrgL += vCharSize; }
    if (titleR) { mrgR += vCharSize; }
    if (titleB) { mrgB += vCharSize; }
    if (titleT) { mrgT += vCharSize; }
    
    /* Add space for X and Y scales. Assume that Y scales are NOT rotated. */
    if (scaleL) { mrgL += hLabelSize; }
    if (scaleR) { mrgR += hLabelSize; }
    if (scaleB) { mrgB += vLabelSize; }
    if (scaleT) { mrgT += vLabelSize; }
    
    /* Make sure that there is enough space for overshoot of scale labels: */
    double hOver = 0.50*hLabelSize; /* Estimated H overshoot of X scale labels */
    double vOver = 0.75*vLabelSize; /* Estimated V overshoot of Y scale labels. */
    if (scaleL || scaleR)
      { if (mrgB < vOver) { mrgB = vOver; }
        if (mrgT < vOver) { mrgT = vOver; }
      }
    if (scaleB || scaleT)
      { if (mrgL < hOver) { mrgL = hOver; }
        if (mrgR < hOver) { mrgR = hOver; }
      }
    
    /* Compute dimensions of plottable area (graph,titles,scales): */
    double hTotSize = hGraphSize + mrgL + mrgR;
    double vTotSize = vGraphSize + mrgB + mrgT;

    /* Create plot object: */
    msm_ps_tools_t *mps = msm_ps_tools_new(wr, name, tag, hTotSize, vTotSize, fontSize, maxXLabChars, maxYLabChars, mrg);
    
    /* Set reference rectangles to graph area only: */
    double hMin = 0 + mrgL, hMax = hTotSize - mrgR;
    double vMin = 0 + mrgB, vMax = vTotSize - mrgT;
    msm_ps_tools_set_epswr_ref_window(mps, hMin, hMax, vMin, vMax);
    msm_ps_tools_set_graph_ref_window(mps, hMin, hMax, vMin, vMax);
    
    return mps;
  }

void msm_ps_tools_get_plot_size(msm_ps_tools_t *mps, double *hSize, double *vSize)
  { *hSize = mps->size[0]; *vSize = mps->size[1]; }
   
void msm_ps_tools_debug_rectangle(char *title, double vMin[], double vMax[])
  { fprintf
      ( stderr, 
        "%s ref window [%10.6f %10.6f]  [%10.6f %10.6f]\n", 
        title, vMin[0], vMax[0], vMin[1], vMax[1]
      );
  }

void msm_ps_tools_set_epswr_ref_window
  ( msm_ps_tools_t *mps,
    double hMin, double hMax, 
    double vMin, double vMax
  )
  { mps->epswrMin[0] = hMin; mps->epswrMax[0] = hMax;
    mps->epswrMin[1] = vMin; mps->epswrMax[1] = vMax;
    msm_ps_tools_debug_rectangle("Epswr", mps->epswrMin, mps->epswrMax);
  }

void msm_ps_tools_shrink_epswr_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  )
  { mps->epswrMin[0] += lMrg; mps-> epswrMax[0] -= rMrg;
    mps->epswrMin[1] += bMrg; mps-> epswrMax[1] -= tMrg;
    msm_ps_tools_debug_rectangle("Epswr", mps->epswrMin, mps->epswrMax);
  }

void msm_ps_tools_set_graph_ref_window
  ( msm_ps_tools_t *mps,
    double xMin, double xMax, 
    double yMin, double yMax
  )
  { mps->graphMin[0] = xMin; mps->graphMax[0] = xMax;
    mps->graphMin[1] = yMin; mps->graphMax[1] = yMax;
    msm_ps_tools_debug_rectangle("Graph", mps->graphMin, mps->graphMax);
  }

void msm_ps_tools_expand_graph_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  )
  { mps->graphMin[0] -= lMrg; mps-> graphMax[0] += rMrg;
    mps->graphMin[1] -= bMrg; mps-> graphMax[1] += tMrg;
    msm_ps_tools_debug_rectangle("Graph", mps->graphMin, mps->graphMax);
  }

double msm_ps_tools_map_x(msm_ps_tools_t *mps, double x)
  { double s = (x - mps->graphMin[0])/(mps->graphMax[0] - mps->graphMin[0]), r = 1-s;
    return r*mps->epswrMin[0] + s*mps->epswrMax[0];
  }
  
double msm_ps_tools_map_y(msm_ps_tools_t *mps, double y)
  { double s = (y - mps->graphMin[1])/(mps->graphMax[1] - mps->graphMin[1]), r = 1-s;
    return r*mps->epswrMin[1] + s*mps->epswrMax[1];
  }

double msm_ps_tools_map_coord(msm_ps_tools_t *mps, epswr_axis_t axis, double coord)
  { demand((axis == epswr_axis_HOR) || (axis == epswr_axis_VER), "bad axis");
    double s = (coord - mps->graphMin[axis])/(mps->graphMax[axis] - mps->graphMin[axis]), r = 1-s;
    return r*mps->epswrMin[axis] + s*mps->epswrMax[axis];
  }
  
void msm_ps_tools_map_coords(msm_ps_tools_t *mps, double x, double y, double *h, double *v)
  { (*h) = msm_ps_tools_map_x(mps, x);
    (*v) = msm_ps_tools_map_y(mps, y);
  }

/* EPSWR -> GRAPH */

double msm_ps_tools_unmap_h(msm_ps_tools_t *mps, double h)
  { double s = (h - mps->epswrMin[0])/(mps->epswrMax[0] - mps->epswrMin[0]), r = 1-s;
    return r*mps->graphMin[0] + s*mps->graphMax[0];
  }
  
double msm_ps_tools_unmap_v(msm_ps_tools_t *mps, double v)
  { double s = (v - mps->epswrMin[1])/(mps->epswrMax[1] - mps->epswrMin[1]), r = 1-s;
    return r*mps->graphMin[1] + s*mps->graphMax[1];
  }

double msm_ps_tools_unmap_coord(msm_ps_tools_t *mps, epswr_axis_t axis, double coord)
  { demand((axis == epswr_axis_HOR) || (axis == epswr_axis_VER), "bad axis");
    double s = (coord - mps->epswrMin[axis])/(mps->epswrMax[axis] - mps->epswrMin[axis]), r = 1-s;
    return r*mps->graphMin[axis] + s*mps->graphMax[axis];
  }
  
void msm_ps_tools_unmap_coords(msm_ps_tools_t *mps, double h, double v, double *x, double *y)
  { (*x) = msm_ps_tools_unmap_h(mps, h);
    (*y) = msm_ps_tools_unmap_v(mps, v);
  }

/* AXES, TICS, ETC */

void msm_ps_tools_draw_segment(msm_ps_tools_t *mps, double xa, double ya, double xb, double yb)
  { /* Get the Epswr coordinate of the points: */
    double ha, va;
    msm_ps_tools_map_coords(mps, xa, ya, &ha, &va);
    double hb, vb;
    msm_ps_tools_map_coords(mps, xb, yb, &hb, &vb);
    /* Plot the segment: */
    epswr_segment(mps->eps, ha, va, hb, vb);
  }

void msm_ps_tools_draw_ref_axis(msm_ps_tools_t *mps, epswr_axis_t axis, double R, double G, double B)
  { epswr_set_pen(mps->eps, R,G,B, 0.10, 0.0, 0.0);
    /* Get the Epswr coordinate of the requested axis along the other axis: */
    double cpos = msm_ps_tools_map_coord(mps, 1-axis, 0.0);
    epswr_axis(mps->eps, axis, cpos, mps->epswrMin[axis], mps->epswrMax[axis]);
  }

void msm_ps_tools_draw_ref_frame(msm_ps_tools_t *mps, double R, double G, double B)
  { epswr_set_pen(mps->eps, R,G,B, 0.10, 0.0, 0.0);
    epswr_rectangle(mps->eps, mps->epswrMin[0], mps->epswrMax[0], mps->epswrMin[1], mps->epswrMax[1], FALSE, TRUE);
  }

void msm_ps_tools_draw_plot_frame(msm_ps_tools_t *mps, double R, double G, double B)
  { epswr_set_pen(mps->eps, R,G,B, 0.10, 0.0, 0.0);
    epswr_rectangle(mps->eps, 0, mps->size[0], 0, mps->size[1], FALSE, TRUE);
  }
  
void msm_ps_tools_draw_tic
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double xt, 
    double yt,
    double ticSize,
    double ticAlign,
    double R, double G, double B,
    char *label,
    double labAlign
  )
  { /* Compute nominal Epswr coordinates {ht,vt} of tic: */
    double ht = msm_ps_tools_map_x(mps, xt);
    double vt = msm_ps_tools_map_y(mps, yt);
    /* Plot the tic proper: */
    epswr_set_pen(mps->eps, R,G,B, 0.10, 0.0, 0.0);
    epswr_tic(mps->eps, axis, ht, vt, ticSize, ticAlign);
    if (label != NULL)
      { /* Decide the vertical and horizontal alignment of label: */
        double hAlign = (axis == epswr_axis_HOR ? 0.5 : labAlign);
        double vAlign = (axis == epswr_axis_VER ? 0.5 : labAlign);
        /* Adjust the position of the label's reference point to account for tic size: */
        double cadj = (1-labAlign-ticAlign)*ticSize; 
        double hb = (axis == epswr_axis_HOR ? ht : ht + cadj);
        double vb = (axis == epswr_axis_VER ? vt : vt + cadj);
        /* Set pen to requested color: */
        epswr_set_pen(mps->eps, 0,0,0, 0.10, 0.0, 0.0);
        epswr_set_fill_color(mps->eps, 0,0,0);
        epswr_label(mps->eps, label, "Rg", hb, vb, 0.0, FALSE, hAlign, vAlign, TRUE, FALSE);
      }
  }

void msm_ps_tools_draw_scale
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double pos,
    double ticSize,
    double ticAlign,
    double ticMinDist, 
    double ticMinStep, 
    double R, double G, double B,
    char *fmt,
    double labAlign,
    double labMinDist,
    double labMinStep
  )
  { 
    /* Oh paranoia, not even Goya could draw'ya: */
    demand((axis == epswr_axis_HOR) || (axis == epswr_axis_VER), "bad axis");
    
    /* Get the raw span of tics in Epswr coords (mm): */
    double epswrMin = mps->epswrMin[axis];
    double epswrMax = mps->epswrMax[axis];
    
    /* Swap if needed to ensure {epswrMin <= epswrMax}: */
    if (epswrMin > epswrMax) { double c = epswrMin; epswrMin = epswrMax; epswrMax = c; }
    
    /* Leave some space {eps} between range limits and first/last tics: */
    double eps = 0.5; /* In mm. */
    epswrMin += eps; epswrMax -= eps;
    
    /* Choose the Graph coords of minor tics (0 if no tics in range): */
    double ztMin, ztMax; /* Min and max Graph coords of minor tics. */
    double ztStep; /* Graph coord increment between minor tics (0 if no tics). */
    msm_ps_tools_choose_tic_coords
      ( mps, axis, epswrMin, epswrMax, ticMinDist, ticMinStep, &ztMin, &ztMax, &ztStep );
    /* fprintf(stderr, "  ztMin = %24.16e ztMax = %24.16e ztStep = %24.16e\n", ztMin, ztMax, ztStep); */
    if (ztStep <= 0) { /* Not enough space for tics: */ return; }
    
    /* Compute number {ntsteps} of minor tic steps: */
    int32_t ntsteps = (int32_t)msm_round((ztMax/ztStep) - (ztMin/ztStep));
    assert(ntsteps >= 0); 
    
    /* Compute minor-to-major tic ratio {labPer} (0 if no major tics) and phase {labSkp}: */
    int32_t labPer;   /* Minor tic steps in one major tic step. */
    int32_t labSkp; /* Minor tics before the first major tic. */
    if (fmt == NULL)
      { /* Graphs wants no major tics: */
        labPer = labSkp = 0;
      }
    else 
      { /* Choose major tics honoring {labSp} and {minStep}: */
        msm_choose_label_coords
          ( mps, axis, ztMin, ztMax, ztStep, labMinDist, labMinStep, &labPer, &labSkp);
      }
      
    int32_t i;
    for (i = 0; i <= ntsteps; i++)
      { /* Get fractions {r:s} of {i} betwen 0 and {ntsteps}: */
        double s = (ntsteps == 0 ? 0.5 : ((double)i)/((double)ntsteps)), r = 1-s;
        /* Compute nominal Graph coordinate {zt} of tic on the given {axis}: */
        double zt = r*ztMin + s*ztMax;
        /* Nominal Graph coords of tic mark are {xt,yt}: */
        double xt = (axis == epswr_axis_HOR ? zt : pos);
        double yt = (axis == epswr_axis_VER ? zt : pos);
        /* Decide which kind of label to draw: */
        if ((labPer > 0) && (imod(i, labPer) == labSkp))
          { /* Draw a major (labeled) tic: */
            assert(fmt != NULL);
            /* Format the text: */
            char *label = NULL;
            asprintf(&label, fmt, zt);
            /* Draw a labeled tic: */
            msm_ps_tools_draw_tic(mps, axis, xt, yt, 1.5*ticSize, ticAlign, R,G,B, label, labAlign);
            free(label);
          }
        else
          { /* Draw a minor (label-less) tic: */
            msm_ps_tools_draw_tic(mps, axis, xt, yt, ticSize, ticAlign, R,G,B, NULL, 0.0);
          }
      }
  }

void msm_ps_tools_choose_tic_coords
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double epswrMin, 
    double epswrMax, 
    double minDist, 
    double minStep, 
    double *zMinP,
    double *zMaxP,
    double *zStepP
  )
  { /* Check required conditions: */
    demand((minDist > 0) || (minStep > 0), "neither minDist nor MinStep were given"); 
    /* Local copies of results: */
    double graphMin, graphMax, zStep;
    /* Any space at all? */
    if (epswrMin > epswrMax) 
      { /* Interval is empty, no tics there: */
        graphMin = +INF; graphMax = -INF; zStep = 0;
      }
    else
      { /* Compute the Graph increment {zStep} equivalent to Epswr increment {labSp}: */
        double zDif = fabs(mps->graphMax[axis] - mps->graphMin[axis]);
        double cDif = fabs(mps->epswrMax[axis] - mps->epswrMin[axis]);
        zStep = minDist*(zDif/cDif);
        if (zStep < minStep) { zStep = minStep; }
        assert(zStep >= 0);
        /* Round {zStep} to a nice value: */
        zStep = epswr_round_to_nice(zStep);
        /* Get range {[graphMin_graphMax]} of Graph coords corresponding to {[epswrMin_epswrMax]}: */
        graphMin = msm_ps_tools_unmap_coord(mps, axis, epswrMin);
        graphMax = msm_ps_tools_unmap_coord(mps, axis, epswrMax);
        if (graphMin > graphMax) { double z = graphMin; graphMin = graphMax; graphMax = z; }
        /* fprintf(stderr, "  graphMin =  %24.16e graphMax =  %24.16e  zStep = %24.16e\n", graphMin, graphMax, zStep); */
        /* See whether there are any multiples of {zStep} in {graphMin,graphMax}: */
        assert(graphMin <= graphMax);
        double qMin = ceil(graphMin/zStep);
        double qMax = floor(graphMax/zStep);
        /* fprintf(stderr, "  qMin =  %24.16e qMax =  %24.16e\n", qMin, qMax); */
        if (qMin <= qMax)
          { /* Synchronize {graphMin,graphMax} to multiples of {zStep}: */
            double zrMin = qMin*zStep; assert(zrMin >= graphMin); 
            double zrMax = qMax*zStep; assert(zrMax <= graphMax);
            /* fprintf(stderr, "  zrMin = %24.16e zrMax = %24.16e\n", zrMin, zrMax); */
            assert(zrMin <= zrMax);
            graphMin = zrMin; graphMax = zrMax;
          }
        else
          { /* No multiples of {zStep} in range: */
            graphMin = +INF; graphMax = -INF; zStep = 0;
          }
      }
    /* fprintf(stderr, "  graphMin =  %24.16e graphMax =  %24.16e  zStep = %24.16e\n", graphMin, graphMax, zStep); */
    (*zMinP) = graphMin; (*zMaxP) = graphMax; (*zStepP) = zStep;
  }

void msm_choose_label_coords
  ( msm_ps_tools_t *mps,
    epswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int32_t *labPerP,
    int32_t *labSkpP
  )
  {
    double eps = 1.0e-12; /* Relative fudge factor to compensate for roundoff errors. */

     /* First, we choose the label period {labPer}: */
    int32_t labPer; 
    
    auto int32_t min_mult(int32_t k, double unit, double min);
      /* Increments the positive integer {k} so that {k*unit >= minv},
        allowing for some roundoff noise in {unit}. */
         
    int32_t min_mult(int32_t k, double unit, double minv)
      { unit = (1+eps)*unit;
        if (k*unit < minv)
          { k = (int32_t)ceil(minv/unit);
            assert(k*unit >= minv);
          }
        return k;
      }
      
    if ((ztStep == 0) || (ztMin > ztMax))
      { /* No tics, no major tics: */
        labPer = 0;
      }
    else
      { /* Set {labPer} to the minmum period that satisfies {labSp}: */
        double ticDist; /* Actual tic spacing {ticDist} in Epswr coords (fudged up): */
        ticDist = fabs(msm_ps_tools_map_coord(mps, axis, ztStep) - msm_ps_tools_map_coord(mps, axis, 0));
        labPer = min_mult(1, ticDist, minDist);
        assert(labPer > 0);

        /* Adjust {labPer} upwards to honor {minStep}: */
        demand(ztStep > 0, "invalid ztStep"); 
        labPer = min_mult(labPer, ztStep, minStep);
        assert(labPer > 0);

        /* Adjust {labPer} so that the label value increment is nice: */
        if (labPer > 1)
          { /* Compute the label step {zbStep} in Graph coords: */
            double zbStep = labPer*ztStep;
            /* Increment {labPer} until {zbStep} is nice (should't be too far): */
            while (TRUE)
              { if (zbStep <= minStep) { break; }
                double zbNice = epswr_round_to_nice(zbStep);
                /* fprintf(stderr, "    labPer = %d zbStep =  %24.16e", labPer, zbStep); */
                /* fprintf(stderr, " zbNice = %24.16e\n", zbNice); */
                if (zbNice == INF) { /* Overflow, give up: */ labPer = 0; break; }
                if (zbNice <= zbStep) { /* Phew! */ break; }
                int32_t labPerNew = min_mult(labPer, ztStep, zbNice);
                if(labPerNew <= labPer) break;
                labPer = labPerNew;
                zbStep = labPer*ztStep;
              }
          }
      }
      
    /* Compute the number {labSkp} of minor tics before first major tic: */
    int32_t labSkp;
    if (labPer == 0)
      { /* Just in case: */ labSkp = 0; }
    else
      { /* Compute coords {zbMin,zbMax} of first and last major tics: */
        double zbStep = labPer*ztStep; /* Label spacing in Graph coords. */
        double qbMin = ceil(ztMin/zbStep);
        double qbMax = floor(ztMax/zbStep);
        /* fprintf(stderr, "  qbMin =  %24.16e qbMax =  %24.16e\n", qbMin, qbMax); */
        if (qbMin <= qbMax)
          { /* There is a multiple of {zbStep} in the range: */
            double zbMin = qbMin*zbStep;
            labSkp = (int32_t)msm_round((zbMin - ztMin)/ztStep);
            /* fprintf(stderr, "  zbMin  =  %24.16e labSkp = %d\n", zbMin, labSkp); */
            assert((labSkp >= 0) && (labSkp < labPer));
          }
        else
          { /* No multiples of {zbStep} in range: */
            labPer = labSkp = 0;
          }
      }
    
    /* Return results: */
    (*labPerP) = labPer; (*labSkpP) = labSkp;
  }

void msm_ps_tools_draw_y_polyline
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int32_t n
  )
  { if (n < 2) { return; }
    double n1 = n - 1;
    double h0, v0; /* Previously plotted point. */
    msm_ps_tools_map_coords(mps, xMin, y[0], &h0, &v0);
    int32_t i;
    for (i = 1; i < n; i++)
      { /* Plot step from point {i-1} to point {i}: */
        double s = i/n1, r = 1-s;
        double x1 = r*xMin + s*xMax;
        double h1, v1;
        msm_ps_tools_map_coords(mps, x1, y[i], &h1, &v1);
        epswr_segment(mps->eps, h0, v0, h1, v1);
        h0 = h1; v0 = v1;
      }
  }

void msm_ps_tools_draw_y_dots
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int32_t n,
    double rad,
    bool_t fill,
    bool_t draw
  )
  { if (n < 1) { return; }
    double n1 = n - 1;
    int32_t i;
    for (i = 0; i < n; i++)
      { /* Plot point {i}: */
        double s = i/n1, r = 1-s;
        double x = r*xMin + s*xMax;
        double h, v;
        msm_ps_tools_map_coords(mps, x, y[i], &h, &v);
        epswr_dot(mps->eps, h, v, rad, fill, draw);
      }
  }

void msm_ps_tools_close(msm_ps_tools_t *mps)
  { epswr_end_figure(mps->eps);
    free(mps);
  }

void msm_ps_tools_compute_data_range
  ( int32_t n, 
    int32_t stride, 
    double z[], 
    double *zMinP, 
    double *zMaxP
  )
  { int32_t i;
    double graphMin = +INF; 
    double graphMax = -INF;
    for (i = 0; i < n; i++)
      { double zi = z[i*stride];
        if (zi < graphMin) { graphMin = zi; }
        if (zi > graphMax) { graphMax = zi; }
      }
    /* Make sure that the range is not empty or trivial: */
    if (graphMin == graphMax)
      { double zAbs = fabs(graphMin);
        if (zAbs == 0) { zAbs = 1.0; }
        graphMin -= 0.1*zAbs;
        graphMax += 0.1*zAbs;
      }
    (*zMinP) = graphMin;
    (*zMaxP) = graphMax;
  }

void msm_ps_tools_draw_graphs
  ( msm_ps_tools_t *mps,
    int32_t nc,
    int32_t nd,
    double x[],
    double start,
    double step,
    double y[],
    double yMin,
    double yMax
  )
  {
    demand(nd > 0, "can't plot empty graph");

    /* Get the {epswr_figure_t} handle for {epswr.h} routines: */
    epswr_figure_t *eps = mps->eps;

    /* Choose the actual X range {[xPlotMin_xPlotMax]} of graphs: */
    /* If {circ} is true, allow for half-steps at each end: */
    double xPlotMin, xPlotMax;
    if (x == NULL)
      { xPlotMin = start; xPlotMax = start + (nd - 1)*step; }
    else
      { demand(x[0] <= x[nd-1], "abscissas out of order");
        xPlotMin = x[0];
        xPlotMax = x[nd-1];
      }
    /* Choose the nominal plot X range {[xMin_xMax]}: */
    double xMag = fmax(fabs(xPlotMin), fabs(xPlotMax)) + 1.0;
    double xDif = fabs(xPlotMax - xPlotMin) + 1.0;
    double xSkosh = fmax(1.0e-12*xMag, 0.05*xDif); /* Extra skosh. */
    double xMin = xPlotMin - xSkosh;
    double xMax = xPlotMax + xSkosh;
    
    if (yMin >= yMax)
      { /* Given range {[yMin _ yMax]} is empty or trivial, compute it from the data: */
        int32_t ncd = nc*nd;
        msm_ps_tools_compute_data_range(ncd, 1, y, &yMin, &yMax); 
        /* Include 0 in the range, if near enough: */
        double yDifRaw = yMax - yMin;
        if ((yMin > 0) && (yMin < +0.25*yDifRaw)) { yMin = 0; }
        if ((yMax < 0) && (yMax > -0.25*yDifRaw)) { yMax = 0; }
      }

    /* Add some skosh to the Y range: */
    { double yMag = fmax(fabs(yMin), fabs(yMax)) + 1.0e-6;
      double yDif = fabs(yMax - yMin) + 1.0e-6;
      double ySkosh = fmax(1.0e-12*yMag, 0.05*yDif); /* Extra skosh. */
      yMin = yMin - ySkosh;
      yMax = yMax + ySkosh;
    }
    
    /* Set Graph reference window, with some skosh: */
    msm_ps_tools_set_graph_ref_window(mps, xMin, xMax, yMin, yMax);

    /* Estimate the character dimensions in mm: */
    double ptPmm = msm_PT_PER_MM;
    double vCharSize = mps->fontSize/ptPmm; /* Estimated character height (mm), with skosh. */
    double hCharSize = 0.8*vCharSize;  /* Estimated character width (mm), with skosh; guess. */
  
    /* Shrink Epswr window to leave space for tics and labels: */
    double ticSz = 1.0;  /* In mm. */
    double yScaleWd = mps->maxYLabChars*hCharSize + (ticSz + 2.0); /* Width of Y scale and tics (mm). */
    double xScaleHt = vCharSize + (ticSz + 2.0); /* Height of X scale and tics (mm). */
    double hLabelSize = mps->maxXLabChars*hCharSize; /* Nominal width of X scale label. */
    double hOver = 0.50*hLabelSize; /* Estimated H overshoot of X scale labels */
    msm_ps_tools_shrink_epswr_ref_window(mps, yScaleWd, hOver, xScaleHt, 0.0);
    
    /* Decide whether to plot dots at individual samples: */
    double xStep; /* Average step size in Graph coords. */
    if (x == NULL)
      { xStep = step; }
    else if (nd >= 2)
      { xStep = (x[nd-1] - x[0])/(nd-1); }
    else
      { xStep = 0; }
    double hStep = msm_ps_tools_map_x(mps, xPlotMin + xStep) - msm_ps_tools_map_x(mps, xPlotMin);
    bool_t show_dots = (hStep >= 1.5);  /* If spaced at least 1.5 mm */

    /* Colors for channels (all with brightness = 0.30): */
    demand(nc <= 3, "too many channels");
    double R[3] = { 1.00, 0.00, 0.00 };
    double G[3] = { 0.00, 0.50, 0.33 };
    double B[3] = { 0.00, 0.00, 1.00 };
    
    /* Draw axes, tics, labels, etc: */
    char *font = "Times-Roman";
    epswr_set_label_font(eps, font, mps->fontSize);
    double D = 0.5; /* Brightness of axis and frame color. */
    msm_ps_tools_draw_ref_axis(mps, epswr_axis_HOR, D,D,D);
    msm_ps_tools_draw_scale(mps, epswr_axis_HOR, yMin, ticSz,1.0, 5.0,1.0, D,D,D, "%.0f",   1.2, 12.0,1.0);
    msm_ps_tools_draw_scale(mps, epswr_axis_VER, xMin, ticSz,1.0, 3.0,0.0, D,D,D, "%+.2f",  1.2,  4.0,0.0);

    /* Decide indices of first and last segment to plot. */
    /* Segment {i} extends from point {i-1} to point {i}. */
    /* If circular, plots segments {0..nd}, clipped to {[xPlotMin _ xPlotMax]}. */
    /* If not circular, plots segments {1..nd-1}, no need to clip. */
    int32_t iIni = (1); /* Index of first segment to plot. */
    int32_t iFin = (nd-1); /* Index of last segment to plot. */

    int32_t c;
    for (c = 0; c < nc; c++)
      { /* Plot channel {c}. */
        /* Set the pen color: */
        epswr_set_pen(eps, R[c], G[c], B[c], 0.25, 0.0, 0.0);
        /* Plot channel {c} of the sequence. */
        
        auto void get_data_point(int32_t ix, double *xi, double *yi);
          /* Obtains the coordinates {*xi,*yi} of data point number {ix},
            taking circularity into account. */
          
        void get_data_point(int32_t ix, double *xi, double *yi)
          { double dx = 0.0;
            assert((ix >= 0) && (ix < nd));
            (*xi) = (x == NULL ? (double)(start + ix*step) : x[ix]) + dx;
            (*yi) = y[c*nd + ix];
          }
        
        auto void get_plot_point(int32_t ix, double *xp, double *yp);
          /* Obtains the coordinates {*xp,*yp} of plot point number
            {ix}, taking circularity into account. For {ix == -1}
            returns {xp==xPlotMin}, for {ix == nd} returns
            {xp==xPlotMax}; in both cases {yp} is obtained by
            interpolation. In the other cases, returns the data point
            number {ix}. */
          
        void get_plot_point(int32_t ip, double *xp, double *yp)
          { if ((ip >= 0) && (ip < nd))
              { get_data_point(ip, xp, yp); }
          }
        
        /* Compute first point {h0,v0} of polyline: */
        double x0, y0; /* Previously plotted point (Graph coords). */
        get_plot_point(iIni-1, &x0, &y0);
        double h0, v0; /* Previously plotted point (plot coords). */
        msm_ps_tools_map_coords(mps, x0, y0, &h0, &v0);
        /* Now plot all steps. */
        int32_t i = iIni;
        while (i <= iFin)
          { /* Plot step from sample {i-1} to sample {i}: */
            double x1, y1;
            get_plot_point(i, &x1, &y1);
            double h1, v1;
            msm_ps_tools_map_coords(mps, x1, y1, &h1, &v1);
            epswr_segment(eps, h0, v0, h1, v1);
            h0 = h1; v0 = v1;
            i++;
          }
        if (show_dots)
          { /* Draw dots at data points: */
            epswr_set_fill_color(eps, R[c], G[c], B[c]);
            for (i = 0; i < nd; i++)
              { double xi, yi;
                get_data_point(i, &xi, &yi);
                double hi, vi;
                msm_ps_tools_map_coords(mps, xi, yi, &hi, &vi);
                epswr_dot(eps, hi, vi, 0.5, TRUE, FALSE);
              }
          }
      }
   
    /* Draw a thin frame in black: */
    msm_ps_tools_draw_ref_frame(mps, D,D,D);
  }

void msm_ps_tools_draw_histogram
  ( msm_ps_tools_t *mps,
    int32_t nd,
    double x[],
    double y[],
    double yMin,
    double yMax
  )
  {
    demand(nd > 0, "can't plot empty histogram");

    /* Get the {epswr_figure_t} handle for {epswr.h} routines: */
    epswr_figure_t *eps = mps->eps;

    /* Choose the actual X range {[xPlotMin_xPlotMax]} of histogram: */
    /* Allow for half-steps at each end: */
    double xPlotMin, xPlotMax;
    if (x == NULL)
      { xPlotMin = (-0.5) - 0.5; xPlotMax = (nd - 1 + 0.5) + 0.5; }
    else
      { demand(x[0] < x[nd], "invalid X range");
        double xStep = (x[nd] - x[0])/nd; /* Average step width. */
        xPlotMin = x[0] - 0.5*xStep;    
        xPlotMax = x[nd] + 0.5*xStep;
      }

    /* Choose the nominal X range {[xMin_xMax]} of plot: */
    double xMag = fmax(fabs(xPlotMin), fabs(xPlotMax));
    double xDif = fabs(xPlotMax - xPlotMin);
    double xSkosh = fmax(1.0e-12*xMag, 0.02*xDif); /* Extra skosh. */
    double xMin = xPlotMin - xSkosh;
    double xMax = xPlotMax + xSkosh;
    
    /* Set Graph reference window, with some skosh: */
    msm_ps_tools_set_graph_ref_window(mps, xMin, xMax, yMin, yMax);
  
    /* Shrink Epswr window to leave space for tics and labels: */
    double ptPmm = msm_PT_PER_MM;
    double ticSz = 1.0;  /* In mm. */
    double yScaleWd = 2*mps->fontSize/ptPmm + (ticSz + 2.0); /* Width of Y scale and tics. */
    double xScaleHt = mps->fontSize/ptPmm + (ticSz + 2.0); /* Height of X scale and tics. */
    msm_ps_tools_shrink_epswr_ref_window(mps, yScaleWd, 0.0, xScaleHt, 0.0);
    
    /* Draw axes, tics, labels, etc: */
    char *font = "Times-Roman";
    epswr_set_label_font(eps, font, mps->fontSize);
    double D = 0.5; /* Lightness of axis and frame */
    msm_ps_tools_draw_ref_axis(mps, epswr_axis_HOR, D,D,D);
    msm_ps_tools_draw_scale(mps, epswr_axis_HOR, yMin, ticSz,1.0, 5.0,1.0, D,D,D, "%.0f",   1.2, 12.0,1.0);
    msm_ps_tools_draw_scale(mps, epswr_axis_VER, xMin, ticSz,1.0, 3.0,0.0, D,D,D, "%+.2f",  1.2,  4.0,0.0);

    /* Set the pen color: */
    epswr_set_pen(eps, 1.00,0.00,0.00, 0.25, 0.0, 0.0);

    /* Prepare for first half-step: */
    double x0 = xPlotMin, y0 = 0.0;
    double h0, v0; /* End of previous step. */
    msm_ps_tools_map_coords(mps, x0, y0, &h0, &v0);

    /* Now plot all steps, including half-steps at each end ({i==-1} and {i==nd}): */
    int32_t i = -1;
    while (i <= nd)
      { /* Get upper right corner {x1,y1} of next bar: */
        double x1 = (i >= nd ? xPlotMax : (x == NULL ? i+0.5 : x[i+1]));
        double y1 = (i < 0 ? 0.0: (i >= nd ? 0.0 : y[i]));
        /* Map to plot coordinates {h1,v1}: */
        double h1, v1;
        msm_ps_tools_map_coords(mps, x1, y1, &h1, &v1);
        
        if (i >= 0)
          { /* Plot left wall of bar {i}: */
            epswr_segment(eps, h0, v0, h0, v1);
          }
        /* Plot horizontal part of bar {i}: */
        epswr_segment(eps, h0, v1, h1, v1);
        h0 = h1; v0 = v1;
        i++;
      }
   
    /* Draw a thin frame in black: */
    msm_ps_tools_draw_ref_frame(mps, D,D,D);
  }

epswr_figure_t *msm_ps_tools_get_eps_figure(msm_ps_tools_t *mps)
  { return mps->eps; }
