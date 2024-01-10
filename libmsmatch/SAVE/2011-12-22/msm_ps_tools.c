/* See msm_ps_tools.h */
/* Last edited on 2009-01-06 03:10:43 by stolfi */

#define msm_ps_tools_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_ps_tools.h>

#include <pswr.h>
#include <affirm.h>
#include <jsmath.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

/* INTERNAL REPRESENTATION */

struct msm_ps_tools_t
  { PSStream *ps;
    /* Usable figure dimensions in mm */
    double size[2]; 
    /* Label font size: */
    double fontSize;
    /* Client-to-PSStream coordinate conversion: */
    double cMin[2], zMin[2];
    double cMax[2], zMax[2];
  };
  /* The Postscript stream {ps} is set up so that the plot area is
    {size[0]} by {size[1]} mm, and the {pswr.h} client coordinates
    (here called /device coordinates/) are measured in mm from the
    lower left corner.
    
    The {msm_ps_tools.h} client coordinates are such that the rectangle
    {[zMin[0] _ zMax[0]] × [zMin[1] _ zMax[1]]} in client coords gets
    plotted to the rectangle {[cMin[0]_cMax[0]] × [cMin[1]_cMax[1]]}
    in device coords. */

msm_ps_tools_t *msm_ps_tools_new
  ( FILE *wr,
    char *name,
    char *tag,
    double hSize, 
    double vSize,
    double fontSize,
    double mrg
  )
  { /* Join {name} and {tag} into a file name prefix: */
    char *prefix = NULL;
    asprintf(&prefix, "%s%s", name, tag);
    if (wr == NULL)
      { /* Add the extension ".eps" and open the file: */
        wr = msm_open_write(name, tag, ".eps", TRUE);
      }
    /* Create an Encapsulated Postscript stream from {wr}: */
    double mm = msm_PT_PER_MM;
    double hPageSize = hSize + 2*mrg;
    double vPageSize = vSize + 2*mrg;
    PSStream *ps = pswr_new_stream(prefix, wr, TRUE, "doc", NULL, hPageSize*mm, vPageSize*mm);
    /* Cannot free the {prefix} here --- it's in use! */
    /* Sets the page layout to a single picture {hSize} by {vSize}: */
    pswr_set_canvas_layout(ps, hSize*mm, vSize*mm, FALSE, mrg*mm, mrg*mm, 0, 1, 1);
    /* Sets the client plot window to {[0 _ hSize] × [0 _ vSize]}, in mm: */
    pswr_new_picture(ps, 0, hSize, 0, vSize); 
    /* Create the {msm_ps_tools_t} object: */
    msm_ps_tools_t *mps = (msm_ps_tools_t *)notnull(malloc(sizeof(msm_ps_tools_t)), "no mem");
    mps->ps = ps;
    mps->size[0] = hSize; mps->size[1] = vSize;
    /* Set default client and device reference window: */
    int ax;
    for (ax = 0; ax < 2; ax++)
      { /* Device ref window is the whole plot area: */
        mps->cMin[ax] = 0; mps->cMax[ax] = mps->size[ax];
        /* Client ref window is set so that cient and device coords are the same: */
        mps->zMin[ax] = mps->cMin[ax]; mps->zMax[ax] = mps->cMax[ax];
      }
    mps->fontSize = fontSize;
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
    int maxLabChars,
    double mrg
  )
  {
    double mm = msm_PT_PER_MM;
    
    /* Estimate the character dimensions in mm: */
    double vCharSize = fontSize/mm; /* Estimated character height, with skosh. */
    double hCharSize = 0.8*vCharSize; /* Reasonable guess? */
    
    /* Compute the max label dimensions in mm: */
    double hLabelSize = maxLabChars*hCharSize;
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
    
    /* Make sure that there is enough space for overshoot of scales: */
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
    msm_ps_tools_t *mps = msm_ps_tools_new(wr, name, tag, hTotSize, vTotSize, fontSize, mrg);
    
    /* Set reference rectangles to graph area only: */
    double hMin = 0 + mrgL, hMax = hTotSize - mrgR;
    double vMin = 0 + mrgB, vMax = vTotSize - mrgT;
    msm_ps_tools_set_device_ref_window(mps, hMin, hMax, vMin, vMax);
    msm_ps_tools_set_client_ref_window(mps, hMin, hMax, vMin, vMax);
    
    return mps;
  }

PSStream *msm_ps_tools_get_ps_stream(msm_ps_tools_t *mps)
  { return mps->ps; }

void msm_ps_tools_get_plot_size(msm_ps_tools_t *mps, double *hSize, double *vSize)
  { *hSize = mps->size[0]; *vSize = mps->size[1]; }
   
void msm_ps_tools_set_device_ref_window
  ( msm_ps_tools_t *mps,
    double hMin, double hMax, 
    double vMin, double vMax
  )
  { mps->cMin[0] = hMin; mps->cMax[0] = hMax;
    mps->cMin[1] = vMin; mps->cMax[1] = vMax;
  }

void msm_ps_tools_set_client_ref_window
  ( msm_ps_tools_t *mps,
    double xMin, double xMax, 
    double yMin, double yMax
  )
  { mps->zMin[0] = xMin; mps->zMax[0] = xMax;
    mps->zMin[1] = yMin; mps->zMax[1] = yMax;
  }

void msm_ps_tools_shrink_device_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  )
  { mps->cMin[0] += lMrg; mps-> cMax[0] -= rMrg;
    mps->cMin[1] += bMrg; mps-> cMax[1] -= tMrg;
  }

void msm_ps_tools_expand_client_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  )
  { mps->zMin[0] -= lMrg; mps-> zMax[0] += rMrg;
    mps->zMin[1] -= bMrg; mps-> zMax[1] += tMrg;
  }

double msm_ps_tools_map_x(msm_ps_tools_t *mps, double x)
  { double s = (x - mps->zMin[0])/(mps->zMax[0] - mps->zMin[0]), r = 1-s;
    return r*mps->cMin[0] + s*mps->cMax[0];
  }
  
double msm_ps_tools_map_y(msm_ps_tools_t *mps, double y)
  { double s = (y - mps->zMin[1])/(mps->zMax[1] - mps->zMin[1]), r = 1-s;
    return r*mps->cMin[1] + s*mps->cMax[1];
  }

double msm_ps_tools_map_coord(msm_ps_tools_t *mps, pswr_axis_t axis, double coord)
  { demand((axis == HOR) || (axis == VER), "bad axis");
    double s = (coord - mps->zMin[axis])/(mps->zMax[axis] - mps->zMin[axis]), r = 1-s;
    return r*mps->cMin[axis] + s*mps->cMax[axis];
  }
  
void msm_ps_tools_map_coords(msm_ps_tools_t *mps, double x, double y, double *h, double *v)
  { (*h) = msm_ps_tools_map_x(mps, x);
    (*v) = msm_ps_tools_map_y(mps, y);
  }

/* DEVICE -> CLIENT */

double msm_ps_tools_unmap_h(msm_ps_tools_t *mps, double h)
  { double s = (h - mps->cMin[0])/(mps->cMax[0] - mps->cMin[0]), r = 1-s;
    return r*mps->zMin[0] + s*mps->zMax[0];
  }
  
double msm_ps_tools_unmap_v(msm_ps_tools_t *mps, double v)
  { double s = (v - mps->cMin[1])/(mps->cMax[1] - mps->cMin[1]), r = 1-s;
    return r*mps->zMin[1] + s*mps->zMax[1];
  }

double msm_ps_tools_unmap_coord(msm_ps_tools_t *mps, pswr_axis_t axis, double coord)
  { demand((axis == HOR) || (axis == VER), "bad axis");
    double s = (coord - mps->cMin[axis])/(mps->cMax[axis] - mps->cMin[axis]), r = 1-s;
    return r*mps->zMin[axis] + s*mps->zMax[axis];
  }
  
void msm_ps_tools_unmap_coords(msm_ps_tools_t *mps, double h, double v, double *x, double *y)
  { (*x) = msm_ps_tools_unmap_h(mps, h);
    (*y) = msm_ps_tools_unmap_v(mps, v);
  }

/* AXES, TICS, ETC */

void msm_ps_tools_draw_segment(msm_ps_tools_t *mps, double xa, double ya, double xb, double yb)
  { /* Get the device coordinate of the points: */
    double ha, va;
    msm_ps_tools_map_coords(mps, xa, ya, &ha, &va);
    double hb, vb;
    msm_ps_tools_map_coords(mps, xb, yb, &hb, &vb);
    /* Plot the segment: */
    pswr_segment(mps->ps, ha, va, hb, vb);
  }

void msm_ps_tools_draw_ref_axis(msm_ps_tools_t *mps, pswr_axis_t axis)
  { pswr_set_pen(mps->ps, 0.25,0.25,0.25, 0.10, 0.0, 0.0);
    /* Get the device coordinate of the requested axis along the other axis: */
    double cpos = msm_ps_tools_map_coord(mps, 1-axis, 0.0);
    pswr_axis(mps->ps, axis, cpos, mps->cMin[axis], mps->cMax[axis]);
  }

void msm_ps_tools_draw_ref_frame(msm_ps_tools_t *mps)
  { pswr_set_pen(mps->ps, 0,0,0, 0.10, 0.0, 0.0);
    pswr_rectangle(mps->ps, mps->cMin[0], mps->cMax[0], mps->cMin[1], mps->cMax[1], FALSE, TRUE);
  }

void msm_ps_tools_draw_plot_frame(msm_ps_tools_t *mps)
  { pswr_set_pen(mps->ps, 0,0,0, 0.10, 0.0, 0.0);
    pswr_rectangle(mps->ps, 0, mps->size[0], 0, mps->size[1], FALSE, TRUE);
  }
  
void msm_ps_tools_draw_tic
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double xt, 
    double yt,
    double ticSize,
    double ticAlign,
    char *label,
    double labAlign
  )
  { /* Compute nominal device coordinates {ht,vt} of tic: */
    double ht = msm_ps_tools_map_x(mps, xt);
    double vt = msm_ps_tools_map_y(mps, yt);
    /* Plot the tic proper: */
    pswr_tic(mps->ps, axis, ht, vt, ticSize, ticAlign);
    if (label != NULL)
      { /* Decide the vertical and horizontal alignment of label: */
        double hAlign = (axis == HOR ? 0.5 : labAlign);
        double vAlign = (axis == VER ? 0.5 : labAlign);
        /* Adjust the position of the label's reference point to account for tic size: */
        double cadj = (1-labAlign-ticAlign)*ticSize; 
        double hb = (axis == HOR ? ht : ht + cadj);
        double vb = (axis == VER ? vt : vt + cadj);
        pswr_label(mps->ps, label, hb, vb, hAlign, vAlign);
      }
  }

void msm_ps_tools_draw_scale
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double pos,
    double ticSize,
    double ticAlign,
    double ticMinDist, 
    double ticMinStep, 
    char *fmt,
    double labAlign,
    double labMinDist,
    double labMinStep
  )
  { 
    /* Oh paranoia, not even Goya could draw'ya: */
    demand((axis == HOR) || (axis == VER), "bad axis");
    
    /* Get the raw span of tics in device coords (mm): */
    double cMin = mps->cMin[axis];
    double cMax = mps->cMax[axis];
    
    /* Swap if needed to ensure {cMin <= cMax}: */
    if (cMin > cMax) { double c = cMin; cMin = cMax; cMax = c; }
    
    /* Leave some space {eps} between range limits and first/last tics: */
    double eps = 0.5; /* In mm. */
    cMin += eps; cMax -= eps;
    
    /* Choose the client coords of minor tics (0 if no tics in range): */
    double ztMin, ztMax; /* Min and max client coords of minor tics. */
    double ztStep; /* Client coord increment between minor tics (0 if no tics). */
    msm_ps_tools_choose_tic_coords
      ( mps, axis, cMin, cMax, ticMinDist, ticMinStep, &ztMin, &ztMax, &ztStep );
    /* fprintf(stderr, "  ztMin = %24.16e ztMax = %24.16e ztStep = %24.16e\n", ztMin, ztMax, ztStep); */
    if (ztStep <= 0) { /* Not enough space for tics: */ return; }
    
    /* Compute number {ntsteps} of minor tic steps: */
    int ntsteps = (int)floor((ztMax/ztStep) - (ztMin/ztStep) + 0.5);
    assert(ntsteps >= 0); 
    
    /* Compute minor-to-major tic ratio {labPer} (0 if no major tics) and phase {labSkp}: */
    int labPer;   /* Minor tic steps in one major tic step. */
    int labSkp; /* Minor tics before the first major tic. */
    if (fmt == NULL)
      { /* Clients wants no major tics: */
        labPer = labSkp = 0;
      }
    else 
      { /* Choose major tics honoring {labSp} and {minStep}: */
        msm_choose_label_coords
          ( mps, axis, ztMin, ztMax, ztStep, labMinDist, labMinStep, &labPer, &labSkp);
      }
      
    /* Set pen to solid thin black: */
    pswr_set_pen(mps->ps, 0,0,0, 0.10, 0.0, 0.0);
    
    int i;
    for (i = 0; i <= ntsteps; i++)
      { /* Get fractions {r:s} of {i} betwen 0 and {ntsteps}: */
        double s = (ntsteps == 0 ? 0.5 : ((double)i)/((double)ntsteps)), r = 1-s;
        /* Compute nominal client coordinate {zt} of tic on the given {axis}: */
        double zt = r*ztMin + s*ztMax;
        /* Nominal client coords of tic mark are {xt,yt}: */
        double xt = (axis == HOR ? zt : pos);
        double yt = (axis == VER ? zt : pos);
        /* Decide which kind of label to draw: */
        if ((labPer > 0) && (imod(i, labPer) == labSkp))
          { /* Draw a major (labeled) tic: */
            assert(fmt != NULL);
            /* Format the text: */
            char *label = NULL;
            asprintf(&label, fmt, zt);
            /* Draw a labeled tic: */
            msm_ps_tools_draw_tic(mps, axis, xt, yt, 1.5*ticSize, ticAlign, label, labAlign);
            free(label);
          }
        else
          { /* Draw a minor (label-less) tic: */
            msm_ps_tools_draw_tic(mps, axis, xt, yt, ticSize, ticAlign, NULL, 0.0);
          }
      }
  }

void msm_ps_tools_choose_tic_coords
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double cMin, 
    double cMax, 
    double minDist, 
    double minStep, 
    double *zMinP,
    double *zMaxP,
    double *zStepP
  )
  { /* Check required conditions: */
    demand((minDist > 0) || (minStep > 0), "neither minDist nor MinStep were given"); 
    /* Local copies of results: */
    double zMin, zMax, zStep;
    /* Any space at all? */
    if (cMin > cMax) 
      { /* Interval is empty, no tics there: */
        zMin = +INF; zMax = -INF; zStep = 0;
      }
    else
      { /* Compute the client increment {zStep} equivalent to device increment {labSp}: */
        double zDif = fabs(mps->zMax[axis] - mps->zMin[axis]);
        double cDif = fabs(mps->cMax[axis] - mps->cMin[axis]);
        zStep = minDist*(zDif/cDif);
        if (zStep < minStep) { zStep = minStep; }
        assert(zStep >= 0);
        /* Round {zStep} to a nice value: */
        zStep = pswr_round_to_nice(zStep);
        /* Get range {[zMin_zMax]} of client coords corresponding to {[cMin_cMax]}: */
        zMin = msm_ps_tools_unmap_coord(mps, axis, cMin);
        zMax = msm_ps_tools_unmap_coord(mps, axis, cMax);
        if (zMin > zMax) { double z = zMin; zMin = zMax; zMax = z; }
        /* fprintf(stderr, "  zMin =  %24.16e zMax =  %24.16e  zStep = %24.16e\n", zMin, zMax, zStep); */
        /* See whether there are any multiples of {zStep} in {zMin,zMax}: */
        assert(zMin <= zMax);
        double qMin = ceil(zMin/zStep);
        double qMax = floor(zMax/zStep);
        /* fprintf(stderr, "  qMin =  %24.16e qMax =  %24.16e\n", qMin, qMax); */
        if (qMin <= qMax)
          { /* Synchronize {zMin,zMax} to multiples of {zStep}: */
            double zrMin = qMin*zStep; assert(zrMin >= zMin); 
            double zrMax = qMax*zStep; assert(zrMax <= zMax);
            /* fprintf(stderr, "  zrMin = %24.16e zrMax = %24.16e\n", zrMin, zrMax); */
            assert(zrMin <= zrMax);
            zMin = zrMin; zMax = zrMax;
          }
        else
          { /* No multiples of {zStep} in range: */
            zMin = +INF; zMax = -INF; zStep = 0;
          }
      }
    /* fprintf(stderr, "  zMin =  %24.16e zMax =  %24.16e  zStep = %24.16e\n", zMin, zMax, zStep); */
    (*zMinP) = zMin; (*zMaxP) = zMax; (*zStepP) = zStep;
  }

void msm_choose_label_coords
  ( msm_ps_tools_t *mps,
    pswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int *labPerP,
    int *labSkpP
  )
  {
    double eps = 1.0e-12; /* Relative fudge factor to compensate for roundoff errors. */

     /* First, we choose the label period {labPer}: */
    int labPer; 
    
    auto int min_mult(int k, double unit, double min);
      /* Increments the positive integer {k} so that {k*unit >= minv},
        allowing for some roundoff noise in {unit}. */
         
    int min_mult(int k, double unit, double minv)
      { unit = (1+eps)*unit;
        if (k*unit < minv)
          { k = (int)ceil(minv/unit);
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
        double ticDist; /* Actual tic spacing {ticDist} in device coords (fudged up): */
        ticDist = fabs(msm_ps_tools_map_coord(mps, axis, ztStep) - msm_ps_tools_map_coord(mps, axis, 0));
        labPer = min_mult(1, ticDist, minDist);
        assert(labPer > 0);

        /* Adjust {labPer} upwards to honor {minStep}: */
        demand(ztStep > 0, "invalid ztStep"); 
        labPer = min_mult(labPer, ztStep, minStep);
        assert(labPer > 0);

        /* Adjust {labPer} so that the label value increment  is nice: */
        if (labPer > 1)
          { /* Compute the label step {zbStep} in client coords: */
            double zbStep = labPer*ztStep;
            /* Increment {labPer} until {zbStep} is nice (should't be too far): */
            while (TRUE)
              { assert(zbStep > minStep);
                double zbNice = pswr_round_to_nice(zbStep);
                /* fprintf(stderr, "    labPer = %d zbStep =  %24.16e", labPer, zbStep); */
                /* fprintf(stderr, " zbNice = %24.16e\n", zbNice); */
                if (zbNice == INF) { /* Overflow, give up: */ labPer = 0; break; }
                if (zbNice <= zbStep) { /* Phew! */ break; }
                int labPerNew = min_mult(labPer, ztStep, zbNice);
                assert(labPerNew > labPer);
                labPer = labPerNew;
                zbStep = labPer*ztStep;
              }
          }
      }
      
    /* Compute the number {labSkp} of minor tics before first major tic: */
    int labSkp;
    if (labPer == 0)
      { /* Just in case: */ labSkp = 0; }
    else
      { /* Compute coords {zbMin,zbMax} of first and last major tics: */
        double zbStep = labPer*ztStep; /* Label spacing in client coords. */
        double qbMin = ceil(ztMin/zbStep);
        double qbMax = floor(ztMax/zbStep);
        /* fprintf(stderr, "  qbMin =  %24.16e qbMax =  %24.16e\n", qbMin, qbMax); */
        if (qbMin <= qbMax)
          { /* There is a multiple of {zbStep} in the range: */
            double zbMin = qbMin*zbStep;
            labSkp = (int)floor((zbMin - ztMin)/ztStep + 0.5);
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
    int n
  )
  { if (n < 2) { return; }
    double n1 = n - 1;
    double h0, v0; /* Previously plotted point. */
    msm_ps_tools_map_coords(mps, xMin, y[0], &h0, &v0);
    int i;
    for (i = 1; i < n; i++)
      { /* Plot step from point {i-1} to point {i}: */
        double s = i/n1, r = 1-s;
        double x1 = r*xMin + s*xMax;
        double h1, v1;
        msm_ps_tools_map_coords(mps, x1, y[i], &h1, &v1);
        pswr_segment(mps->ps, h0, v0, h1, v1);
        h0 = h1; v0 = v1;
      }
  }

void msm_ps_tools_draw_y_dots
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int n,
    double rad,
    bool_t fill,
    bool_t draw
  )
  { if (n < 1) { return; }
    double n1 = n - 1;
    int i;
    for (i = 0; i < n; i++)
      { /* Plot point {i}: */
        double s = i/n1, r = 1-s;
        double x = r*xMin + s*xMax;
        double h, v;
        msm_ps_tools_map_coords(mps, x, y[i], &h, &v);
        pswr_dot(mps->ps, h, v, rad, fill, draw);
      }
  }

void msm_ps_tools_close(msm_ps_tools_t *mps)
  { pswr_close_stream(mps->ps);
    free(mps);
  }

void msm_ps_tools_compute_data_range(int n, int stride, double z[], double *zMinP, double *zMaxP)
  { int i;
    double zMin = +INF; 
    double zMax = -INF;
    for (i = 0; i < n; i++)
      { double zi = z[i*stride];
        if (zi < zMin) { zMin = zi; }
        if (zi > zMax) { zMax = zi; }
      }
    /* Make sure that the range is not empty or trivial: */
    if (zMin == zMax)
      { double zAbs = fabs(zMin);
        if (zAbs == 0) { zAbs = 1.0; }
        zMin -= 0.1*zAbs;
        zMax += 0.1*zAbs;
      }
    (*zMinP) = zMin;
    (*zMaxP) = zMax;
  }

void msm_ps_tools_draw_graphs
  ( msm_ps_tools_t *mps,
    int nc,
    int nd,
    bool_t circ,
    double x[],
    double y[],
    double yMin,
    double yMax
  )
  {
    demand(nd > 0, "can't plot empty graph");

    /* Get the {PSStream} handle for {pswr.h} routines: */
    PSStream *ps = msm_ps_tools_get_ps_stream(mps);

    /* Choose the actual X range {[xPlotMin_xPlotMax]} of graphs: */
    /* If {circ} is true, allow for half-steps at each end: */
    double xPlotMin, xPlotMax;
    if (x == NULL)
      { if (circ)
          { xPlotMin = 0 - 0.5; xPlotMax = nd - 1 + 0.5; }
        else
          { xPlotMin = 0; xPlotMax = nd - 1; }
      }
    else
      { if (circ)
          { demand(x[0] < x[nd], "X period must be non-zero");
            double xLast = x[nd] - x[nd-1]; /* Length of last step. */
            xPlotMin = x[0] - 0.5*xLast;    
            xPlotMax = x[nd-1] + 0.5*xLast; 
          }
        else
          { demand(x[0] <= x[nd-1], "abscissas out of order");
            xPlotMin = x[0];
            xPlotMax = x[nd-1];
          }
      }
    /* Choose the nominal plot X range {[xMin_xMax]}: */
    double xMag = fmax(fabs(xPlotMin), fabs(xPlotMax));
    double xDif = fabs(xPlotMax - xPlotMin);
    double xSkosh = fmax(1.0e-12*xMag, 0.05*xDif); /* Extra skosh. */
    double xMin = xPlotMin - xSkosh;
    double xMax = xPlotMax + xSkosh;
    
    if (yMin >= yMax)
      { /* Given range {[yMin _ yMax]} is empty or trivial, compute it from the data: */
        int ncd = nc*nd;
        msm_ps_tools_compute_data_range(ncd, 1, y, &yMin, &yMax); 
        /* Include 0 in the range, if near enough: */
        double yDif = yMax - yMin;
        if ((yMin > 0) && (yMin < +0.25*yDif)) { yMin = 0; }
        if ((yMax < 0) && (yMax > -0.25*yDif)) { yMax = 0; }
      }
    
    /* Set client reference window, with some skosh: */
    msm_ps_tools_set_client_ref_window(mps, xMin, xMax, yMin, yMax);
  
    /* Shrink device window to leave space for tics and labels: */
    double mm = msm_PT_PER_MM;
    double ticSz = 1.0;  /* In mm. */
    double yScaleWd = 2*mps->fontSize + (ticSz + 2.0)/mm; /* Width of Y scale and tics. */
    double xScaleHt = mps->fontSize + (ticSz + 2.0)/mm; /* Height of X scale and tics. */
    msm_ps_tools_shrink_device_ref_window(mps, yScaleWd, 0.0, xScaleHt, 0.0);
    
    /* Decide whether to plot dots at individual samples: */
    double xStep; /* Average step size in client coords. */
    if (x == NULL)
      { xStep = 1.0; }
    else if (circ)
      { xStep = (x[nd] - x[0])/nd; }
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
    pswr_set_label_font(ps, font, mps->fontSize);
    msm_ps_tools_draw_ref_axis(mps, HOR);
    msm_ps_tools_draw_scale(mps, HOR, yMin, ticSz,1.0, 5.0,1.0, "%.0f",   1.2, 12.0,1.0);
    msm_ps_tools_draw_scale(mps, VER, xMin, ticSz,1.0, 3.0,0.0, "%+.2f",  1.2,  4.0,0.0);

    /* Decide indices of first and last segment to plot. */
    /* Segment {i} extends from point {i-1} to point {i}. */
    /* If circular, plots segments {0..nd}, clipped to {[xPlotMin _ xPlotMax]}. */
    /* If not circular, plots segments {1..nd-1}, no need to clip. */
    int iIni = (circ ? 0 : 1); /* Index of first segment to plot. */
    int iFin = (circ ? nd : nd-1); /* Index of last segment to plot. */

    int c;
    for (c = 0; c < nc; c++)
      { /* Plot channel {c}. */
        /* Set the pen color: */
        pswr_set_pen(ps, R[c], G[c], B[c], 0.25, 0.0, 0.0);
        /* Plot channel {c} of the sequence. */
        
        auto void get_data_point(int ix, double *xi, double *yi);
          /* Obtains the coordinates {*xi,*yi} of data point number {ix},
            taking circularity into account. */
          
        void get_data_point(int ix, double *xi, double *yi)
          { double dx = 0.0;
            if (circ)
              { double xPeriod = (x == NULL ? nd : x[nd] - x[0]);
                while (ix < 0)   { ix += nd; dx -= xPeriod; }
                while (ix >= nd) { ix -= nd; dx += xPeriod; }
              }
            assert((ix >= 0) && (ix < nd));
            (*xi) = (x == NULL ? (double)ix : x[ix]) + dx;
            (*yi) = y[c*nd + ix];
          }
        
        auto void get_plot_point(int ix, double *xp, double *yp);
          /* Obtains the coordinates {*xp,*yp} of plot point number
            {ix}, taking circularity into account. For {ix == -1}
            returns {xp==xPlotMin}, for {ix == nd} returns
            {xp==xPlotMax}; in both cases {yp} is obtained by
            interpolation. In the other cases, returns the data point
            number {ix}. */
          
        void get_plot_point(int ip, double *xp, double *yp)
          { 
            if ((ip >= 0) && (ip < nd))
              { get_data_point(ip, xp, yp); }
            else
              { assert(circ);
              double x0, y0, x1, y1; /* Adjacent data points. */
                if (ip == -1)
                  { (*xp) = xPlotMin;
                    get_data_point(-1, &x0, &y0);
                    get_data_point( 0, &x1, &y1);
                  }
                else if (ip == nd) 
                  { (*xp) = xPlotMax;
                    get_data_point(nd-1, &x0, &y0);
                    get_data_point(nd+0, &x1, &y1);
                  }
                else
                  { assert(FALSE); }
                /* Interpolate between data points {ip} and {ip+1}: */
                assert(x0 <= (*xp));
                assert((*xp) <= x1);  
                double s = (x0 >= x1 ? 1.0 : ((*xp) - x0)/(x1 - x0));
                double r = 1 - s;
                (*yp) = r*y0 + s*y1;
              }
          }
        
        /* Compute first point {h0,v0} of polyline: */
        double x0, y0; /* Previously plotted point (client coords). */
        get_plot_point(iIni-1, &x0, &y0);
        double h0, v0; /* Previously plotted point (plot coords). */
        msm_ps_tools_map_coords(mps, x0, y0, &h0, &v0);
        /* Now plot all steps. */
        int i = iIni;
        while (i <= iFin)
          { /* Plot step from sample {i-1} to sample {i}: */
            double x1, y1;
            get_plot_point(i, &x1, &y1);
            double h1, v1;
            msm_ps_tools_map_coords(mps, x1, y1, &h1, &v1);
            pswr_segment(ps, h0, v0, h1, v1);
            h0 = h1; v0 = v1;
            i++;
          }
        if (show_dots)
          { /* Draw dots at data points: */
            pswr_set_fill_color(ps, R[c], G[c], B[c]);
            for (i = 0; i < nd; i++)
              { double xi, yi;
                get_data_point(i, &xi, &yi);
                double hi, vi;
                msm_ps_tools_map_coords(mps, xi, yi, &hi, &vi);
                pswr_dot(ps, hi, vi, 0.5, TRUE, FALSE);
              }
          }
      }
   
    /* Draw a thin frame in black: */
    msm_ps_tools_draw_ref_frame(mps);
  }

void msm_ps_tools_draw_histogram
  ( msm_ps_tools_t *mps,
    int nd,
    double x[],
    double y[],
    double yMin,
    double yMax
  )
  {
    demand(nd > 0, "can't plot empty histogram");

    /* Get the {PSStream} handle for {pswr.h} routines: */
    PSStream *ps = msm_ps_tools_get_ps_stream(mps);

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
    
    /* Set client reference window, with some skosh: */
    msm_ps_tools_set_client_ref_window(mps, xMin, xMax, yMin, yMax);
  
    /* Shrink device window to leave space for tics and labels: */
    double mm = msm_PT_PER_MM;
    double ticSz = 1.0;  /* In mm. */
    double yScaleWd = 2*mps->fontSize + (ticSz + 2.0)/mm; /* Width of Y scale and tics. */
    double xScaleHt = mps->fontSize + (ticSz + 2.0)/mm; /* Height of X scale and tics. */
    msm_ps_tools_shrink_device_ref_window(mps, yScaleWd, 0.0, xScaleHt, 0.0);
    
    /* Draw axes, tics, labels, etc: */
    char *font = "Times-Roman";
    pswr_set_label_font(ps, font, mps->fontSize);
    msm_ps_tools_draw_ref_axis(mps, HOR);
    msm_ps_tools_draw_scale(mps, HOR, yMin, ticSz,1.0, 5.0,1.0, "%.0f",   1.2, 12.0,1.0);
    msm_ps_tools_draw_scale(mps, VER, xMin, ticSz,1.0, 3.0,0.0, "%+.2f",  1.2,  4.0,0.0);

    /* Set the pen color: */
    pswr_set_pen(ps, 1.00,0.00,0.00, 0.25, 0.0, 0.0);

    /* Prepare for first half-step: */
    double x0 = xPlotMin, y0 = 0.0;
    double h0, v0; /* End of previous step. */
    msm_ps_tools_map_coords(mps, x0, y0, &h0, &v0);

    /* Now plot all steps, including half-steps at each end ({i==-1} and {i==nd}): */
    int i = -1;
    while (i <= nd)
      { /* Get upper right corner {x1,y1} of next bar: */
        double x1 = (i >= nd ? xPlotMax : (x == NULL ? i+0.5 : x[i+1]));
        double y1 = (i < 0 ? 0.0: (i >= nd ? 0.0 : y[i]));
        /* Map to plot coordinates {h1,v1}: */
        double h1, v1;
        msm_ps_tools_map_coords(mps, x1, y1, &h1, &v1);
        
        if (i >= 0)
          { /* Plot left wall of bar {i}: */
            pswr_segment(ps, h0, v0, h0, v1);
          }
        /* Plot horizontal part of bar {i}: */
        pswr_segment(ps, h0, v1, h1, v1);
        h0 = h1; v0 = v1;
        i++;
      }
   
    /* Draw a thin frame in black: */
    msm_ps_tools_draw_ref_frame(mps);
  }
