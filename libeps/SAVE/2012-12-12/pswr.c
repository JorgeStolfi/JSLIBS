/* See pswr.h */
/* Last edited on 2012-12-11 11:42:36 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <pswr.h>
#include <pswr_def.h>
#include <pswr_aux.h>
#include <pswr_vis.h>

#include <bool.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>

#define mm (72.0/25.4)
  /* One millimeter in Postscript points. */

#define INF INFINITY
  /* IEEE plus infinity. */

/* INTERNAL PROTOTYPES */

void pswr_compute_joining_arc_radius
  ( double ptx0, double pty0, 
    double ptx1, double pty1, 
    double ptx2, double pty2, 
    double rad,
    double *arad
  );
  /* Checks whether the corner between points {pt0=(ptx0,pty0)},
    {pt1=(ptx1,pty1)}, and {pt2=(ptx2,pty2)} can be rounded with an arc
    of radius {rad}. If so, sets {*arad} to {rad}. Otherwise sets
    {*arad} to zero. */

/* STREAM CREATION AND FINALIZATION */

PSStream *pswr_new_stream
  ( char *prefix,       /* Output filename prefix. */
    FILE *file,         /* Optional open file handle. */
    bool_t eps,         /* TRUE for EPS figures, FALSE for PS document. */
    char *docName,      /* Document name for PS output. */
    char *paperSize,    /* Paper size ({letter}, {a3}, etc.). */
    double hCanvasSize, /* Total canvas width (in pt). */
    double vCanvasSize  /* Total canvas height (in pt). */
  )
  {
    PSStream *ps = (PSStream *)notnull(malloc(sizeof(PSStream)), "out of mem");
    
    ps->verbose = FALSE; /* For now. */
    ps->eps = eps;
    /* Make a fresh private copy of the {prefix}, or of an empty string: */
    ps->prefix = txtcat((prefix == NULL ? "" : prefix), "");
    ps->file = file;
    
    if (eps || (paperSize == NULL) || (strlen(paperSize) == 0))
      { /* Ignore the {paperSize} string, use the given page dimensions: */
        ps->paperSize = NULL;
      }
    else
      { /* Make a fresh private copy of the {paperSize} string: */
        ps->paperSize = txtcat(paperSize, "");
        /* Get page dimensions from the {paperSize} string: */
        pswr_get_paper_dimensions(paperSize, &(hCanvasSize), &(vCanvasSize));
      }
    if (hCanvasSize <= 0)
      { fprintf(stderr, "hCanvasSize = %6.1f pt\n", hCanvasSize);
        affirm(FALSE, "invalid hCanvasSize");
      }
    if (vCanvasSize <= 0)
      { fprintf(stderr, "vCanvasSize = %6.1f pt\n", vCanvasSize);
        affirm(FALSE, "invalid vCanvasSize");
      }
    ps->hCanvasSize = hCanvasSize;
    ps->vCanvasSize = vCanvasSize;
    
    /* Clear font list: */
    ps->fonts = NULL; ps->nFonts = 0;
    
    if (! eps)
      { /* Ensure file is open: */
        if (ps->file == NULL) 
          { char *fileName = NULL;
            asprintf(&(fileName), "%s%s.ps", ps->prefix, docName);
            ps->file = open_write(fileName, TRUE);
            free(fileName);
          }
        /* Write file preamble: */
        pswr_begin_document(ps);
      }
    
    /* Define the stream's state: */
    ps->curCanv = 0;  
    ps->curSlot = INT_MAX; /* Indicates {between_canvas} state */
    
    /* Default canvas layout: */
    if (eps)
      { pswr_set_canvas_layout(ps, 0, 0, TRUE,  4.0,  4.0, 0, 1, 1); }
    else
      { pswr_set_canvas_layout(ps, 0, 0, TRUE, 72.0, 72.0, 5, 1, 1); }

    return ps;
  }
  
void pswr_close_stream(PSStream *ps)
  { if (pswr_within_canvas(ps)) { pswr_end_canvas(ps); }
    if (! ps->eps) { pswr_end_document(ps); }  
    free(ps->paperSize);
    free(ps->prefix);
    free(ps);
  }

bool_t pswr_is_eps(PSStream *ps)
  {
    return ps->eps; 
  }

/* DEFINING THE PAGE LAYOUT */

void pswr_set_canvas_size
  ( PSStream *ps,       /* Picture stream. */
    double hCanvasSize,   /* Width of canvas (in pt). */
    double vCanvasSize    /* Height of canvas (in pt). */
  )  
  {
    if (ps->eps) 
      { /* If canvas in progress, flush it: */
        if (pswr_within_canvas(ps)) { pswr_end_canvas(ps); }
        ps->hCanvasSize = hCanvasSize;
        ps->vCanvasSize = vCanvasSize;
        pswr_set_canvas_layout(ps, 0, 0, TRUE,  4.0,  4.0, 0, 1, 1);
      }
    else
      { fprintf(stderr, "pswr_set_canvas_size ignored (not EPS stream)\n"); }
  }

void pswr_get_canvas_size
  ( PSStream *ps,       /* Picture stream. */
    double *hCanvasSize,  /* OUT: Width of canvas (in pt). */
    double *vCanvasSize   /* OUT: Height of canvas (in pt). */
  )  
  { *hCanvasSize = ps->hCanvasSize;
    *vCanvasSize = ps->vCanvasSize;
  }

void pswr_set_canvas_layout
  ( PSStream *ps,     /* Picture stream. */
    double hPicSize,      /* Width of each picture (in pt). */
    double vPicSize,      /* Height of each picture (in pt). */
    bool_t adjustPicSize, /* TRUE to fit {hPicSize,vPicSize} to canvas size. */
    double hPicMargin,    /* Left/right margin for each picture (in pt). */
    double vPicMargin,    /* Top/bottom margin for each picture (in pt). */
    int captionLines,     /* Number of caption lines below each picture. */
    int hPicCount,        /* Number of pictures in each row. */   
    int vPicCount         /* Number of pictures in each column. */
  )    
  {  
    /* If canvas in progress, flush it: */
    if (pswr_canvas_dirty(ps)) { pswr_end_canvas(ps); }
    
    /* Get canvas dimensions: */
    double hCanvasSize = ps->hCanvasSize;
    double vCanvasSize = ps->vCanvasSize;
    
    /* Compute and save requested caption space: */
    if (captionLines < 0) { captionLines = 0; }
    double vCapSize = 10.0*((double)captionLines); /* Caption height (in pt). */
    ps->captionLines = captionLines;
    
    /* Compute the top and bottom picture margins (including captions): */
    double tPicMargin = vPicMargin;
    double bPicMargin = vPicMargin + vCapSize;
    
    /* Save picture margins: */
    ps->hPicMargin = hPicMargin;
    ps->vPicMargin = vPicMargin;

    /* Choose default canvas margin: */
    double hCanvasMargin, vCanvasMargin;
    if (ps->eps)
      { hCanvasMargin = vCanvasMargin =  4.0; }
    else
      { hCanvasMargin = vCanvasMargin = 72.0; }
      
    /* Compute area of canvas that is available for picture slots, */
    /* allowing picture margins and bottom caption to invade canvas margins: */
    double hUsableSize, vUsableSize;
    { /* Compute the extra canvas margins (sides, top, bottom) needed beyond the picture margin: */
      double hXCanvasMargin = (hCanvasMargin < hPicMargin ? 0 : hCanvasMargin - hPicMargin);
      double tXCanvasMargin = (vCanvasMargin < tPicMargin ? 0 : vCanvasMargin - tPicMargin);
      double bXCanvasMargin = (vCanvasMargin < bPicMargin ? 0 : vCanvasMargin - bPicMargin);
      hUsableSize = hCanvasSize - 2*hXCanvasMargin;
      vUsableSize = vCanvasSize - tXCanvasMargin - bXCanvasMargin;
    }
    
    /* (Re)compute slot counts and dimensions, if necessary: */
    double hSlotSize, vSlotSize;
    pswr_fix_slot_size(hUsableSize, hPicMargin, 0.0,      &hPicCount, &hPicSize, &hSlotSize);
    pswr_fix_slot_size(vUsableSize, vPicMargin, vCapSize, &vPicCount, &vPicSize, &vSlotSize);

    /* Adjust sizes to fit canvas, if so requested: */
    if (adjustPicSize)
      { /* Area that is available for pictures (minus pic margins and captions): */
        double hPlottableSize = hUsableSize - hPicCount*2*hPicMargin;
        double vPlottableSize = vUsableSize - vPicCount*(tPicMargin + bPicMargin);
        /* Maximum picture magnification factors in each direction: */
        double hMag = hPlottableSize / (hPicCount * hPicSize);
        double vMag = vPlottableSize / (vPicCount * vPicSize);
        /* Picture magnification factor: */
        double mag = (hMag < vMag ? hMag : vMag); 
        /* Adjust picture size: */
        hPicSize *= mag;
        vPicSize *= mag;
        /* Adjust slot size (preserving margins and captions): */
        hSlotSize = hPicSize + 2*hPicMargin;
        vSlotSize = vPicSize + (tPicMargin + bPicMargin);
      }
      
    if (ps->verbose)
      { fprintf(stderr, "-- canvas layout --\n");
        fprintf(stderr, "  canvasSize   = %6.1f × %6.1f\n", hCanvasSize, vCanvasSize);
        fprintf(stderr, "  canvasMargin = %6.1f × %6.1f\n", hCanvasMargin, vCanvasMargin);
        fprintf(stderr, "  slotSize   = %6.1f × %6.1f\n", hSlotSize, vSlotSize);
        fprintf(stderr, "  picSize    = %6.1f × %6.1f\n", hPicSize, vPicSize);
        fprintf(stderr, "  picMargin  = %6.1f × %6.1f\n", hPicMargin, vPicMargin);
        fprintf(stderr, "  capSize    = %6.1f\n", vCapSize);
        fprintf(stderr, "  picCount   = %6d × %6d\n", hPicCount, vPicCount);
      }
      
    /* Check for canvas overflow: */
    if (hSlotSize*hPicCount - 2*hPicMargin > hCanvasSize) 
      { fprintf(stderr, "warning -- pictures too wide for canvas\n"); }
    if (vSlotSize*vPicCount - 2*vPicMargin > vCanvasSize) 
      { fprintf(stderr, "warning -- pictures too tall for canvas\n"); }
    
    /* Compute top left corner of first slot in canvas: */
    double hFirst = (hCanvasSize - hPicCount*hSlotSize)/2.0;
    double hFirstMin = (ps->eps ? 0.0 : 36.0) - hPicMargin;
    if (hFirst < hFirstMin) { hFirst = hFirstMin; }
    double vFirst = vCanvasSize - (vCanvasSize - vPicCount*vSlotSize)/2.0;
    double vFirstMax = vCanvasSize - ((ps->eps ? 0.0 : 36.0) - hPicMargin);
    if (vFirst > vFirstMax) { vFirst = vFirstMax; }

    /* Save layout parameters: */
    ps->hPicSize = hPicSize;
    ps->vPicSize = vPicSize;
    ps->hSlotSize = hSlotSize;
    ps->vSlotSize = vSlotSize;
    ps->hFirst = hFirst;
    ps->vFirst = vFirst;
    ps->hPicCount = hPicCount;
    ps->vPicCount = vPicCount;
  }

/* START OF NEW CANVAS */

void pswr_new_canvas(PSStream *ps, const char *pageName)
  { if (pswr_within_canvas(ps)) { pswr_end_canvas(ps); }
    pswr_begin_canvas(ps, pageName);
    affirm(pswr_within_canvas(ps), "duh?");
    affirm(! pswr_canvas_dirty(ps), "duh?");
  }

void pswr_sync_canvas(PSStream *ps, const char *pageName)
  { if (pswr_canvas_dirty(ps)) { pswr_new_canvas(ps, pageName); }
    affirm(! pswr_canvas_dirty(ps), "duh?");
  }

void pswr_fill_row(PSStream *ps)
  { if ((pswr_within_canvas(ps)) && (ps->curSlot >= 0))
      { int leftSlot = (ps->curSlot / ps->hPicCount) * ps->hPicCount;
        ps->curSlot = leftSlot + ps->hPicCount - 1;
      }
  }

void pswr_fill_canvas(PSStream *ps)
  { if (pswr_within_canvas(ps))
      { ps->curSlot = ps->hPicCount*ps->vPicCount - 1;
        affirm(! pswr_next_slot_available(ps), "duh?");
      }
  }

/* START OF NEW PICTURE */

void pswr_new_picture
  ( PSStream *ps,              /* Postscript picture stream. */
    double xMin, double xMax,  /* Client X plotting range. */
    double yMin, double yMax   /* Client Y plotting range. */
  )
  {
    /* Make sure that there is a current canvas: */
    if (! pswr_within_canvas(ps)) { pswr_begin_canvas(ps, NULL); }
    /* Make sure that there is a next slot in the canvas: */
    if (! pswr_next_slot_available(ps)) 
      { /* Canvas exhausted, start a new one: */
        pswr_end_canvas(ps);
        pswr_begin_canvas(ps, NULL);
        affirm(pswr_next_slot_available(ps), "duh?");
      }

    /* Advance to that slot: */
    int nextSlot = ps->curSlot+1;
    affirm((nextSlot >= 0) && (nextSlot < ps->hPicCount*ps->vPicCount), "duh?");
    /* Get row and column indices of slot: */
    int hIndex = nextSlot % ps->hPicCount;
    int vIndex = nextSlot / ps->hPicCount;
    /* Get actual plot window {[hMin_hMax]×[vMin_vMax]} for this slot: */
    double hPicSize = ps->hPicSize;
    double vPicSize = ps->vPicSize;
    double hMin = ps->hFirst + hIndex*ps->hSlotSize + ps->hPicMargin;
    double vMax = ps->vFirst - vIndex*ps->vSlotSize - ps->vPicMargin;
    double hMax = hMin + hPicSize;
    double vMin = vMax - vPicSize;
    /* Adjust plot rectangle to fit client aspect ratio: */
    double hScale = hPicSize/(xMax - xMin);
    double vScale = vPicSize/(yMax - yMin);
    double rh = hPicSize/2*(hScale <= vScale ? 1.0 : vScale/hScale);
    double rv = vPicSize/2*(vScale <= hScale ? 1.0 : hScale/vScale);
    double hCtr = (hMax + hMin)/2;
    double vCtr = (vMax + vMin)/2;
    pswr_set_window
      ( ps, 
        xMin,    xMax,    yMin,    yMax,
        hCtr-rh, hCtr+rh, vCtr-rv, vCtr+rv
      );
    
    /* Fix the current slot, which was clobbered by {pswr_set_window}: */
    ps->curSlot = nextSlot;
    /* Reset graphics state: */
    pswr_set_pen(ps, 0.0,0.0,0.0,  0.15,  0.0,0.0);
    pswr_set_label_font(ps, "Courier", 10.0);
    if (ps->verbose)
      { fprintf(stderr, "-- plot ranges --\n");
        fprintf
          ( stderr, "  client [%8.4f _ %8.4f] x [%8.4f _ %8.4f]\n",
            xMin,    xMax,    yMin,    yMax
          );
        fprintf
          ( stderr, "  device [%8.4f _ %8.4f] x [%8.4f _ %8.4f]\n",
            hCtr-rh, hCtr+rh, vCtr-rv, vCtr+rv
          );
      }
  }
  
/* PLOTTING WINDOW */

void pswr_set_canvas_window
  ( PSStream *ps,
    double hMin, double hMax,
    double vMin, double vMax
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    /* Tell {pswr_new_picture} that the canvas is all full: */
    ps->curSlot = ps->hPicCount*ps->vPicCount - 1;
    /* Save the window parameters in the {PSStream} record: */
    pswr_save_window_data
      ( ps,
        hMin, hMax, vMin, vMax,
        hMin, hMax, vMin, vMax
      );
    
    /* Write the window setup commands to the Postcript file: */
    FILE *file = ps->file;
    pswr_write_window_setup_cmds(file, hMin, hMax, hMin, hMax);
    fflush(file);

    /* Reset the grid to {1 × 1} cells: */
    pswr_set_grid(ps, 1, 1);
  }

void pswr_set_client_window
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax
  )
  { /* Get the current canvas window: */
    double hMin = ps->hMin, hMax = ps->hMax;
    double vMin = ps->vMin, vMax = ps->vMax;
    /* Save the window parameters in the {PSStream} record: */
    pswr_save_window_data
      ( ps,
        xMin, xMax, yMin, yMax,
        hMin, hMax, vMin, vMax
      );
    /* Write the window setup commands to the Postcript file: */
    FILE *file = ps->file;
    pswr_write_window_setup_cmds(file, hMin, hMax, vMin, vMax);
    fflush(file);
  }

void pswr_set_window
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax,

    double hMin, double hMax,
    double vMin, double vMax
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    /* Tell {pswr_new_picture} that the canvas is all full: */
    ps->curSlot = ps->hPicCount*ps->vPicCount - 1;
    /* Save the window parameters in the {PSStream} record: */
    pswr_save_window_data
      ( ps,
        xMin, xMax, yMin, yMax,
        hMin, hMax, vMin, vMax
      );
    
    /* Write the window setup commands to the Postcript file: */
    FILE *file = ps->file;
    pswr_write_window_setup_cmds(file, hMin, hMax, vMin, vMax);
    fflush(file);

    /* Reset the grid to {1 × 1} cells: */
    pswr_set_grid(ps, 1, 1);
  }

void pswr_set_grid(PSStream *ps, int xn, int yn)
  {
    /* Save the grid parameters in the {PSStream} record: */
    pswr_save_grid_data(ps, xn, yn);
    
    /* Write the grid setup commands in the Postcript file: */
    FILE *file = ps->file;
    pswr_write_grid_setup_cmds(file, xn, yn);
    fflush(file);
  }

/* DRAWING COMMANDS */

void pswr_set_pen
  ( PSStream *ps,
    double R, double G, double B,
    double width,
    double dashLength,
    double dashSpace
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    FILE *file = ps->file;
    demand((!isnan(R)) && (!isnan(G)) && (!isnan(B)), "invalid pen color");
    if (R < 0.0) { R = 0.0; } else if (R > 1.0) { R = 1.0; }
    if (G < 0.0) { G = 0.0; } else if (G > 1.0) { G = 1.0; }
    if (B < 0.0) { B = 0.0; } else if (B > 1.0) { B = 1.0; }
    fprintf(file, "%5.3f %5.3f %5.3f setrgbcolor\n", R, G, B);
    fprintf(file, "%.3f setlinewidth\n", width * mm);
    if ((dashLength == 0.0) | (dashSpace == 0.0))
      { fprintf(file, " [ ] 0 setdash\n"); }
    else
      { fprintf(file,
          " [ %.3f %.3f ] 0 setdash\n",
          dashLength * mm, dashSpace * mm
        );
      }
  }

void pswr_segment
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    double psxa = ps->hMin + ps->xScale * (xa - ps->xMin);
    double psya = ps->vMin + ps->yScale * (ya - ps->yMin);
    double psxb = ps->hMin + ps->xScale * (xb - ps->xMin);
    double psyb = ps->vMin + ps->yScale * (yb - ps->yMin);
    if (pswr_segment_is_invisible(ps, psxa, psya, psxb, psyb))
      { return; }
    FILE *file = ps->file;
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxa, psya, psxb, psyb
    );
  }

void pswr_curve
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    double psxa = ps->hMin + ps->xScale * (xa - ps->xMin);
    double psya = ps->vMin + ps->yScale * (ya - ps->yMin);
    double psxb = ps->hMin + ps->xScale * (xb - ps->xMin);
    double psyb = ps->vMin + ps->yScale * (yb - ps->yMin);
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psxd = ps->hMin + ps->xScale * (xd - ps->xMin);
    double psyd = ps->vMin + ps->yScale * (yd - ps->yMin);
    if (pswr_curve_is_invisible(ps, psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f arcd\n",
      psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd
    );
  }

void pswr_coord_line(PSStream *ps, pswr_axis_t axis, double pos)
  { affirm(pswr_within_canvas(ps), "bad state");
    double pspos;
    if (axis == HOR)
      { pspos = ps->hMin + ps->xScale * (pos - ps->xMin); }
    else if (axis == VER)
      { pspos = ps->vMin + ps->yScale * (pos - ps->yMin); }
    else
      { affirm(FALSE, "invalid axis"); pspos = 0.0; }
    FILE *file = ps->file;
    fprintf(file, "%6.1f %sgrd\n", pspos, (axis == HOR ? "x" : "y"));
  }

void pswr_axis(PSStream *ps, pswr_axis_t axis, double pos, double lo, double hi)
  { 
    affirm(pswr_within_canvas(ps), "bad state");
    if (axis == HOR)
      { pswr_segment(ps, lo, pos, hi, pos); }
    else if (axis == VER)
      { pswr_segment(ps, pos, lo, pos, hi); }
    else
      { affirm(FALSE, "invalid axis"); }
  }

void pswr_grid_lines(PSStream *ps)
  { affirm(pswr_within_canvas(ps), "bad state");
    FILE *file = ps->file;
    fprintf(file, "gridlines\n");
  }

void pswr_frame(PSStream *ps)
  { affirm(pswr_within_canvas(ps), "bad state");
    FILE *file = ps->file;
    fprintf(file, "wframe\n");
  }

/* CLOSED FIGURES */

void pswr_set_fill_color(PSStream *ps, double R, double G, double B)
  { 
    affirm(pswr_within_canvas(ps), "bad state");
    if (isnan(R) || (fabs(R) == INFINITY) || (R < 0.0)) { R = G = B = -1.0; }
    double *fc = ps->fillColor;
    if ((R != fc[0]) || (G != fc[1]) || (B != fc[2]))
      { FILE *file = ps->file;
        fc[0] = R; fc[1] = G; fc[2] = B;
        pswr_write_fill_color_set_cmd(file, fc); 
      }
  }

void pswr_rectangle
  ( PSStream *ps,
    double xlo, double xhi,
    double ylo, double yhi,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    FILE *file = ps->file;
    double psxlo = ps->hMin + ps->xScale * (xlo - ps->xMin);
    double psxhi = ps->hMin + ps->xScale * (xhi - ps->xMin);
    double psylo = ps->vMin + ps->yScale * (ylo - ps->yMin);
    double psyhi = ps->vMin + ps->yScale * (yhi - ps->yMin);
    if (pswr_rectangle_is_invisible(ps, psxlo, psxhi, psylo, psyhi))
      { return; }
    fprintf(file, "%d %d %6.1f %6.1f  %6.1f %6.1f",
      draw, fill, psxlo, psxhi, psylo, psyhi
    );
    fprintf(file, " rec\n");
    fflush(file);
  }

void pswr_triangle
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxa = ps->hMin + ps->xScale * (xa - ps->xMin);
    double psya = ps->vMin + ps->yScale * (ya - ps->yMin);
    double psxb = ps->hMin + ps->xScale * (xb - ps->xMin);
    double psyb = ps->vMin + ps->yScale * (yb - ps->yMin);
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    if (pswr_triangle_is_invisible(ps, psxa, psya, psxb, psyb, psxc, psyc))
      { return; }
    FILE *file = ps->file;
    fprintf(file,
      "%d %d %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f ",
      draw, fill, psxa, psya, psxb, psyb, psxc, psyc
    );
    fprintf(file, " tri\n");
    fflush(file);
  }

void pswr_quadrilateral
  ( PSStream *ps,
    double x00, double y00,
    double x01, double y01,
    double x10, double y10,
    double x11, double y11,
    bool_t fill, bool_t draw
  )
  { if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double x[4], y[4];
    x[0] = x00; y[0] = y00;
    x[1] = x01; y[1] = y01;
    x[2] = x11; y[2] = y11;
    x[3] = x10; y[3] = y10;
    pswr_polygon(ps, TRUE, x, y, 4, fill, draw, TRUE);
  }

void pswr_polygon
  ( PSStream *ps,
    bool_t closed,
    double x[], double y[],
    int n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double *psx = (double *)malloc(n*sizeof(double));
    double *psy = (double *)malloc(n*sizeof(double));
    int i;
    /* Map points to device coordinates: */
    for (i=0; i<n; i++)
      { psx[i] = ps->hMin + ps->xScale * (x[i] - ps->xMin);
        psy[i] = ps->vMin + ps->yScale * (y[i] - ps->yMin);
      }
    if (! pswr_polygon_is_invisible(ps, psx, psy, n))
      { /* Plot them: */
        FILE *file = ps->file;
        fprintf(file, "%d %d %d %d", draw, fill, evenOdd, closed);
        if (n > 6) { fprintf(file, "\n"); }
        /* Write the sides in the reverse order: */
        int nplin = 0; /* Number of points in current line. */
        for (i = n-1; i >= 0; i--)
          { if (nplin >= 6) { fprintf(file, "\n"); nplin = 0; }
            fprintf(file, "  %6.1f %6.1f", psx[i], psy[i]);
          }
        fprintf(file, "  %d", n);
        fprintf(file, " %s\n", "pol");
        fflush(file);
      }
    free(psx); free(psy);
  }

void pswr_rounded_polygon
  ( PSStream *ps,
    bool_t closed,
    double x[], double y[],
    int n,
    double rad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    /* Corners: */
    double *psx = (double *)malloc(n*sizeof(double));
    double *psy = (double *)malloc(n*sizeof(double));
    int i;
    /* Map points and radius to device coordinates: */
    for (i=0; i<n; i++)
      { psx[i] = ps->hMin + ps->xScale * (x[i] - ps->xMin);
        psy[i] = ps->vMin + ps->yScale * (y[i] - ps->yMin);
      }
    double psrad = ps->yScale * rad; /* Assumes equal scales. */
    if (! pswr_polygon_is_invisible(ps, psx, psy, n))
      { /* Compute the radii to use at each corner: */
        double *arad = (double *)malloc(n*sizeof(double));
        for (i=0; i<n; i++)
          { int j = (i + 1) % n;
            int k = (i + 2) % n;
            /* Adjust the rounding radius at corner {j = i+1}: */
            if ((! closed) && ((j == 0) || (j == n-1)))
              { arad[j] = 0; }
            else
              { pswr_compute_joining_arc_radius
                  ( psx[i],psy[i], psx[j],psy[j], psx[k],psy[k], psrad, &(arad[j]) );
              }
          }
        /* Output the plot command: */
        FILE *file = ps->file;
        fprintf(file, "%d %d %d %d", draw, fill, evenOdd, closed);
        if (n > 3) { fprintf(file, "\n"); }
        /* Write the corners and radii in reverse order: */
        int nplin = 0; /* Number of points in current line. */
        int ii;
        for (ii = n; ii >= 0; ii--)
          { if (nplin >= 3) { fprintf(file, "\n"); nplin = 0; }
            i = ii % n;
            fprintf(file, "  %6.1f %6.1f %8.3f", psx[i], psy[i], arad[i]);
            nplin++;
          }
        fprintf(file, "\n");
        /* Write the number of points and the starting point: */
        fprintf(file, "  %d", n);
        double pmdx = (closed ? (psx[0]+psx[n-1])/2 : psx[0]); 
        double pmdy = (closed ? (psy[0]+psy[n-1])/2 : psy[0]); 
        fprintf(file, "  %6.1f %6.1f", pmdx, pmdy);
        fprintf(file, " %s\n", "cirpol");
        fflush(file);
        free(arad);
      }
    free(psx); free(psy);
  }

void pswr_bezier_polygon
  ( PSStream *ps,
    bool_t closed,
    double x[], double y[],
    int n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    int np = 4*n; /* Number of points. */
    double *psx = (double *)malloc(np*sizeof(double));
    double *psy = (double *)malloc(np*sizeof(double));
    int i, j;
    /* Map points to device coordinates: */
    for (i = 0; i < np; i++)
      { psx[i] = ps->hMin + ps->xScale * (x[i] - ps->xMin);
        psy[i] = ps->vMin + ps->yScale * (y[i] - ps->yMin);
      }
    /* 
      The following test assumes that if the straight polygon with
      those {np} vertices is invisible, the Bézier polygon is
      invisible too. This is true for simple visibility tests (such as
      bounding box), but not for more precise ones. To be correct in
      any case, we should checl the convex hull of those points
      instead.
    */
    if (! pswr_polygon_is_invisible(ps, psx, psy, np))
      { /* Plot it: */
        FILE *file = ps->file;
        fprintf(file, "%d %d %d %d", draw, fill, evenOdd, closed);
        fprintf(file, "\n");
        /* Write the arcs in the reverse order: */
        for (i = n-1; i >= 0; i--)
          { int k0 = 4*i;              /* Start of arc number {i} */
            for (j = 0; j < 4; j++)
              { int j1 = (j + 1) % 4;
                double psxi = psx[k0+j1];
                double psyi = psy[k0+j1];
                fprintf(file, "  %6.1f %6.1f", psxi, psyi);
              }
            fprintf(file, "\n");
          }
        fprintf(file, "  %d", n);
        fprintf(file, " %s\n", "bzpol");
        fflush(file);
      }
    free(psx); free(psy);
  }

void pswr_circle
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = ps->yScale * rad;
    if (pswr_circle_is_invisible(ps, psxc, psyc, psrad))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f  %6.2f", draw, fill, psxc, psyc, psrad); 
    fprintf(file, " cir\n");
    fflush(file);
  }

void pswr_lune
  ( PSStream *ps,
    double xc, double yc, double rad, double tilt,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = ps->yScale * rad;
    double pstilt = tilt * 180.0 / M_PI;
    if (pswr_lune_is_invisible(ps, psxc, psyc, psrad, pstilt))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f  %6.1f %6.2f",
      draw, fill, psxc, psyc, psrad, pstilt); 
    fprintf(file, " lun\n");
    fflush(file);
  }

void pswr_slice
  ( PSStream *ps,
    double xc, double yc, double rad, 
    double start, double stop,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = ps->yScale * rad;
    double psstart = start;
    double psstop  = stop;
    if (pswr_slice_is_invisible(ps, psxc, psyc, psrad, psstart, psstop))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f  %6.1f  %6.2f %6.2f",
      draw, fill, psxc, psyc, psrad, psstart, psstop); 
    fprintf(file, " pie\n");
    fflush(file);
  }

/* PLOT MARKS */

void pswr_tic
  ( PSStream *ps, 
    pswr_axis_t axis, 
    double xc, double yc, 
    double ticSize,
    double align 
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    FILE *file = ps->file;
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psxa, psya, psxb, psyb;
    if (axis == HOR)
      { psxa = psxb = psxc; 
        psya = psyc - align*ticSize*mm; 
        psyb = psyc + (1-align)*ticSize*mm; 
      }
    else if (axis == VER)
      { psxa = psxc - align*ticSize*mm; 
        psxb = psxc + (1-align)*ticSize*mm;
        psya = psyb = psyc; 
      }
    else
      { affirm(FALSE, "invalid axis"); psxa = psxb = psya = psyb = 0.0; }
    if (pswr_segment_is_invisible(ps, psxa, psya, psxb, psyb))
      { return; }
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxa, psya, psxb, psyb
    );
  }   
 
void pswr_dot
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = rad * mm;
    if (pswr_circle_is_invisible(ps, psxc, psyc, psrad))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f  %6.1f", draw, fill, psxc, psyc, psrad); 
    fprintf(file, " cir\n");
    fflush(file);
  }

void pswr_cross
  ( PSStream *ps, 
    double xc, double yc, double rad, bool_t diag,
    bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (! draw) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = rad * mm;
    if (pswr_rectangle_is_invisible(ps, psxc-psrad, psxc+psrad, psyc-psrad, psyc+psrad))
      { return; }
    double psxA, psyA, psxB, psyB, psxC, psyC, psxD, psyD;
    if (diag)
      { double d = psrad*0.7071067812;
        psxA = psxc - d; psyA = psyc - d;
        psxB = psxc + d; psyB = psyc + d;
        psxC = psxc - d; psyC = psyc + d;
        psxD = psxc + d; psyD = psyc - d;
      }
    else
      { psxA = psxc - psrad; psyA = psyc;
        psxB = psxc + psrad; psyB = psyc;
        psxC = psxc;         psyC = psyc - psrad;
        psxD = psxc;         psyD = psyc + psrad;
      }
    FILE *file = ps->file;
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyA, psxB, psyB
    );
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyC, psxD, psyD
    );
  }

void pswr_asterisk
  ( PSStream *ps, 
    double xc, double yc, double rad,
    bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (!draw) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psrad = rad * mm;
    double ct = 0.92387953 * psrad;
    double st = 0.38268343 * psrad;
    if (pswr_rectangle_is_invisible(ps, psxc-ct, psxc+ct, psyc-ct, psyc+ct))
      { return; }
    double psxA = psxc - st; double psyA = psyc - ct;
    double psxB = psxc + st; double psyB = psyc + ct;
    double psxC = psxc - ct; double psyC = psyc - st;
    double psxD = psxc + ct; double psyD = psyc + st;
    FILE *file = ps->file;
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyA, psxB, psyB
    );
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyB, psxB, psyA
    );
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyC, psxD, psyD
    );
    fprintf(file,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyD, psxD, psyC
    );
  }

void pswr_square
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double d = 0.7071067812 * rad * mm;
    double psxlo = psxc - d;
    double psxhi = psxc + d;
    double psylo = psyc - d;
    double psyhi = psyc + d;
    if (pswr_rectangle_is_invisible(ps, psxlo, psxhi, psylo, psyhi))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f %6.1f %6.1f", draw, fill, psxlo, psxhi, psylo, psyhi); 
    fprintf(file, " rec\n");
    fflush(file);
  }

void pswr_diamond
  ( PSStream *ps, 
    double xc, double yc,
    double xRad, double yRad,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    double psxrad = xRad * mm;
    double psyrad = yRad * mm;
    double psxlo = psxc - psxrad;
    double psxhi = psxc + psxrad;
    double psylo = psyc - psyrad;
    double psyhi = psyc + psyrad;
    if (pswr_rectangle_is_invisible(ps, psxlo, psxhi, psylo, psyhi))
      { return; }
    FILE *file = ps->file;
    bool_t evenOdd = TRUE;
    bool_t closed = TRUE;
    fprintf(file, "%d %d %d %d", draw, fill, evenOdd, closed);
    /* Write the sides in the reverse order: */
    fprintf(file, "  %6.1f %6.1f", psxlo, psyc); 
    fprintf(file, "  %6.1f %6.1f", psxc,  psylo); 
    fprintf(file, "  %6.1f %6.1f", psxhi, psyc); 
    fprintf(file, "  %6.1f %6.1f", psxc,  psyhi); 
    fprintf(file, "  %d", 4);
    fprintf(file, " %s\n", "pol");
    fflush(file);
  }   

void pswr_arrowhead 
  ( PSStream *ps,
    double xa, double ya, double xb, double yb,
    double width, double length, 
    double fraction,
    bool_t fill, bool_t draw
  )
  {
    affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    
    /* Seg endpoints in device coords: */
    double psxa = ps->xScale * (xa - ps->xMin);
    double psya = ps->yScale * (ya - ps->yMin);
    double psxb = ps->xScale * (xb - ps->xMin);
    double psyb = ps->yScale * (yb - ps->yMin);
    
    /* Unit direction vector: */
    double dx = psxb - psxa;
    double dy = psyb - psya;
    double d = sqrt(dx*dx + dy*dy);
    dx /= d; dy /= d;
    
    /* Arrow dimensions in device units: */
    double psw = width/2.0 * mm;
    double psh = length * mm;
    
    /* Corners of triangle: */
    double noitcarf = 1.0 - fraction;
    double psxt = ps->hMin + noitcarf*psxa + fraction*psxb;
    double psyt = ps->vMin + noitcarf*psya + fraction*psyb;
    double psxu = psxt - psh * dx + psw * dy;
    double psyu = psyt - psh * dy - psw * dx;
    double psxv = psxt - psh * dx - psw * dy;
    double psyv = psyt - psh * dy + psw * dx;
    
    if (pswr_triangle_is_invisible(ps, psxt, psyt, psxu, psyu, psxv, psyv))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f",
      draw, fill, psxt, psyt, psxu, psyu, psxv, psyv); 
    fprintf(file, " tri\n");
    fflush(file);
  } 

void pswr_grid_cell
  ( PSStream *ps, 
    int xi, int yi,
    bool_t fill, bool_t draw
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    FILE *file = ps->file;
    fprintf(file, "%d %d %3d %3d", draw, fill, xi, yi);
    fprintf(file, " cel\n");
    fflush(file);
  }
    
/* TEXT PRINTING */

void pswr_set_label_font(PSStream *ps, const char *font, double size)
  { affirm(pswr_within_canvas(ps), "bad state");
    pswr_add_font(font, &(ps->nFonts), &(ps->fonts));
    FILE *file = ps->file;
    pswr_write_label_font_setup_cmds(file, font, size);
  }

void pswr_label
  ( PSStream *ps, 
    const char *text,
    double x, double y, 
    double xAlign, double yAlign
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    double psx = ps->hMin + ps->xScale * (x - ps->xMin);
    double psy = ps->vMin + ps->yScale * (y - ps->yMin);
    FILE *file = ps->file;
    pswr_write_ps_string(file, text, "\\267");
    fprintf(file, " ");
    fprintf(file, "  %5.3f %5.3f  %6.1f %6.1f  lbsh\n",
      xAlign, yAlign, psx, psy
    );
    fflush(file);
  }

/* Buffer size for {pswr_add_caption}: */
#define pswr_CAPBUFSIZE 50
    
void pswr_add_caption(PSStream *ps, const char *txt, double xAlign)
  { affirm(pswr_within_canvas(ps), "bad state");
    if (ps->captionLines > 0)
      { FILE *file = ps->file;
        char showcmd[pswr_CAPBUFSIZE], nlcmd[pswr_CAPBUFSIZE];
        snprintf(showcmd, pswr_CAPBUFSIZE, " %5.3f capsh", xAlign);
        snprintf(nlcmd, pswr_CAPBUFSIZE, ") %s\nnl (", showcmd);
        fprintf(file, "nl ");
        pswr_write_ps_string(file, txt, nlcmd);
        fprintf(file, showcmd);
        fprintf(file, "\n");
        fflush(file);
      }
  }

/* MISCELLANEOUS */ 

void pswr_set_verbose(PSStream *ps, const bool_t verbose)
  { ps->verbose = verbose; }

void pswr_comment(PSStream *ps, const char *title)
  { FILE *file = ps->file;
    affirm(file != NULL, "no file");
    fprintf(ps->file, "\n%% [%s]\n", title);
    if (ps->verbose) { fprintf(stderr, "[%s]\n", title); }
    fflush(file);
  }
 
void pswr_flush (PSStream *ps)
  { FILE *file = ps->file;
    affirm(file != NULL, "no file");
    fflush(file);
  }

double pswr_round_to_nice(double x)
  { /* Special cases for 0 and infinity: */
    if ((x == 0) || (x == INF)) { return x; }
    /* Remove the sign of {x} and save it in {s}: */
    double s = (x < 0 ? -1 : +1);
    x = fabs(x);
    /* Compensate for rounding errors in the algorithm below: */
    x = 0.9999999 * x;
    /* If {x} is too small, don't bother: */
    if (x < 1.0e-300) { return s*1.0e-300; }
    /* Finds the smallest power {z} of 10 such that {z > x}: */
    double z = 1.0;
    while (z/10.0 > x) { z /= 10.0; }
    while (z <= x ) { z *= 10.0; }
    /* Check for overflow: */
    if (z == INF) { return s * INF; }
    /* Finds a nice multiple of that power of 10: */
    if (x < 0.20 * z)
      { return 0.20 * s * z; }
    else if (x < 0.25 * z)
      { return 0.25 * s * z; }
    else if (x < 0.5 * z)
      { return 0.5 * s * z; }
    else
      { return s * z; }
  }

void pswr_compute_joining_arc_radius
  ( double ptx0, double pty0, 
    double ptx1, double pty1, 
    double ptx2, double pty2, 
    double rad,
    double *arad
  )
  {
    /* Compute vector {ux,uy} and distance {du} from point 1 to point 0: */
    double ux = ptx0 - ptx1, uy = pty0 - pty1;
    double du = hypot(ux,uy);
    /* Compute vector {vx,vy} and distance {dv} from point 1 to point 2: */
    double vx = ptx2 - ptx1, vy = pty2 - pty1;
    double dv = hypot(vx,vy);
    if ((du < 1.0e-8) || (dv < 1.0e-8))
      { /* Side with zero length, no rounding: */
        (*arad) = 0.0;
      }
    else
      { /* Sides with nonzero length. */
        /* Normalize {ux,uy} and {vx,vy}: */
        ux /= du; uy /= du;
        vx /= dv; vy /= dv;
        /* Compute unit vector {wx,wy} of bisector: */
        double wx = ux + vx;
        double wy = uy + vy;
        double dw = hypot(wx,wy);
        if (dw == 0)
          { /* Angle is 180 degrees, return arc with zero length: */
            (*arad) = 0.0;
          }
        else
          { /* Angle is not 180. */
            /* Normalize {wx,wy}: */
            wx /= dw; wy /= dw;
            /* Get cosine {ct} and sine {st} of half-angle: */
            double ct = wx*ux + wy*uy;
            double st = fabs(wx*uy - wy*ux);
            if (st < 1.0e-8)
              { /* Angle is essentially zero or 360. */
                (*arad) = 0.0;
              }
            else
              { /* Get distance {ds} from point 1 to ends of arc: */
                double dm = rad/st;
                double ds = dm * ct;
                double dsmax = fmin(du,dv);
                if (ds > 1.01*dsmax)
                  { /* Length removed by rounding is excessive: */
                    (*arad) = 0.0;
                  }
                else if (ds < 0.99*dsmax)
                  { /* Seems OK: */
                    (*arad) = rad;
                  }
                else 
                  { /* Reduce the radius slightly to ensure some bit of side remains: */
                    (*arad) = 0.99*rad*(dsmax/ds);
                  }
              }
          }
      }
  }
  
/* PAPER SIZES */

void pswr_get_paper_dimensions(const char *paperSize, double *xpt, double *ypt)
  { /* US sizes are defined in inches, which are 72 "pt" exactly: */
    /* ISO paper sizes are defined in "mm", and need conversion: */
    if ((strcmp(paperSize, "letter") == 0) || (strcmp(paperSize, "Letter") == 0))
      { (*xpt) = 612.0; (*ypt) = 792.0; return; /* 8,5 x 11 in, 216 x 279 mm */  }
    if ((strcmp(paperSize, "a4") == 0) || (strcmp(paperSize, "A4") == 0))
      { (*xpt) = 210 * mm; (*ypt) = 297 * mm; return; }
    if ((strcmp(paperSize, "a3") == 0) || (strcmp(paperSize, "A3") == 0))
      { (*xpt) = 297 * mm; (*ypt) = 420 * mm; return; }
    if ((strcmp(paperSize, "legal") == 0) || (strcmp(paperSize, "Legal") == 0))
      { (*xpt) = 612.0; (*ypt) = 1008.0; return; /* 8,5 x 14 in, 216 x 356 mm */ }
    if ((strcmp(paperSize, "executive") == 0) || (strcmp(paperSize, "Executive") == 0))
      { (*xpt) = 540.0; (*ypt) = 720.0; return; /* 7,5 x 10 in, 190 x 254 mm */ }
    if ((strcmp(paperSize, "ledger") == 0) || (strcmp(paperSize, "Ledger") == 0)
      || (strcmp(paperSize, "tabloid") == 0) || (strcmp(paperSize, "Tabloid") == 0))
      { (*xpt) = 792.0; (*ypt) = 1224.0; return; /* 11 x 17 in, 279 x 432 mm */ }
    if ((strcmp(paperSize, "a10") == 0) || (strcmp(paperSize, "A10") == 0))
      { (*xpt) = 26 * mm; (*ypt) = 37 * mm; return; }
    if ((strcmp(paperSize, "a9") == 0) || (strcmp(paperSize, "A9") == 0))
      { (*xpt) = 37 * mm; (*ypt) = 52 * mm; return; }
    if ((strcmp(paperSize, "a8") == 0) || (strcmp(paperSize, "A8") == 0))
      { (*xpt) = 52 * mm; (*ypt) = 74 * mm; return; }
    if ((strcmp(paperSize, "a7") == 0) || (strcmp(paperSize, "A7") == 0))
      { (*xpt) = 74 * mm; (*ypt) = 105 * mm; return; }
    if ((strcmp(paperSize, "a6") == 0) || (strcmp(paperSize, "A6") == 0))
      { (*xpt) = 105 * mm; (*ypt) = 148 * mm; return; }
    if ((strcmp(paperSize, "a5") == 0) || (strcmp(paperSize, "A5") == 0))
      { (*xpt) = 148 * mm; (*ypt) = 210 * mm; return; }
    if ((strcmp(paperSize, "a2") == 0) || (strcmp(paperSize, "A2") == 0))
      { (*xpt) = 420 * mm; (*ypt) = 594 * mm; return; }
    if ((strcmp(paperSize, "a1") == 0) || (strcmp(paperSize, "A1") == 0))
      { (*xpt) = 594 * mm; (*ypt) = 841 * mm; return; }
    if ((strcmp(paperSize, "a0") == 0) || (strcmp(paperSize, "A0") == 0))
      { (*xpt) = 841 * mm; (*ypt) = 1189 * mm; return; }
    if ((strcmp(paperSize, "2a0") == 0) || (strcmp(paperSize, "2A0") == 0))
      { (*xpt) = 1189 * mm; (*ypt) = 1682 * mm; return; }
    if ((strcmp(paperSize, "4a0") == 0) || (strcmp(paperSize, "4A0") == 0))
      { (*xpt) = 1682 * mm; (*ypt) = 2378 * mm; return; }
    affirm(0, "unkown paper size");
  }
