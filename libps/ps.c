/* See ps.h */
/* Last edited on 2012-12-08 23:35:36 by stolfilocal */

#include <ps.h>
#include <ps_compat.h>
#include <pswr.h>
#include <affirm.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* PS STREAM TABLES */

PSStream *ps_new_file
  ( FILE *psFile,
    bool_t eps,
    const char *paperSize,
    double hSize, double vSize
  );
  /* Creates a new PS stream for file {psFile} of the given type,
    and adds it to the stream table. The parameters {hmin,hmax,vmin,vmax}
    are used only for EPS files. */

void ps_delete_file(FILE *psFile);
  /* Closes the stream associated with the given {psFile},
    frees its storage, and deletes it from the stream table. */

/*** IMPLEMENTATIONS ***/

void ps_begin_figure
  ( FILE *psFile,
    double hmin, double hmax,
    double vmin, double vmax
  )
  { PSStream *pswr = ps_new_file(psFile, TRUE, NULL, hmax, vmax);
    /* Oh paranoia, not even Goya could draw ya... */
    affirm(pswr->file == psFile, "bad stream record");
    affirm(pswr->eps, "bad stream type");
    pswr_new_canvas(pswr, NULL);
  }
  
void ps_end_figure(FILE *psFile)
  { ps_delete_file(psFile);
  }

void ps_begin_document(FILE *psFile, const char *paperSize)
  { PSStream *pswr = ps_new_file(psFile, FALSE, paperSize, 0, 0);
    /* Oh paranoia, not even Goya could draw ya... */
    affirm(pswr->file == psFile, "bad stream record"); /* Paranoia... */
  }

void ps_end_document(FILE *psFile)
  { ps_delete_file(psFile);
  }

void ps_begin_page(FILE *psFile, const char *page)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_new_canvas(pswr, page);
  }
  
void ps_end_page(FILE *psFile)
  { /* No-op. */
  }

void ps_set_window
  ( FILE *psFile,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_window(pswr, xmin, xmax, ymin, ymax, hmin, hmax, vmin, vmax);
    pswr_set_grid(pswr, xn, yn);
  }

void ps_get_paper_dimensions(const char *papersize, double *xpt, double *ypt)
  { pswr_get_paper_dimensions(papersize, FALSE, xpt, ypt); }

void ps_comment(FILE *psFile, const char *title)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_comment(pswr, title);
  }

void ps_add_caption(FILE *psFile, const char *txt)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_add_caption(pswr, txt, 0.0);
  }
  
void ps_set_label_font(FILE *psFile, const char *font, double size)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_label_font(pswr, font, size);
  }
  
void ps_put_label
  ( FILE *psFile, 
    const char *text, 
    double x, double y, 
    double xalign, double yalign
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_label(pswr, text, x, y, 0.0, xalign, yalign);
  }

void ps_draw_frame (FILE *psFile)
  { PSStream *pswr = ps_get_stream(psFile);
    ps_comment(psFile, "Draw frame around plot area");
    pswr_frame(pswr);
  }

void ps_set_pen
  ( FILE *psFile,
    double R, double G, double B,
    double width,
    double dashlength,
    double dashspace
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_pen(pswr, R, G, B, width, dashlength, dashspace);
  }

void ps_draw_segment
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_segment(pswr, xa, ya, xb, yb);
  }

void ps_draw_tic
  ( FILE *psFile, 
    char axis, 
    double xc, double yc, 
    double ticsz,
    double align 
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_axis_t psaxis;
    if ((axis == 'x') || (axis == 'X'))
      { psaxis = HOR; }
    else if ((axis == 'y') || (axis == 'Y'))
      { psaxis = VER; }
    else
      { fatalerror("ps_draw_tic: invalid axis"); psaxis = HOR; }
    pswr_tic(pswr, psaxis, xc, yc, ticsz, align);
  }

void ps_draw_curve
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_curve(pswr, xa, ya, xb, yb, xc, yc, xd, yd);
  }
  
void ps_draw_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_rectangle(pswr, xlo, xhi, ylo, yhi, FALSE, TRUE);
  }

void ps_fill_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_rectangle(pswr, xlo, xhi, ylo, yhi, TRUE, FALSE);
  }

void ps_fill_and_draw_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_rectangle(pswr, xlo, xhi, ylo, yhi, TRUE, TRUE);
  }

void ps_draw_polygon
  ( FILE *psFile,
    double x[], double y[],
    int npoints
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_polygon(pswr, TRUE, x, y, npoints, FALSE, TRUE, FALSE);
  }
    
void ps_fill_polygon
  ( FILE *psFile,
    double x[], double y[],
    int npoints,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_polygon(pswr, TRUE, x, y, npoints, TRUE, FALSE, FALSE);
  }
    
void ps_fill_and_draw_polygon
  ( FILE *psFile,
    double x[], double y[],
    int npoints,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_polygon(pswr, TRUE, x, y, npoints, TRUE, TRUE, FALSE);
  }
    
void ps_fill_circle
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_circle(pswr, xc, yc, radius, TRUE, FALSE);
  }
  
void ps_draw_circle
  ( FILE *psFile,
    double xc, double yc, double radius
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_circle(pswr, xc, yc, radius, FALSE, TRUE);
  }
  
void ps_fill_and_draw_circle
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_circle(pswr, xc, yc, radius, TRUE, TRUE);
  }
 
void ps_fill_dot
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_dot(pswr, xc, yc, radius, TRUE, FALSE);
  }
  
void ps_draw_dot
  ( FILE *psFile,
    double xc, double yc, double radius
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_dot(pswr, xc, yc, radius, FALSE, TRUE);
  }
  
void ps_fill_and_draw_dot
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_dot(pswr, xc, yc, radius, TRUE, TRUE);
  }
  
void ps_fill_and_draw_lune
  ( FILE *psFile,
    double xc, double yc, double radius, double tilt,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_lune(pswr, xc, yc, radius, tilt, TRUE, TRUE);
  }
  
void ps_fill_and_draw_slice
  ( FILE *psFile,
    double xc, double yc, double radius, double start, double stop,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_slice(pswr, xc, yc, radius, start, stop, TRUE, TRUE);
  }
  
void ps_fill_triangle
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double R, double G, double B
  )
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_triangle(pswr, xa, ya, xb, yb, xc, yc, TRUE, FALSE);
  }

void ps_fill_grid_cell(FILE *psFile, int xi, int yi, double R, double G, double B)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_set_fill_color(pswr, R, G, B);
    pswr_grid_cell(pswr, xi, yi, TRUE, FALSE);
  }

void ps_draw_coord_line (FILE *psFile, char axis, double coord)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_axis_t psaxis;
    if ((axis == 'x') || (axis == 'X'))
      { psaxis = HOR; }
    else if ((axis == 'y') || (axis == 'Y'))
      { psaxis = VER; }
    else
      { fatalerror("ps_draw_coord-line: invalid axis"); psaxis = HOR; }
    pswr_coord_line(pswr, psaxis, coord);
  }

void ps_draw_grid_lines(FILE *psFile)
  { PSStream *pswr = ps_get_stream(psFile);
    pswr_grid_lines(pswr);
  }
  
/* STREAM RECORD HANDLING */

/* Originally the {PSStream} information was kept in static variables
  internal to this module. That became a problem when some clients 
  needed to create two or more Postscript files in parallel. 
  
  Ideally, we should have changed the {psFile} parameter of exported
  procedures: instead of a naked {FILE *}, it should point to an
  object contaning the Postscript {PSStream} variables and the file
  itself. However, for compatibility with old clients, we cannot
  implement this change now. So, for the time being, we internally map
  {FILE *} pointers to {PSStream} records, through the following tables. */

#define MAXSTATES 10
static int ps_nfiles = 0;              /* Number of currently open PS files */
static PSStream *ps_stream[MAXSTATES]; /* The states of those files. */

PSStream *ps_get_stream(FILE *psFile)
  { { int i;
      for (i = 0; i < ps_nfiles; i++)
        { PSStream *pswri = ps_stream[i];
          if (pswri->file == psFile) { return pswri; }
        }
    }
    affirm(0, "Postscript file was not initialized?");
    return NULL;
  }

PSStream *ps_new_file
  ( FILE *psFile,
    bool_t eps,
    const char *paperSize,
    double hSize, double vSize
  )
  { /* Sanity check: */
    { int i;
      for (i = 0; i < ps_nfiles; i++)
        { PSStream *pswri = ps_stream[i];
          affirm(pswri->file != psFile, "Postscript file initialized twice?");
        }
    }
    /* Allocate record and store in table: */
    affirm(ps_nfiles < MAXSTATES, "too many simultaneous files");
    { PSStream *pswr = 
        pswr_new_stream("ps", psFile, eps, "doc", "letter", FALSE, hSize, vSize);
      ps_stream[ps_nfiles] = pswr;
      ps_nfiles++;
      return pswr;
    }
  }

void ps_delete_file(FILE *psFile)
  { int i;
    /* Find record in table: */
    for (i = 0; i < ps_nfiles; i++)
      { PSStream *pswri = ps_stream[i];
        if (pswri->file == psFile) { break; }
      }
    affirm(i < ps_nfiles, "Postscript file not initialized?");
    /* Remove from table and free record: */
    pswr_close_stream(ps_stream[i]);
    while (i+1 < ps_nfiles)
      { ps_stream[i] = ps_stream[i+1]; i++; }
    ps_nfiles--;
  }
 
