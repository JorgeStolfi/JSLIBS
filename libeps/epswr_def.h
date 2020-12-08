/* Private definitions for epswr.h. */
/* Last edited on 2020-10-27 18:51:37 by jstolfi */

#ifndef epswr_def_H
#define epswr_def_H

#include <stdio.h> 

typedef struct epswr_def_figure_t  
  { /* LOCAL ATTRIBUTES NOT AVAILABLE IN THE POSTSCRIPT FILE */
    
    /* Miscellaneous: */
    FILE *wr;              /* EPS file handle. */
    int verbose;           /* TRUE prints a trace of some commands. */
    
    /* Client plotting window: */
    double xMin, xMax;     /* Client X plotting range (in Client coords). */
    double yMin, yMax;     /* Client Y plotting range (in Client coords). */
    
    /* Drawing parameters: */
    double fillColor[3];   /* Fill color for subsequent graphics ops. */
    
    /* Text layout parameters: */
    double hCtrText;       /* Device {h} coordinate of center of nominal text rectangle. */
    double vCtrText;       /* Device {v} coordinate of center of nominal text rectangle. */
    double hSizeText;      /* Width of text area before rotation. */
    double rotText;        /* Rotation of text area (degrees, counterclockwise).*/
    double vTopText;       /* Vert position of top of free area (bottom of prev text line) rel center. */

    /* ATTRIBUTES THAT MIRROR POSTSCRIPT FILE VARIABLES: */
    
    /* Total figure size (Postscript variables {hSize,vSize}: */
    double hSize;          /* Total figure width (in pt). */ /* !!! Should be mirrored !!! */
    double vSize;          /* Total figure height (in pt). */ /* !!! Should be mirrored !!! */
    
    /* Device plotting window (Postscript variables {hMin,hMax,vMin,vMax}: */
    double hMin, hMax;     /* X plotting range (in pt, Device coords). */
    double vMin, vMax;     /* Y plotting range (in pt, Device coords). */
    
    /* Label font parameters (Postscript variable {labelFont}, includes scaling): */
    const char *labelFont; /* Name of font for labels. */
    double labelFontSize;  /* Nominal height of font for labels (in pt). */
    
    /* !!! Make a copy of font names, release at end !!! */
    
    /* Text font parameters (Postscript variable {textFont}, includes scaling): */
    const char *textFont;  /* Name of font for text.  */
    double textFontSize;   /* Nominal height of font for text (in pt). */
    
    /* OTHER INFORMATION: */

    int nFonts;            /* Number of fonts used in file. */
    char **fonts;          /* Fonts used (alloc size is {2^ceil(log2(nFonts))}. */

  } epswr_def_figure_t;
  /* An {epswr_def_figure_t} is a handle to a Postscript file plus some
    state. The fields {hSize,vSize} give the total dimensions of the
    Encapsulated Postscript figure, including margins and
    captions. */

void epswr_def_set_device_window
  ( epswr_def_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    const char *requestor
  );
  /* Saves the given data to the fields {.hMin, .hMax, .vMin, .vMax}
    of {eps} (the current plotting window in Device coordinates).
    
    If {eps->verbose} is true, prints the values to {stderr}. The {requestor} argument 
    should be the name of the function that called this one.
    
    WARNING: Does not change the Client coordinates of plot window, which
    may become invalid. */

void epswr_def_undefine_client_window(epswr_def_figure_t *eps);
  /* Sets the Client coordinates of the plot window ({.xMin,.xMax,.yMin,.yMax})
    to all {NAN}s. */

void epswr_def_reset_client_window(epswr_def_figure_t *eps);
  /* Sets the Client coordinates{.xMin,.xMax,.yMin,.yMax} of the plot window 
    be {0,dh,0,dv} where {dh} and {dv} are the width and height of the
    plot window in Device coordinates.  This makes Client coordinates
    to be in Device units but relative to the low corner
    of the plot window rather than the low corner of the figure. */

void epswr_def_set_client_window
  ( epswr_def_figure_t *eps,
    double xMin, double xMax,
    double yMin, double yMax,
    const char *requestor
  );
  /* Saves the given data to the fields {.xMin, .xMax, .yMin, .yMax}
    of {eps} (the current plotting window in Client coordinates).
    
    If {eps->verbose} is true, prints the values to {stderr}. The {requestor} argument 
    should be the name of the function that called this one. */

void epswr_def_get_client_window
  ( epswr_def_figure_t *eps,
    double *xMinP, double *xMaxP,
    double *yMinP, double *yMaxP
  );
  /* Sets the variables {*xMinP,*xMaxP,*yMinP,*yMaxP} to the bounds
    of the current plot window (in Client coordinates). */
    
/* DEBUGGING */

void epswr_def_report_window
  ( epswr_def_figure_t *eps,
    const char *requestor,
    const char *which,
    const char *action,
    double hxMin, double hxMax,
    double vyMin, double vyMax
  );
  /* Prints the window {[hxMin _ hxMax] x [vyMin _ vyMax]} to {stderr}.
    The parameter {requestor} should be the name of a function that 
    called this one (directly or indirectly).  The {which} string
    should be "Device" or "Client". The {action} string could be
    "set to", "requested", etc. */

#endif
