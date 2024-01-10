/* Private definitions for pswr.h. */
/* Last edited on 2009-08-24 21:56:28 by stolfi */

#ifndef pswr_def_H
#define pswr_def_H

#include <stdio.h>

typedef struct PSStreamPrivate  /* A PSStream is a Postscript file plus some state: */
  { int verbose;           /* TRUE prints a trace of some commands. */
    /* Stream parameters: */ 
    int eps;               /* Kind of Postscript file (TRUE = encapsulated). */
    char *prefix;           /* File name prefix. */
    /* Canvas dimensions: */
    char *paperSize;       /* Paper size ("letter", "a3", etc.) or NULL. */
    double hCanvasSize;    /* Canvas width (in pt). */
    double vCanvasSize;    /* Canvas height (in pt). */
    /* Canvas layout (picture slot dimensions and counts): */
    double hFirst;         /* Left X coord of first slot. */
    double vFirst;         /* Top Y coord of first slot. */
    double hPicSize;       /* Picture width (in pt). */
    double vPicSize;       /* Picture height (in pt). */
    double hPicMargin;     /* Left/right margin for each picture (in pt). */
    double vPicMargin;     /* Top/bottom margin for each picture (in pt). */
    int captionLines;      /* Number of caption lines (including pic name). */
    int hPicCount;         /* Number of pictures in each row. */
    int vPicCount;         /* Number of pictures in each column. */
    /* Derived parameters: */
    double hSlotSize;      /* Total picture slot width, incl. margins. */
    double vSlotSize;      /* Total picture slot height, incl. margins and caption. */
    /* Pagination state: */
    FILE *file;            /* Current PS file, or NULL if none. */
    int curCanv;           /* Current canvas number. */
    int curSlot;           /* Index of current picture slot in canvas. */
    /* Local graphics state (not available in Postscript file): */
    double xMin, xMax;     /* Client X plotting range (in client coords). */
    double yMin, yMax;     /* Client Y plotting range (in client coords). */
    double xScale, yScale; /* Scale factors (actual_pt/client_unit). */
    double fillColor[3];   /* Fill color for subsequent graphics ops. */
    /* These values mirror Postscript variables: */
    double hMin, hMax;     /* Actual X plotting range (in pt, canvas coords). */
    double vMin, vMax;     /* Actual Y plotting range (in pt, canvas coords). */
    int xGridCount;        /* Number of grid cells in X direction. */
    int yGridCount;        /* Number of grid cells in Y direction. */
    /* Font list: */
    int nFonts;            /* Number of fonts used in file. */
    char **fonts;          /* Fonts used (alloc size is {2^ceil(log2(nfonts))}. */
  } PSStreamPrivate;
  /* The fields {hCanvasSize,vCanvasSize} give the total dimensions
    of the physical page, or of the Encapsulated Postscript
    figure (virtual page).
    
    Individual pictures are arranged in each canvas (real or virtual) as
    a grid of {hPicCount} by {vPicCount} rectangular slots. Each slot
    measures {hSlotSize} by {vSlotSize} points, and comprises: a
    plottable area measuring {hPicSize} by {vPicSize}, blank margins
    of width {hPicMargin} and {vPicMargin} all around, plus enough
    space below the picture for {captionLines} lines of 10pt caption
    text.
    
    Slots are numbered from 0 to {hPicCount*vPicCount-1}, in Western reading
    order. After each call to {pswr_begin_picture}, the current picture slot
    is the one with index {curSlot} in the canvas {curCanv}.
    
    If {curSlot} is -1, then a new canvas has been started, but no 
    {pswr_begin_picture} has been issued on it. The plot window is 
    set to the whole canvas. */

#endif
