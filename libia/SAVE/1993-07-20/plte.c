/* See plte.h */

#include "plte.h"
#include "pstools.h"
#include "foimisc.h"
#include "iomisc.h"
#include <stdio.h>

/*** PROTOTYPES FOR INTERNAL FUNCTIONS ***/

void plte_aux_set_coords(
    FILE *psfile,
    char axis,
    double ulo, double uhi,
    double wlo, double whi
  );
  
void plte_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray,
    char *operator
  );
  
/*** IMPLEMENTATIONS ***/

void plte_aux_set_coords(
    FILE *psfile,
    char axis,
    double ulo, double uhi,
    double wlo, double whi
  )
  {
    double scale = (uhi - ulo)/(whi - wlo);
    pst_set_coords(psfile, axis, ulo, uhi, wlo, whi);
    fprintf(psfile, "/%cmin %.6f def\n", axis, ulo);
    fprintf(psfile, "/%cmax %.6f def\n", axis, uhi);
    fprintf(psfile, "/%cpt %f def\n", axis, scale);
    fprintf(psfile, "/%cin %cpt 72.0 mul def\n", axis, axis);
    fprintf(psfile, "/%cmm %cin 25.4 div def\n", axis, axis);
    fprintf(psfile, "\n");
  }

static double plte_xmin, plte_xmax, plte_xscale;
static double plte_ymin, plte_ymax, plte_yscale;

#define plte_hmin (1.25 * 72.0)
#define plte_hmax (7.25 * 72.0)
#define plte_vmin (3.50 * 72.0)
#define plte_vmax (9.50 * 72.0)

void plte_begin_figure(
    FILE *psfile,
    double xmin, double xmax,
    double ymin, double ymax,
    int xn, int yn
  )
  {
    double fxn = xn;
    double fyn = yn;
    char *date = today();

    fprintf(stderr, "date = %s\n", date);

    fprintf(psfile, "%%!PS-Adobe-2.0\n");
    fprintf(psfile, "%%%%CreationDate: %s\n", date);
    fprintf(psfile, "%%%%BoundingBox: %f %f %f %f\n", 
      plte_hmin, plte_vmin, plte_hmax, plte_vmax
    );
    fprintf(psfile, "%%%%Pages: 0\n");
    fprintf(psfile, "%%%%EndComments\n");
    
    fprintf(psfile, "/$pltedict 6400 dict def \n");
    fprintf(psfile, "$pltedict begin\n");
    fprintf(psfile, "$pltedict /mtrx matrix put\n");
    fprintf(psfile, "/l {lineto} bind def\n");
    fprintf(psfile, "/m {moveto} bind def\n");
    fprintf(psfile, "/s {stroke} bind def\n");
    fprintf(psfile, "/n {newpath} bind def\n");
    fprintf(psfile, "/gs {gsave} bind def\n");
    fprintf(psfile, "/gr {grestore} bind def\n");
    fprintf(psfile, "/clp {closepath} bind def\n");
    fprintf(psfile, "end\n");
    fprintf(psfile, "%%%%EndProlog\n");

    fprintf(psfile, "%%%%Page: 1\n");
    fprintf(psfile, "$pltedict begin\n");
    fprintf(psfile, "/savedstate save def\n");

    fprintf(psfile, "%% Round joints and caps:\n");
    fprintf(psfile, "1 setlinecap 1 setlinejoin\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Black thin lines:\n");
    fprintf(psfile, "0 setlinewidth 0 setgray [ ] 0 setdash\n");
    fprintf(psfile, "\n");

    pst_set_coords(psfile, 'x', 0.0, 1.0, plte_hmin, plte_hmax);
    plte_xmin = xmin;
    plte_xmax = xmax;
    plte_xscale = 1.0/(xmax - xmin);

    fprintf(psfile, "/xmin 0.0 def  %% min plottable x\n");
    fprintf(psfile, "/xmax 1.0 def  %% max plottable x\n");
    fprintf(psfile, "/xn %d def     %% grid cells along x axis\n", xn);
    fprintf(psfile, "/xstep %f def  %% x-size of grid cell\n", 1.0/fxn);

    pst_set_coords(psfile, 'y', 0.0, 1.0, plte_vmin, plte_vmax);
    plte_ymin = ymin;
    plte_ymax = ymax;
    plte_yscale = 1.0/(ymax - ymin);

    fprintf(psfile, "/ymin 0.0 def  %% min plottable y\n");
    fprintf(psfile, "/ymax 1.0 def  %% max plottable y\n");
    fprintf(psfile, "/yn %d def     %% grid cells along y axis\n", yn);
    fprintf(psfile, "/ystep %f def  %% y-size of grid cell\n", 1.0/fyn);

    assert (
      ((plte_hmax - plte_hmin) == (plte_vmax - plte_vmin)),
      "plte_begin_figure: unequal Postscript scales"
    );

    fprintf(psfile, "%% Units of measure:\n");
    fprintf(psfile, "/pt %f def\n", 1.0/(plte_vmax - plte_vmin));
    fprintf(psfile, "/in pt 72.0 mul def \n");
    fprintf(psfile, "/mm pt 72.0 25.4 div mul def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Set clipping path to boundary of plot area:\n");
    fprintf(psfile, "newpath\n");
    fprintf(psfile, "  xmin ymin moveto\n");
    fprintf(psfile, "  xmax ymin lineto\n");
    fprintf(psfile, "  xmax ymax lineto\n");
    fprintf(psfile, "  xmin ymax lineto\n");
    fprintf(psfile, "  xmin ymin lineto\n");
    fprintf(psfile, "clip\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Rectangle draw operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ recd --> \n");
    fprintf(psfile, "/recd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    3 index 2 index moveto\n");
    fprintf(psfile, "    2 index 2 index lineto\n");
    fprintf(psfile, "    2 index 1 index lineto\n");
    fprintf(psfile, "    3 index 1 index lineto\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Rectangle fill operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ /gray/ recf --> \n");
    fprintf(psfile, "/recf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setgray\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    3 index 2 index moveto\n");
    fprintf(psfile, "    2 index 2 index lineto\n");
    fprintf(psfile, "    2 index 1 index lineto\n");
    fprintf(psfile, "    3 index 1 index lineto\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Rectangle fill and stroke operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ /gray/ recfd --> \n");
    fprintf(psfile, "/recfd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setgray\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    3 index 2 index moveto\n");
    fprintf(psfile, "    2 index 2 index lineto\n");
    fprintf(psfile, "    2 index 1 index lineto\n");
    fprintf(psfile, "    3 index 1 index lineto\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    gsave fill grestore\n");
    fprintf(psfile, "    0 setgray\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Cell fill operator:\n");
    fprintf(psfile, "%%   /xi/ /yi/ celf --> \n");
    fprintf(psfile, "/celf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  3 1 roll \n");
    fprintf(psfile, "  exch dup \n");
    fprintf(psfile, "  xstep mul xmin add exch 1 add xstep mul xmin add\n");
    fprintf(psfile, "  3 2 roll dup\n");
    fprintf(psfile, "  ystep mul ymin add exch 1 add ystep mul ymin add\n");
    fprintf(psfile, "  5 4 roll \n");
    fprintf(psfile, "  recf\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Segment draw operator:\n");
    fprintf(psfile, "%%   /xa/ /ya/ /xb/ /yb/ segd --> \n");
    fprintf(psfile, "/segd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    moveto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Triangle fill operator:\n");
    fprintf(psfile, "%%   /xa/ /ya/ /xb/ /yb/ /xc/ /yc/ /gray/ trif --> \n");
    fprintf(psfile, "/trif\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setgray\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    moveto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Draw an X-value grid line:\n");
    fprintf(psfile, "%%   /x/ xgrd --> \n");
    fprintf(psfile, "/xgrd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    dup ymin moveto\n");
    fprintf(psfile, "    ymax lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Draw an Y-value grid line:\n");
    fprintf(psfile, "%%   /y/ ygrd --> \n");
    fprintf(psfile, "/ygrd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    dup xmin exch moveto\n");
    fprintf(psfile, "    xmax exch lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fflush(psfile);
  }

void plte_end_figure(FILE *psfile)
  {
    fprintf(psfile, "savedstate restore\n");
    fprintf(psfile, "%% Now we are back to the standard coord system.\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "end\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void plte_begin_section(FILE *psfile, char *title)
  {
    fprintf(psfile, "%%%s\n", title);
    fprintf(stderr, "[%s]\n", title);
    fflush(psfile);
  }

void plte_end_section(FILE *psfile)
  {
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void plte_draw_frame (FILE *psfile)
  {
    plte_begin_section(psfile, "Draw frame around plot area");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "%% Assumes xmax, xmin, ymax, ymin are defined.\n");
    fprintf(psfile, "  initclip\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "  xmin ymin moveto\n");
    fprintf(psfile, "  xmax ymin lineto\n");
    fprintf(psfile, "  xmax ymax lineto\n");
    fprintf(psfile, "  xmin ymax lineto\n");
    fprintf(psfile, "  xmin ymin lineto\n");
    fprintf(psfile, "  closepath stroke\n");
    fprintf(psfile, "grestore\n");
    plte_end_section(psfile);
  }

void plte_set_pen(
    FILE *psfile,
    double gray,
    double width,
    double dashlength,
    double dashspace
  )
  {
    fprintf(psfile, "%5.3f setgray\n", gray);
    fprintf(psfile, "mm %.3f mul setlinewidth\n", width);
    if ((dashlength == 0.0) | (dashspace == 0.0))
      { fprintf(psfile, "[ ] 0 setdash\n"); }
    else
      { fprintf(psfile,
          "[ mm %.3f mul mm %.3f mul ] 0 setdash\n",
          dashlength, dashspace
        );
      }
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void plte_draw_segment(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb
  )
  {
    double psxa = plte_xscale * (xa - plte_xmin);
    double psya = plte_yscale * (ya - plte_ymin);
    double psxb = plte_xscale * (xb - plte_xmin);
    double psyb = plte_yscale * (yb - plte_ymin);
    fprintf(psfile,
      "%6.4f %6.4f  %6.4f %6.4f segd\n",
      psxa, psya, psxb, psyb
    );
    fflush(psfile);
  }

void plte_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray,
    char *operator
  )
  {
    double psxlo = plte_xscale * (xlo - plte_xmin);
    double psxhi = plte_xscale * (xhi - plte_xmin);
    double psylo = plte_yscale * (ylo - plte_ymin);
    double psyhi = plte_yscale * (yhi - plte_ymin);
    fprintf(psfile, "%6.4f %6.4f  %6.4f %6.4f",
      psxlo, psxhi, psylo, psyhi
    );
    if (gray >= 0.0) fprintf(psfile, "  %5.3f", gray);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
    putc('R', stderr);
  }

void plte_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi
  )
  {
    plte_aux_rectangle(psfile, xlo, xhi, ylo, yhi, -1.0, "recd");
  }

void plte_fill_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  )
  {
    plte_aux_rectangle(psfile, xlo, xhi, ylo, yhi, gray, "recf");
  }

void plte_fill_and_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  )
  {
    plte_aux_rectangle(psfile, xlo, xhi, ylo, yhi, gray, "recfd");
  }

void plte_fill_triangle(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double gray
  )
  {
    double psxa = plte_xscale * (xa - plte_xmin);
    double psya = plte_yscale * (ya - plte_ymin);
    double psxb = plte_xscale * (xb - plte_xmin);
    double psyb = plte_yscale * (yb - plte_ymin);
    double psxc = plte_xscale * (xc - plte_xmin);
    double psyc = plte_yscale * (yc - plte_ymin);
    fprintf(psfile,
      "%6.4f %6.4f  %6.4f %6.4f  %6.4f %6.4f  %5.3f trif\n",
      psxa, psya, psxb, psyb, psxc, psyc, gray
    );
    fflush(psfile);
    putc('V', stderr);
  }

void plte_fill_grid_cell(FILE *psfile, int xi, int yi, double gray)
  {
    fprintf(psfile, "%3d %3d  %5.3f celf\n", xi, yi, gray);
    fflush(psfile);
    putc('o', stderr);
  }

void plte_draw_coord_line (FILE *psfile, char axis, double coord)
  {
    double pscoord;
    if (axis == 'x')
      { pscoord = plte_xscale * (coord - plte_xmin); }
    else if (axis == 'y')
      { pscoord = plte_yscale * (coord - plte_ymin); }
    else
      { error("plte_draw_coord-line: invalid axis"); }
    fprintf(psfile, "%6.4f %cgrd\n", pscoord, axis);
  }

void plte_draw_grid_lines(FILE *psfile)
  {
    fprintf(psfile, "%% Grid lines:\n");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "  initclip\n");
    fprintf(psfile, "  0 1 xn {\n");
    fprintf(psfile, "    xstep mul xmin add xgrd\n");
    fprintf(psfile, "  } for\n");
    fprintf(psfile, "  0 1 yn {\n");
    fprintf(psfile, "    ystep mul ymin add ygrd\n");
    fprintf(psfile, "  } for\n");
    fprintf(psfile, "grestore\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }
