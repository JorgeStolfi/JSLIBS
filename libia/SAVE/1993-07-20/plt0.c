/* See plt0.h */

#include "plt0.h"
#include "pstools.h"
#include "foimisc.h"
#include "iomisc.h"
#include <stdio.h>

/*** PROTOTYPES FOR INTERNAL FUNCTIONS ***/

void plt0_aux_set_coords(
    FILE *psfile,
    char axis,
    double ulo, double uhi,
    double wlo, double whi
  );
  
void plt0_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray,
    char *operator
  );
  
/*** IMPLEMENTATIONS ***/

void plt0_begin_file(FILE *psfile)
  {
    pst_begin_file(psfile);
  }

void plt0_end_file(FILE *psfile)
  {
    pst_end_file(psfile);
  }

void plt0_aux_set_coords(
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

static double plt0_xmin, plt0_xmax, plt0_xscale;
static double plt0_ymin, plt0_ymax, plt0_yscale;

#define plt0_hmin (1.25 * 72.0)
#define plt0_hmax (7.25 * 72.0)
#define plt0_vmin (3.50 * 72.0)
#define plt0_vmax (9.50 * 72.0)

void plt0_begin_page(
    FILE *psfile,
    int page,
    double xmin, double xmax,
    double ymin, double ymax,
    int xn, int yn
  )
  {
    double fxn = xn;
    double fyn = yn;
    pst_begin_page(psfile, page);
    fprintf(psfile, "30 dict begin\n");
    fprintf(psfile, "gsave\n");

    fprintf(psfile, "%% Round joints and caps:\n");
    fprintf(psfile, "1 setlinecap 1 setlinejoin\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Black thin lines:\n");
    fprintf(psfile, "0 setlinewidth 0 setgray [ ] 0 setdash\n");
    fprintf(psfile, "\n");

    pst_set_coords(psfile, 'x', 0.0, 1.0, plt0_hmin, plt0_hmax);
    plt0_xmin = xmin;
    plt0_xmax = xmax;
    plt0_xscale = 1.0/(xmax - xmin);

    fprintf(psfile, "/xmin 0.0 def  %% min plottable x\n");
    fprintf(psfile, "/xmax 1.0 def  %% max plottable x\n");
    fprintf(psfile, "/xn %d def     %% grid cells along x axis\n", xn);
    fprintf(psfile, "/xstep %f def  %% x-size of grid cell\n", 1.0/fxn);

    pst_set_coords(psfile, 'y', 0.0, 1.0, plt0_vmin, plt0_vmax);
    plt0_ymin = ymin;
    plt0_ymax = ymax;
    plt0_yscale = 1.0/(ymax - ymin);

    fprintf(psfile, "/ymin 0.0 def  %% min plottable y\n");
    fprintf(psfile, "/ymax 1.0 def  %% max plottable y\n");
    fprintf(psfile, "/yn %d def     %% grid cells along y axis\n", yn);
    fprintf(psfile, "/ystep %f def  %% y-size of grid cell\n", 1.0/fyn);

    assert (
      ((plt0_hmax - plt0_hmin) == (plt0_vmax - plt0_vmin)),
      "plt0_begin_page: unequal Postscript scales"
    );

    fprintf(psfile, "%% Units of measure:\n");
    fprintf(psfile, "/pt %f def\n", 1.0/(plt0_vmax - plt0_vmin));
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

    fprintf(psfile, "%% Caption text cursor:\n");
    fprintf(psfile, "/xtext xmin def\n");
    fprintf(psfile, "/ytext ymin def\n");
    fprintf(psfile, "/dytext 10 pt mul def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Caption font setup:\n");
    fprintf(psfile, "/Courier findfont\n");
    fprintf(psfile, "dytext scalefont setfont\n");
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

    fprintf(psfile, "%% Operator to move to new caption line:\n");
    fprintf(psfile, "%%   nl --> \n");
    fprintf(psfile, "/nl\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  /ytext ytext dytext sub def\n");
    fprintf(psfile, "  xtext ytext moveto\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Operator to print string at CP without clipping:\n");
    fprintf(psfile, "%%   /s/ shw --> \n");
    fprintf(psfile, "/shw\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave initclip show grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fflush(psfile);
  }

void plt0_end_page(FILE *psfile)
  {
    fprintf(psfile, "grestore\n");
    fprintf(psfile, "%% Now we are back to the standard coord system.\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "end\n");
    fprintf(psfile, "\n");
    pst_end_page(psfile);
  }

void plt0_add_caption(FILE *psfile, char *txt)
  {
    fprintf(psfile, "nl ");
    pst_put_text(psfile, txt, ") shw\nnl (");
    fprintf(psfile, " shw\n");
    fflush(psfile);
  }

void plt0_begin_section(FILE *psfile, char *title)
  {
    fprintf(psfile, "%%%s\n", title);
    fprintf(stderr, "[%s]\n", title);
    fflush(psfile);
  }

void plt0_end_section(FILE *psfile)
  {
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void plt0_draw_frame (FILE *psfile)
  {
    plt0_begin_section(psfile, "Draw frame around plot area");
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
    plt0_end_section(psfile);
  }

void plt0_set_pen(
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

void plt0_draw_segment(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb
  )
  {
    double psxa = plt0_xscale * (xa - plt0_xmin);
    double psya = plt0_yscale * (ya - plt0_ymin);
    double psxb = plt0_xscale * (xb - plt0_xmin);
    double psyb = plt0_yscale * (yb - plt0_ymin);
    fprintf(psfile,
      "%6.4f %6.4f  %6.4f %6.4f segd\n",
      psxa, psya, psxb, psyb
    );
    fflush(psfile);
  }

void plt0_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray,
    char *operator
  )
  {
    double psxlo = plt0_xscale * (xlo - plt0_xmin);
    double psxhi = plt0_xscale * (xhi - plt0_xmin);
    double psylo = plt0_yscale * (ylo - plt0_ymin);
    double psyhi = plt0_yscale * (yhi - plt0_ymin);
    fprintf(psfile, "%6.4f %6.4f  %6.4f %6.4f",
      psxlo, psxhi, psylo, psyhi
    );
    if (gray >= 0.0) fprintf(psfile, "  %5.3f", gray);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
    putc('R', stderr);
  }

void plt0_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi
  )
  {
    plt0_aux_rectangle(psfile, xlo, xhi, ylo, yhi, -1.0, "recd");
  }

void plt0_fill_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  )
  {
    plt0_aux_rectangle(psfile, xlo, xhi, ylo, yhi, gray, "recf");
  }

void plt0_fill_and_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  )
  {
    plt0_aux_rectangle(psfile, xlo, xhi, ylo, yhi, gray, "recfd");
  }

void plt0_fill_triangle(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double gray
  )
  {
    double psxa = plt0_xscale * (xa - plt0_xmin);
    double psya = plt0_yscale * (ya - plt0_ymin);
    double psxb = plt0_xscale * (xb - plt0_xmin);
    double psyb = plt0_yscale * (yb - plt0_ymin);
    double psxc = plt0_xscale * (xc - plt0_xmin);
    double psyc = plt0_yscale * (yc - plt0_ymin);
    fprintf(psfile,
      "%6.4f %6.4f  %6.4f %6.4f  %6.4f %6.4f  %5.3f trif\n",
      psxa, psya, psxb, psyb, psxc, psyc, gray
    );
    fflush(psfile);
    putc('V', stderr);
  }

void plt0_fill_grid_cell(FILE *psfile, int xi, int yi, double gray)
  {
    fprintf(psfile, "%3d %3d  %5.3f celf\n", xi, yi, gray);
    fflush(psfile);
    putc('o', stderr);
  }

void plt0_draw_coord_line (FILE *psfile, char axis, double coord)
  {
    double pscoord;
    if (axis == 'x')
      { pscoord = plt0_xscale * (coord - plt0_xmin); }
    else if (axis == 'y')
      { pscoord = plt0_yscale * (coord - plt0_ymin); }
    else
      { error("plt0_draw_coord-line: invalid axis"); }
    fprintf(psfile, "%6.4f %cgrd\n", pscoord, axis);
  }

void plt0_draw_grid_lines(FILE *psfile)
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
