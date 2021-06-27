/* Control schematic of the Mazoni-Zabini IG-UNICAMP multifocus microscope */
/* Last edited on 2021-06-26 18:57:53 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsfile.h>

#include <epswr.h>
#include <epswr_dim.h>

/* All dimensions are in mm. */

void miguec_draw_diagram(char *fname);
  /* Writes an EPS file called {fname} with the diagram.  */

void miguec_draw_box
  ( epswr_figure_t *epsf, 
    double xCtr, double yCtr,
    double xSize, double ySize,
    char *text
  );
  /* Draws on {epsf} a framed rectangle with center at {(xCtr,yCtr)}, 
    and dimensions {(xSize,ySize)}, with the {text} centered in it. */

void miguec_draw_arrow
  ( epswr_figure_t *epsf, 
    int n,
    double xv[], double yv[],
    double xText1, double yText1,
    char *text1,
    double xText2, double yText2,
    char *text2
  );
  /* Draws on {epsf} an arrow with polygonal stem
    with vertices {(xv[i],yv[i])} for {i} in {0..n-1}.
    Then, if {text1} is not {NULL}, draws the {text1} centered at 
    the point {(xText1,yText1)}.  Ditto for {text2,xText2,yText2}. */
  
/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  {
    miguec_draw_diagram("out/control_diagram.eps");

    return 0;
  }

void miguec_draw_diagram(char *fname)
  {
    /* Device coordinate extent (pt): */
    double hSize = 360;
    double vSize = 200;
    double hvMarg = 4.0;

    FILE *wr = open_write(fname, TRUE);
    bool_t verbose = FALSE;
    epswr_figure_t *epsf = epswr_new_figure(wr, hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, verbose);
    epswr_set_label_font(epsf, "Helvetica", 8.0);
    epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
    
    /* Size of boxes(mm): */
    double xBSize = 24.0;    
    double yBSize = 20.0; 

    /* Coordinates of boxes: */
    int xN = 4, yN = 3;
    double xBMin = 0.75*xBSize, yBMin = 0.75*yBSize;   /* Center of lower left box. */
    double xBStep = 2.500*xBSize, yBStep = 2.0*yBSize;
    double xBMax = xBMin + (xN-1)*xBStep,  yBMax = yBMin + (yN-1)*yBStep; /* Center of upper right box. */

    /* Define Client coords of the plot window (mm): */
    double xMin = 0.0, xMax = xBMax + xBMin;
    double yMin = 0.0, yMax = yBMax + yBMin;
    epswr_set_client_window(epsf, xMin, xMax, yMin, yMax);

    miguec_draw_box(epsf, xBMin + 0*xBStep, yBMin + 1*yBStep, xBSize, yBSize, "computer");
    miguec_draw_box(epsf, xBMin + 1*xBStep, yBMin + 1*yBStep, xBSize, yBSize, "Arduino");
    miguec_draw_box(epsf, xBMin + 1*xBStep, yBMin + 2*yBStep, xBSize, yBSize, "camera");
    miguec_draw_box(epsf, xBMin + 2*xBStep, yBMin + 0*yBStep, xBSize, yBSize, "encoder");
    miguec_draw_box(epsf, xBMin + 2*xBStep, yBMin + 1*yBStep, xBSize, yBSize, "A4899");
    miguec_draw_box(epsf, xBMin + 2*xBStep, yBMin + 2*yBStep, xBSize, yBSize, "multiplex");
    miguec_draw_box(epsf, xBMin + 3*xBStep, yBMin + 1*yBStep, xBSize, yBSize, "motor");
    miguec_draw_box(epsf, xBMin + 3*xBStep, yBMin + 2*yBStep, xBSize, yBSize, "lights");

    double xv[10], yv[10];
    double xText, yText;
    double xADisp = xBSize/8, yADisp = yBSize/8; /* Displacement of arrows when there are 2 arrows. */
    double xTDisp = xBSize/6, yTDisp = yBSize/6; /* Displacement of label rel to arrow. */
    
    xv[0] = xBMin + 0*xBStep - xADisp, yv[0] = yBMin + 1*yBStep + 0.5*yBSize;
    xv[1] = xv[0], yv[1] = yBMin + 2*yBStep + yADisp;
    xv[2] = xBMin + 1*xBStep - 0.5*xBSize, yv[2] = yv[1];
    xText = (xv[1] + xv[2])/2; yText = yv[2] + yTDisp;
    miguec_draw_arrow(epsf, 3, xv, yv, xText, yText, "settings", 0,0,NULL);

    xv[2] = xBMin + 0*xBStep + xADisp, yv[2] = yBMin + 1*yBStep + 0.5*yBSize;
    xv[1] = xv[2], yv[1] = yBMin + 2*yBStep - yADisp;
    xv[0] = xBMin + 1*xBStep - 0.5*xBSize, yv[0] = yv[1];
    xText = (xv[1] + xv[0])/2; yText = yv[0] - yTDisp;
    miguec_draw_arrow(epsf, 3, xv, yv, xText, yText, "image", 0,0,NULL);

    xv[0] = xBMin + 0*xBStep + 0.5*xBSize, yv[0] = yBMin + 1*yBStep + yADisp;
    xv[1] = xBMin + 1*xBStep - 0.5*xBSize, yv[1] = yv[0];
    xText = (xv[1] + xv[0])/2; yText = yv[0] + yTDisp;
    miguec_draw_arrow(epsf, 2, xv, yv, xText, yText, "command", 0,0,NULL);

    xv[1]= xBMin + 0*xBStep + 0.5*xBSize, yv[1] = yBMin + 1*yBStep - yADisp;
    xv[0] = xBMin + 1*xBStep - 0.5*xBSize, yv[0] = yv[1];
    xText = (xv[1] + xv[0])/2; yText = yv[0] - yTDisp;
    miguec_draw_arrow(epsf, 2, xv, yv, xText, yText, "position", xText, yText-yTDisp, "done");

    xv[0] = xBMin + 1*xBStep + 0.5*xBSize, yv[0] = yBMin + 1*yBStep - yADisp;
    xv[1] = xBMin + 2*xBStep - 0.5*xBSize, yv[1] = yv[0];
    xText = (xv[1] + xv[0])/2; yText = yv[0] + yTDisp;
    miguec_draw_arrow(epsf, 2, xv, yv, xText, yText, "motor ctl", 0,0,NULL);

    xv[0] = xBMin + 2*xBStep + 0.5*xBSize, yv[0] = yBMin + 1*yBStep;
    xv[1] = xBMin + 3*xBStep - 0.5*xBSize, yv[1] = yv[0];
    xText = (xv[1] + xv[0])/2; yText = yv[0] + yTDisp;
    miguec_draw_arrow(epsf, 2, xv, yv, xText, yText, "motor pwr", 0,0,NULL);

    xv[0] = xBMin + 1*xBStep + 0.5*xBSize, yv[0] = yBMin + 1*yBStep + yADisp;
    xv[1] = xv[0] + 2*xADisp, yv[1] = yv[0];
    xv[2] = xv[1], yv[2] = yBMin + 2*yBStep;
    xv[3] = xBMin + 2*xBStep - 0.5*xBSize, yv[3] = yv[2];
    xText = (xv[2] + xv[3])/2; yText = yv[3] + yTDisp;
    miguec_draw_arrow(epsf, 4, xv, yv, xText, yText, "light ctl", 0,0,NULL);

    int aN = 6; /* Number of light arrows. */
    for (int k = 0; k < aN; k++) 
      { if ((k < 3) || (k == aN-1))
          { double dy = 2*yADisp*(1 - 2*((double)k)/((double)aN-1));
            xv[0] = xBMin + 2*xBStep + 0.5*xBSize, yv[0] = yBMin + 2*yBStep + dy;
            xv[1] = xBMin + 3*xBStep - 0.5*xBSize, yv[1] = yv[0];
            xText = (xv[1] + xv[0])/2; yText = yv[0] + yTDisp;
            char *text;
            if (k == 0)
              { text = "light pwr"; }
            else if (k == aN-1)
              { text = "..."; }
            else 
              { text = NULL; }
            miguec_draw_arrow(epsf, 2, xv, yv, xText, yText, text, 0,0,NULL);
          }
      }

    xv[2] = xBMin + 1*xBStep, yv[2] = yBMin + 1*yBStep - 0.5*yBSize;
    xv[1] = xv[2], yv[1] = yBMin + 0*yBStep;
    xv[0] = xBMin + 2*xBStep - 0.5*xBSize, yv[0] = yv[1];
    xText = (xv[1] + xv[0])/2; yText = yv[0] + yTDisp;
    miguec_draw_arrow(epsf, 3, xv, yv, xText, yText, "position", 0,0,NULL);

   epswr_end_figure(epsf);
  }

void miguec_draw_box
  ( epswr_figure_t *epsf, 
    double xCtr, double yCtr,
    double xSize, double ySize,
    char *text
  )
  {
    if ((text != NULL) && (strlen(text) > 0))
      { double xLo = xCtr - 0.5*xSize, xHi = xCtr + 0.5*xSize; 
        double yLo = yCtr - 0.5*ySize, yHi = yCtr + 0.5*ySize; 
        epswr_set_fill_color(epsf, 1.000,1.000,0.950);
        epswr_rectangle(epsf, xLo, xHi, yLo, yHi, TRUE, TRUE);
        epswr_set_fill_color(epsf, 0.000,0.000,0.000);
        epswr_label(epsf, text, text, xCtr, yCtr, 0.0, TRUE, 0.5, 0.5, TRUE, FALSE);
      }
  }

void miguec_draw_arrow
  ( epswr_figure_t *epsf, 
    int n,
    double xv[], double yv[],
    double xText1, double yText1,
    char *text1,
    double xText2, double yText2,
    char *text2
  )
  {
    epswr_polygon(epsf, FALSE, xv, yv, n, FALSE, TRUE, FALSE);
    epswr_set_fill_color(epsf, 0.000,0.000,0.000);
    epswr_arrowhead(epsf, xv[n-2], yv[n-2], xv[n-1], yv[n-1], 1.0, 2.5,  1.0, TRUE, TRUE);
    if ((text1 != NULL) && (strlen(text1) > 0))
      { epswr_label(epsf, text1, text1, xText1, yText1, 0.0, TRUE, 0.5, 0.5, TRUE, FALSE); }
    if ((text2 != NULL) && (strlen(text2) > 0))
      { epswr_label(epsf, text2, text2, xText2, yText2, 0.0, TRUE, 0.5, 0.5, TRUE, FALSE); }
  }
    
