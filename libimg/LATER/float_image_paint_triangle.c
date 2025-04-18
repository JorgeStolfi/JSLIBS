/* See {float_image_paint.h}. */
/* Last edited on 2024-12-04 23:21:44 by stolfi */

#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <r2.h>
#include <jsmath.h>
#include <affirm.h>

#include <float_image.h>
#include <float_image_func.h>
#include <float_image_paint.h>

#define Pr fprintf
#define Er stderr

double float_image_paint_triangle
  ( float_image_t *A, 
    int32_t c,           /* Channel. */
    double xa, double ya,    /* Triangle corner. */
    double xb, double yb,    /* Triangle corner. */
    double xc, double ycc,   /* Triangle corner. */
    double hwd,              /* Radius of pen tip. */
    float vfill,             /* Ink value for filling. */                              
    float vdraw,             /* Ink value for stroking. */  
    int32_t m                    /* Subsampling parameter. */
  )
  {
    /* If the drawing ink is invisible, ignore the line width: */
    if (isnan(vdraw)) { hwd = 0.0; }
    
    /* If no fill and no draw, there is nothing to do: */
    if (isnan(vfill) && (hwd <= 0)) { /* Nothing to do: */ return 0.0; }
    
    int32_t nx = (int32_t)A->sz[1];
    int32_t ny = (int32_t)A->sz[2];
    
    /* The sign of {hwd} is irrelevant: */
    hwd = fabs(hwd);
    
    /* Determine the indices of affected columns and rows: */
    double xmin = fmin(xa, fmin(xb, xc));
    double xmax = fmax(xa, fmax(xb, xc));
    double ymin = fmin(ya, fmin(yb, yc));
    double ymax = fmax(ya, fmax(yb, yc));
    int32_t xLo = (int32_t)imax((int32_t)floor(xmin - hwd), 0);
    int32_t xHi = (int32_t)imin((int32_t)floor(xmax + hwd), nx-1);
    int32_t yLo = (int32_t)imax((int32_t)floor(ymin - hwd), 0);
    int32_t yHi = (int32_t)imin((int32_t)floor(ymax + hwd), ny-1);
    if ((xLo > xHi) || (yLo > yHi)) { return 0.0; }
    
    /* Determine the equations of the sides: */
    double Wab = xa*yb - xb*ya, Xab = ya - yb, Yab = xb - xa;
    double Wbc = xb*yc - xc*yb, Xbc = yb - yc, Ybc = xc - xb;
    double Wca = xc*ya - xa*yc, Xca = yc - ya, Yca = xa - xc;
    
    /* Determine the sign of the determinant: */
    

    auto float func_rect(double x, double y);

    float func_tri(double x, double y)
      { if ((x >= xmin+hwd) && (x <= xmax-hwd) && (y >= ymin+hwd) && (y <= ymax-hwd))
          { return vfill; }
        else if ((x >= xmin-hwd) && (x <= xmax+hwd) && (y >= ymin-hwd) && (y <= ymax+hwd))
          { return vdraw; }
        else 
          { return NAN; }
      }

    double w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_rect, NULL, m);
    return w_tot;
  }
  
