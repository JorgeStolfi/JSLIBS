/* See {float_image_paint.h}. */
/* Last edited on 2023-04-23 11:28:58 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <r2.h>
#include <jsmath.h>
#include <affirm.h>
#include <ellipse_crs.h>
#include <ellipse_ouv.h>
#include <ellipse_aligned.h>
#include <gauss_bell.h>
#include <float_image.h>
#include <float_image_func.h>
#include <float_image_paint.h>

#define Pr fprintf
#define Er stderr

double float_image_paint_sample
  ( float_image_t *A, 
    int c, 
    int ix, 
    int iy, 
    float_image_func_t *func, 
    float_image_func_t *mask, 
    int m 
  )
  {
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    
    /* Clip the arguments to the valid range: */
    if ((ix < 0) || (ix >= nx)) { return 0.0; }
    if ((iy < 0) || (iy >= ny)) { return 0.0; }
    if (m < 0) { m = 0; }

    /* Get coordinates {xp,yp} of pixel center: */
    double xp = ix + 0.5;
    double yp = iy + 0.5;
        
    /* Sample the new image {func} over the pixel: */
    double sum_w = 0.0;   /* Sum of antialias weights of all subsamples. */
    double sum_wp = 0.0;  /* Sum of antialias weights times mask values. */
    double sum_wpv = 0.0; /* Sum of antialias weights times mask values times paint value. */
    int dx, dy;
    for (dx = -m; dx <= +m; dx++)
      { double xs = ((double)dx)/(m+1); /* X displacement of subpixel sample. */
        double xa = fabs(xs);
        double wx = (xa < 0.5 ? 1 - 2*xa*xa : 2*(1 - xa)*(1 - xa)); /* X weight factor. */
        for (dy = -m; dy <= m; dy++)
          { double ys = ((double)dy)/(m+1); /* Y displacement of subpixel sample. */
            double ya = fabs(ys);
            double wy = (ya < 0.5 ? 1 - 2*ya*ya : 2*(1 - ya)*(1 - ya)); /* Y weight factor. */
            double w = wx*wy; /* Subsample weight. */
            assert(w > 0.0);
            double v = func(xp + xs, yp + ys);  /* New image value. */
            /* Get the overlay's opacity {p} at the sample point: */
            double p = (isnan(v) ? 0.0 : (mask == NULL ? 1.0 : mask(xp + xs, yp + ys)));
            assert(! isnan(p));
            p = fmax(0.0, fmin(1.0, p));
            if (p != 0.0)
              { sum_wp += w*p;
                sum_wpv += w*p*v;
              }
            sum_w += w;
          }
      }
    
    double w_new; /* Opacity of overlay pixel sample. */
    if (sum_wp > 0)
      { /* Overlay is not totally transparent, mix it into the image sample: */
        assert(sum_w > 0.0);
        w_new = sum_wp/sum_w;
        double w_old = 1.0 - w_new;   
        float *sp = float_image_get_sample_address(A, c, ix, iy);
        double v_old = (*sp);
        (*sp) = (float)(w_old * v_old + sum_wpv/sum_w);
      }
    else
      { w_new = 0.0; }
    return w_new;
  }

double float_image_paint_samples
  ( float_image_t *A, 
    int c, 
    int xLo, 
    int xHi,
    int yLo,
    int yHi, 
    float_image_func_t *func, 
    float_image_func_t *mask,
    int m 
  )
  {
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    /* Clip the arguments to the valid range: */
    if (xLo < 0) { xLo = 0; }
    if (xHi >= nx ) { xHi = nx-1; }
    if (yLo < 0) { yLo = 0; }
    if (yHi >= ny ) { yHi = ny-1; }
    if (m < 0) { m = 0; }
    /* Paint each pixel: */
    double w_tot = 0.0; /* Total overlaid opacity. */
    int ix, iy;
    for (ix = xLo; ix <= xHi; ix++)
      { for (iy = yLo; iy <=yHi; iy++)
          { double w = float_image_paint_sample(A, c, ix, iy, func, mask, m);
            w_tot += w;
          }
      }
    return w_tot;
  }

double float_image_paint_dot
  ( float_image_t *A, 
    int c, 
    double xctr, 
    double yctr, 
    double rad,
    double hwd,
    bool_t round, 
    bool_t diagonal, 
    float vfill,
    float vdraw,
    int m
  )
  {
    /* If the drawing ink is invisible, ignore the line width: */
    if (isnan(vdraw)) { hwd = 0.0; }
    
    /* If no fill and no draw, there is nothing to do: */
    if (isnan(vfill) && (hwd <= 0)) { /* Nothing to do: */ return 0.0; }
    
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    
    /* The signs of {rad} and {hwd} are irrelevant: */
    rad = fabs(rad);
    hwd = fabs(hwd);
    
    /* Determine the max extent {rMax} along any axis: */
    double rMax;
    if ((round) || (! diagonal))
      { rMax = rad + hwd; }
    else
      { rMax = rad*M_SQRT2 + hwd; }
    
    /* Determine the indices of affected columns and rows: */
    int xLo, xHi; float_image_func_get_index_range(xctr, rMax, nx, &xLo, &xHi);
    int yLo, yHi; float_image_func_get_index_range(yctr, rMax, ny, &yLo, &yHi);
    double w_tot; /* Total opacity of dot. */
    
    if (round)
      { /* Round dot. */
        double orad = rad + hwd;
        double irad = fmax(0.0, rad - hwd);
        double orad2 = orad*orad; /* Outer radius, squared. */
        double irad2 = irad*irad; /* Inner radius, squared. */ 
        
        auto float func_round(double x, double y);
    
        float func_round(double x, double y)
          { double dx = x - xctr;
            double dy = y - yctr;
            double dr2 = dx*dx + dy*dy;
            if (dr2 < irad2) 
              { return vfill; }
            else if (dr2 < orad2)
              { return vdraw; }
            else
              { return NAN; }
          }
        
        w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_round, NULL, m);
      }
    else
      { /* Square dot. */
        double orad = rad + hwd;            /* Outer half-side. */
        double irad = fmax(0.0, rad - hwd); /* Inner half-side. */
        double hwd2 = hwd*hwd;              /* Corner radius, squared. */
        
        auto float func_square(double x, double y);
    
        float func_square(double x, double y)
          { double dx = x - xctr;
            double dy = y - yctr;
            if (diagonal)
              { /* Rotate by 45 degrees: */
                double tx = M_SQRT1_2 * (dx + dy);
                double ty = M_SQRT1_2 * (dy - dx); 
                dx = tx; dy = ty;
              }
            dx = fabs(dx);
            dy = fabs(dy);
            if ((dx < irad) && (dy < irad))
              { return vfill; }
            else if ((dx > orad) || (dy > orad))
              { return NAN; }
            else if ((dx < rad) || (dy < rad))
              { return vdraw; }
            else 
              { /* Round corner: */
                dx -= rad;
                dy -= rad;
                double dr2 = dx*dx + dy*dy;
                if (dr2 < hwd2) 
                  { return vdraw; }
                else
                  { return NAN; }
              }
          }
        
        w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_square, NULL, m);
      }
    return w_tot;
  }
  
double float_image_paint_smudge
  ( float_image_t *A, 
    int c,          /* Channel. */                                
    double xctr,    /* Center's X coordinate. */                  
    double yctr,    /* Center's Y coordinate. */                  
    double xdev,    /* Standard deviation in X direction. */   
    double ydev,    /* Standard deviation in Y direction. */                      
    float vfill,     /* Ink value at center. */ 
    int m           /* Subsampling parameter. */
  )
  {
    /* If the center ink is invisible, there is nothing to do: */
    if (isnan(vfill)) { return 0.0; }
    
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    
    /* Determine the affected ranges : */
    double max_devs = gauss_bell_BIG_ARG; /* Max devs that is worth considering. */
    int xLo, xHi; float_image_func_get_index_range(xctr, max_devs*xdev, nx, &xLo, &xHi);
    int yLo, yHi; float_image_func_get_index_range(yctr, max_devs*ydev, ny, &yLo, &yHi);
    
    auto float func_unif(double x, double y);
      /* Procedural image with uniform color {vfill}. */

    float func_unif(double x, double y)
      { return vfill; }

    auto float mask_bell(double x, double y);
      /* Procedural mask with gaussian bell shape. */

    float mask_bell(double x, double y)
      { double gx = gauss_bell_eval(x, xctr, xdev);
        double gy = gauss_bell_eval(y, yctr, ydev);
        return (float)(gx*gy);
      }

    double w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_unif, &mask_bell, m);
    return w_tot;
  }

double float_image_paint_cross
  ( float_image_t *A, 
    int c, 
    double xctr, 
    double yctr, 
    double rad,
    bool_t empty,
    double hwd, 
    bool_t diagonal, 
    float vdraw,
    int m
  )
  {
    /* If the drawing ink is invisible, there is nothing to do: */
    if (isnan(vdraw)) { return 0.0; }
    
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    
    /* The signs of {rad} and {hwd} are irrelevant: */
    rad = fabs(rad);
    hwd = fabs(hwd);
    
    /* Determine the max extent {rMax} along any axis: */
    double rMax;
    if (! diagonal)
      { rMax = rad + hwd; }
    else
      { rMax = rad*M_SQRT1_2 + hwd; }
    
    /* Determine the indices of affected columns and rows: */
    int xLo, xHi; float_image_func_get_index_range(xctr, rMax, nx, &xLo, &xHi);
    int yLo, yHi; float_image_func_get_index_range(yctr, rMax, ny, &yLo, &yHi);
    
    double orad = rad + hwd;  /* Outer radius of cross (incuding {hwd}). */
    double hwd2 = hwd*hwd;    /* Tip radius, squared. */
    
    double rlo = (empty ? 0.50*rad : -0.1); /* Inner radius of cross (not counting {hwd}). */
    double irad = rlo - hwd; /* Inner radius of cross (including {hwd}). */

    auto float func_cross(double x, double y);

    float func_cross(double x, double y)
      { double dx = x - xctr;
        double dy = y - yctr;
        if (diagonal)
          { /* Rotate by 45 degrees: */
            double tx = M_SQRT1_2 * (dx + dy);
            double ty = M_SQRT1_2 * (dy - dx); 
            dx = tx; dy = ty;
          }
        dx = fabs(dx);
        dy = fabs(dy);
        if (dx < dy) { double t = dx; dx = dy; dy = t; }
        if ((dx > rlo) && (dx < rad) && (dy < hwd))
          { /* Rectangular part of arm: */
            return vdraw;
          }
        else if ((dx < irad) || (dx > orad) || (dy > hwd))
          { /* Not in arm: */
            return NAN;
          }
        else
          { /* Round inner or outer tip: */
            if (dx >= rad) { dx = dx - rad; } else { dx = rlo - dx; }
            double dr2 = dx*dx + dy*dy;
            if (dr2 < hwd2) 
              { return vdraw; }
            else
              { return NAN; }
          }
      }

    double w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_cross, NULL, m);
    return w_tot;
  }

double float_image_paint_rectangle
  ( float_image_t *A, 
    int c,           /* Channel. */                                
    double xmin,     /* Min X coordinate. */                  
    double xmax,     /* Max X coordinate. */                  
    double ymin,     /* Min Y coordinate. */                  
    double ymax,     /* Max Y coordinate. */                  
    double hwd,      /* Radius of pen tip. */  
    float vfill,     /* Ink value for filling. */                              
    float vdraw,     /* Ink value for stroking. */  
    int m            /* Subsampling parameter. */
  )
  {
    /* If the drawing ink is invisible, ignore the line width: */
    if (isnan(vdraw)) { hwd = 0.0; }
    
    /* If no fill and no draw, there is nothing to do: */
    if (isnan(vfill) && (hwd <= 0)) { /* Nothing to do: */ return 0.0; }
    
    int nx = (int)A->sz[1];
    int ny = (int)A->sz[2];
    
    /* The signs of {rad} and {hwd} are irrelevant: */
    hwd = fabs(hwd);
    
    /* Determine the indices of affected columns and rows: */
    int xLo = (int)imax((int)floor(xmin - hwd), 0);
    int xHi = (int)imin((int)floor(xmax + hwd), nx-1);
    int yLo = (int)imax((int)floor(ymin - hwd), 0);
    int yHi = (int)imin((int)floor(ymax + hwd), ny-1);
    if ((xLo > xHi) || (yLo > yHi)) { return 0.0; }

    auto float func_rect(double x, double y);

    float func_rect(double x, double y)
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

double float_image_paint_ellipse_crs
  ( float_image_t *A,
    int c,            /* Channel. */
    ellipse_crs_t *E, /* Ellipse parameters. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  )
  {
    assert(E->rad >= 0); /* Paranoia: */
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    double w_tot = float_image_paint_ellipse_ouv(A, c, &(E->ctr), &F, hwd, vfill, vdraw, m);
    return w_tot;
  }

double float_image_paint_ellipse_ouv
  ( float_image_t *A,
    int c,            /* Channel. */
    r2_t *ctr,        /* Center coords. */
    ellipse_ouv_t *F, /* Ellipse parameters, relative to center. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  )
  {
    /* If the drawing ink is invisible, ignore the line width: */
    if (isnan(vdraw)) { hwd = 0.0; }
    
    /* If no fill and no draw, there is nothing to do: */
    if (isnan(vfill) && (hwd <= 0)) { /* Nothing to do: */ return 0.0; }

    /* Paranoia: */
    assert((F->a >= 0) && (F->b >= 0));
    
    /* Get bounding box with some skosh: */
    int xLo, xHi, yLo, yHi;
    ellipse_ouv_int_bbox(ctr, F, hwd + 1.0, &xLo, &xHi, &yLo, &yHi);
    
    auto float func_ellip(double x, double y);
      /* Returns {vdraw} if point {x,y} of {A} is 
        close to the border, else {vfill} if inside the ellipse,
        else {NAN}. */
      
    double w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_ellip, NULL, m);
    return w_tot;
    
    /* IMPLEMENTATION OF INTERNAL PROCS */
    
    float func_ellip(double x, double y)
      { 
        /* Compute point coordinates {p} relative to center: */
        r2_t p = (r2_t){{ x - ctr->c[0], y - ctr->c[1] }};
        float ink;
        if (hwd <= 0)
          { /* No outline; just check inside/outside: */
            ink = (ellipse_ouv_inside(F, &p) ? vfill : NAN);
          }
        else
          { /* Has outline; check whether we are close to the border: */
            double f = ellipse_ouv_border_position(F, hwd, &p);
            if (f >= +1.0)
              { ink = NAN; }
            else if (f <= -1.0)
              { ink = vfill; }
            else
              { ink = vdraw; }
          }
        return ink;
      }
  }

double float_image_paint_ellipse_aligned
  ( float_image_t *A,
    int c,            /* Channel. */
    r2_t *ctr,        /* Center coords. */
    r2_t *rad,        /* Radii in X and Y. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  )
  {
    bool_t debug = FALSE;
    
    /* If the drawing ink is invisible, ignore the line width: */
    if (isnan(vdraw)) { hwd = 0.0; }
    
    /* If no fill and no draw, there is nothing to do: */
    if (isnan(vfill) && (hwd <= 0)) { /* Nothing to do: */ return 0.0; }

    double cx = ctr->c[0];
    double cy = ctr->c[1];
    double rx = rad->c[0];
    double ry = rad->c[1];
    
    if (debug) { Pr(Er, "  ctr = ( %24.16e %24.16e ) rad = ( %24.16e %24.16e )\n", cx, cy, rx, ry); }
    
    /* Paranoia: */
    assert((rx >= 0) && (ry >= 0));
    
    /* Get bounding box with some skosh: */
    int xLo, xHi, yLo, yHi;
    ellipse_aligned_int_bbox(cx, cy, rx, ry, hwd + 1.0, &xLo, &xHi, &yLo, &yHi);
    if (debug) { Pr(Er, "  int bbox = {%d..%d}×{%d..%d}\n", xLo, xHi, yLo, yHi); }
    
    auto float func_ellip(double x, double y);
      /* Returns {vdraw} if point {x,y} of {A} is 
        close to the border, else {vfill} if inside the ellipse,
        else {NAN}. */
      
    double w_tot = float_image_paint_samples(A, c, xLo, xHi, yLo, yHi, &func_ellip, NULL, m);
    return w_tot;
    
    /* IMPLEMENTATION OF INTERNAL PROCS */
    
    float func_ellip(double x, double y)
      { 
        /* Compute point coordinates {p} relative to center: */
        double xp = x - cx;
        double yp = y - cy;
        float ink;
        if (hwd <= 0)
          { /* No outline; just check inside/outside: */
            ink = (ellipse_aligned_inside(rx, ry, xp, yp) ? vfill : NAN);
          }
        else
          { /* Has outline; check whether we are close to the border: */
            double f = ellipse_aligned_border_position(rx, ry, hwd, xp, yp);
            if (f >= +1.0)
              { ink = NAN; }
            else if (f <= -1.0)
              { ink = vfill; }
            else
              { ink = vdraw; }
          }
        return ink;
      }
  }
