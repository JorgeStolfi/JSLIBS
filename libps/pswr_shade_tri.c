/* See pswr_shade_tri.h */
/* Last edited on 2009-08-25 01:07:58 by stolfi */

#define _GNU_SOURCE
#include <assert.h>

#include <pswr.h>
#include <pswr_aux.h>
#include <pswr_vis.h>
#include <bool.h>
#include <affirm.h>

#include <pswr_shade_tri.h>

void pswr_shade_triangle_subd
  ( PSStream *ps,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  );
  /* Same as pswr_shade_triangle, for {ns > 0}. */

void pswr_shade_triangle_gouraud
  ( PSStream *ps,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc
  );
  /* Same as pswr_shade_triangle, but uses true Gouraud shading. */
  
void pswr_shade_triangle
  ( PSStream *ps,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  )
  { if (ns < 0)
      { /* Gouraud-shaded triangle: */
        pswr_shade_triangle_gouraud
          ( ps,
            xa, ya, Ra, Ga, Ba,
            xb, yb, Rb, Gb, Bb,
            xc, yc, Rc, Gc, Bc
          );
      }
    else if (ns == 0)
      { /* Solid triangle: */
        double Rt = (Ra + Rb + Rc)/3;
        double Gt = (Ga + Gb + Gc)/3;
        double Bt = (Ba + Bb + Bc)/3;
        pswr_set_fill_color(ps, Rt, Gt, Bt);
        pswr_triangle(ps, xa, ya, xb, yb, xc, yc, TRUE, FALSE);
      }
    else
      { /* Shading by subdivision: */
        pswr_shade_triangle_subd
          ( ps,
            xa, ya, Ra, Ga, Ba,
            xb, yb, Rb, Gb, Bb,
            xc, yc, Rc, Gc, Bc,
            ns
          );
      }
  }
  
void pswr_shade_triangle_subd
  ( PSStream *ps,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  )
  {
    assert(ns > 0);
    int nt = ns + 1; /* Number of intervals on each side. */
    
    int ib, ic;
    /* Scan the chip rows: */
    for (ic = 0; ic < nt; ic++)
      { /* Rel. contrib of input corner {c} to chip corners in rows {ic} and {ic+1}: */
        double uc = ((double)ic + 0.0)/((double)nt);
        double vc = ((double)ic + 1.0)/((double)nt);
        
        /* Rel contrib of input corner {c} to colors of chip below those corners: */
        double sc = ((double)ic + 1.0/3.0)/((double)nt); /* Up chips */
        double tc = ((double)ic + 2.0/3.0)/((double)nt); /* Down chips */
        
        for (ib = 0; ib < nt-ic; ib++)
          { 
            /* Rel. contrib of input corner {b} to chip corners in cols {ib} and {ib+1}: */
            double ub = ((double)ib + 0.0)/((double)nt);
            double vb = ((double)ib + 1.0)/((double)nt);
            
            /* Rel contrib of input corner {b} to colors of chips below those corners: */
            double sb = ((double)ib + 1.0/3.0)/((double)nt); /* Up chips */
            double tb = ((double)ib + 2.0/3.0)/((double)nt); /* Down chips */
        
            /* Corners and colors of up chip: */
            double xsa = (1 - ub - uc)*xa + ub*xb + uc*xc;
            double ysa = (1 - ub - uc)*ya + ub*yb + uc*yc;
                                     
            double xsb = (1 - vb - uc)*xa + vb*xb + uc*xc;
            double ysb = (1 - vb - uc)*ya + vb*yb + uc*yc;
                                     
            double xsc = (1 - ub - vc)*xa + ub*xb + vc*xc;
            double ysc = (1 - ub - vc)*ya + ub*yb + vc*yc;
            
            double Rs = (1 - sb - sc)*Ra + sb*Rb + sc*Rc;
            double Gs = (1 - sb - sc)*Ga + sb*Gb + sc*Gc;
            double Bs = (1 - sb - sc)*Ba + sb*Bb + sc*Bc;
            
            /* Paint up chip: */
            pswr_set_fill_color(ps, Rs, Gs, Bs);
            pswr_triangle(ps, xsa, ysa, xsb, ysb, xsc, ysc, TRUE, FALSE);
            
            if (ib + ic + 2 <= nt)
              { /* Corners and colors of down chip: */
                double xta = (1 - vb - vc)*xa + vb*xb + vc*xc;
                double yta = (1 - vb - vc)*ya + vb*yb + vc*yc;

                double xtb = xsc;
                double ytb = ysc;

                double xtc = xsb;
                double ytc = ysb;

                double Rt = (1 - tb - tc)*Ra + tb*Rb + tc*Rc;
                double Gt = (1 - tb - tc)*Ga + tb*Gb + tc*Gc;
                double Bt = (1 - tb - tc)*Ba + tb*Bb + tc*Bc;

                /* Paint down chip: */
                pswr_set_fill_color(ps, Rt, Gt, Bt);
                pswr_triangle(ps, xta, yta, xtb, ytb, xtc, ytc, TRUE, FALSE);
              }
          }
      }

  }

void pswr_shade_triangle_gouraud
  ( PSStream *ps,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc
  )
  { affirm(pswr_within_canvas(ps), "bad state");
    if ((Ra < 0.0) || (Rb < 0.0) || (Rc < 0.0)) { return; }
    double psxa = ps->hMin + ps->xScale * (xa - ps->xMin);
    double psya = ps->vMin + ps->yScale * (ya - ps->yMin);
    double psxb = ps->hMin + ps->xScale * (xb - ps->xMin);
    double psyb = ps->vMin + ps->yScale * (yb - ps->yMin);
    double psxc = ps->hMin + ps->xScale * (xc - ps->xMin);
    double psyc = ps->vMin + ps->yScale * (yc - ps->yMin);
    if (pswr_triangle_is_invisible(ps, psxa, psya, psxb, psyb, psxc, psyc))
      { return; }
    FILE *file = ps->file;
    fprintf(file, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psxa, psya, Ra, Ga, Ba);
    fprintf(file, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psxb, psyb, Rb, Gb, Bb);
    fprintf(file, "%6.1f %6.1f  %5.3f %5.3f %5.3f",   psxc, psyc, Rc, Gc, Bc);
    fprintf(file, " gstrif\n");
    fflush(file);
  }
