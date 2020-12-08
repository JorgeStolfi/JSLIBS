/* See epswr_shade_tri.h */
/* Last edited on 2009-08-25 01:07:58 by stolfi */

#define _GNU_SOURCE
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <epswr.h>
#include <epswr_def.h>
#include <epswr_dev.h>
#include <epswr_vis.h>

#include <epswr_shade_tri.h>

void epswr_shade_triangle_subd
  ( epswr_figure_t *epsf,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  );
  /* Same as epswr_shade_triangle, for {ns > 0}. */

void epswr_shade_triangle_gouraud
  ( epswr_figure_t *epsf,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc
  );
  /* Same as epswr_shade_triangle, but uses true Gouraud shading. */
  
void epswr_shade_triangle
  ( epswr_figure_t *epsf,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  )
  { if (ns < 0)
      { /* Gouraud-shaded triangle: */
        epswr_shade_triangle_gouraud
          ( epsf,
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
        epswr_set_fill_color(epsf, Rt, Gt, Bt);
        epswr_triangle(epsf, xa, ya, xb, yb, xc, yc, TRUE, FALSE);
      }
    else
      { /* Shading by subdivision: */
        epswr_shade_triangle_subd
          ( epsf,
            xa, ya, Ra, Ga, Ba,
            xb, yb, Rb, Gb, Bb,
            xc, yc, Rc, Gc, Bc,
            ns
          );
      }
  }
  
void epswr_shade_triangle_subd
  ( epswr_figure_t *epsf,
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
            epswr_set_fill_color(epsf, Rs, Gs, Bs);
            epswr_triangle(epsf, xsa, ysa, xsb, ysb, xsc, ysc, TRUE, FALSE);
            
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
                epswr_set_fill_color(epsf, Rt, Gt, Bt);
                epswr_triangle(epsf, xta, yta, xtb, ytb, xtc, ytc, TRUE, FALSE);
              }
          }
      }

  }

void epswr_shade_triangle_gouraud
  ( epswr_figure_t *epsf,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc
  )
  { if ((Ra < 0.0) || (Rb < 0.0) || (Rc < 0.0)) { return; }
    double psxa; epswr_x_to_h_coord(epsf, xa, &(psxa));
    double psya; epswr_y_to_v_coord(epsf, ya, &(psya));
    double psxb; epswr_x_to_h_coord(epsf, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(epsf, yb, &(psyb));
    double psxc; epswr_x_to_h_coord(epsf, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(epsf, yc, &(psyc));
    if (epswr_triangle_is_invisible(epsf, psxa, psya, psxb, psyb, psxc, psyc))
      { return; }
    FILE *wr = epsf->wr;
    fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psxa, psya, Ra, Ga, Ba);
    fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psxb, psyb, Rb, Gb, Bb);
    fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f",   psxc, psyc, Rc, Gc, Bc);
    fprintf(wr, " gstrif\n");
    fflush(wr);
  }
