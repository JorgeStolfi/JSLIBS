#define PROG_NAME "testplot2D"
#define PROG_DESC "test of {epswr_plot_2D.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-27 20:14:57 by stolfi */

#define testplot_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jsfile.h>
#include <rn.h>

#include <epswr.h>
#include <epswr_iso.h>
#include <epswr_color.h>
#include <epswr_plot_2D.h>

#define DDIM (2) 
  /* Dimension of domain. */

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */

#define XSZ (12.0)
#define YSZ (12.0)
  /* Subfigure dimensions in client units. */

#define HMG (4.0)
#define VMG (4.0)
  /* Minimum margins for EPS figures (pt). */

int32_t main (int32_t argc, char **argv);
void SetPlotWindow(epswr_figure_t *eps, double scale, int32_t ix, int32_t iy, int32_t nx, int32_t ny);
void DoTests(int32_t nx, int32_t ny, bool_t quad);
void DoPaintings(epswr_figure_t *eps, double scale, int32_t nx, int32_t ny, bool_t quad);
void PlotFunc2D
  ( epswr_figure_t *eps, 
    double scale,
    int32_t ix,
    int32_t iy,
    int32_t nx,
    int32_t ny,
    bool_t messy,
    int32_t nf,
    bool_t isolines,
    bool_t bands,
    bool_t quad,
    int32_t ns
  );
  /* Draws a function plot in slot {ix,iy} of {nx} by {ny}
    slot array.
    
    If {quad} is TRUE, tests the quadrangular region plots.
    If {quad} is FALSE, tests the triangular region plots.
    
    The function will compute {nf} values (must be 2 or more).
    It is plotted with the given{isolines} and {bands} prameters,
    and {ns} steps along each side. */ 

epswr_plot_2D_style_t *BuildStyle
  ( int32_t nf,
    bool_t isolines,
    bool_t bands
  );
  /* Builds a style parameter for {}. */

int32_t main (int32_t argc, char **argv)
  { DoTests(4,5, FALSE);
    DoTests(4,5, TRUE);
    return 0;
  }

void DoTests(int32_t nx, int32_t ny, bool_t quad)
  { 
    epswr_figure_t *eps = NULL;
    /* Choose scale: */
    double hscale = 8.0*72/(nx*XSZ);
    double vscale = 6.0*72/(ny*YSZ);
    double scale = fmin(hscale, vscale);
    /* EPS figure size proportional to total image size. */
    double eps_hsz = nx*(scale*XSZ + HMG) - HMG;
    double eps_vsz = ny*(scale*YSZ + VMG) - VMG;
    char *fname = NULL;
    asprintf(&fname, "out/plot2D_%s.eps", (quad ? "qua" : "tri"));
    FILE *wr = open_write(fname, TRUE);
    eps = epswr_new_figure(wr, eps_hsz, eps_vsz, HMG, HMG, VMG, VMG, FALSE);
    DoPaintings(eps, scale, nx, ny, quad);
    epswr_end_figure(eps);
  }

double FA(double x, double y);
double FB(double x, double y);

void DoPaintings(epswr_figure_t *eps, double scale, int32_t nx, int32_t ny, bool_t quad)
  { 
    int32_t nsBig = 25; /* Number of plot steps for messy functions. */
    int32_t nsSma = 3;  /* Number of plot steps for smooth functions. */
    
    bool_t messy, isolines, bands;
    
    messy = TRUE;
    isolines = FALSE;
    bands = FALSE;
    PlotFunc2D(eps, scale, 0, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 0, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 0, 2, nx, ny, messy, 4, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 0, 3, nx, ny, messy, 5, isolines, bands, quad, nsBig);
                                                  
    messy = TRUE;
    isolines = TRUE;
    bands = FALSE;
    PlotFunc2D(eps, scale, 1, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 1, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 1, 2, nx, ny, messy, 4, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 1, 3, nx, ny, messy, 5, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 1, 4, nx, ny, messy, 6, isolines, bands, quad, nsBig);
                                                  
    messy = TRUE;
    isolines = TRUE;
    bands = TRUE;
    PlotFunc2D(eps, scale, 2, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 2, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    
    messy = TRUE;
    isolines = FALSE;
    bands = TRUE;
    PlotFunc2D(eps, scale, 2, 3, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(eps, scale, 2, 4, nx, ny, messy, 3, isolines, bands, quad, nsBig);
                                                  
    messy = FALSE;
    isolines = TRUE;
    bands = FALSE;
    PlotFunc2D(eps, scale, 3, 0, nx, ny, messy, 2, isolines, bands, quad, nsSma);
    PlotFunc2D(eps, scale, 3, 1, nx, ny, messy, 3, isolines, bands, quad, nsSma);
    PlotFunc2D(eps, scale, 3, 2, nx, ny, messy, 4, isolines, bands, quad, nsSma);
    PlotFunc2D(eps, scale, 3, 3, nx, ny, messy, 5, isolines, bands, quad, nsSma);
    PlotFunc2D(eps, scale, 3, 4, nx, ny, messy, 6, isolines, bands, quad, nsSma);
  }

void PlotFunc2D
  ( epswr_figure_t *eps,
    double scale,
    int32_t ix,
    int32_t iy,
    int32_t nx, 
    int32_t ny,
    bool_t messy,
    int32_t nf,
    bool_t isolines,
    bool_t bands,
    bool_t quad,
    int32_t ns
  )
  { double wx = XSZ;
    double wy = YSZ;
    /* Usable area {[0 _ wx]×[0 _ wy]} */
    
    SetPlotWindow(eps, scale, ix, iy, nx, ny);
    
    epswr_set_pen(eps, 0.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    /* int32_t ticlo = 0; */
    /* int32_t tichi = (int32_t)floor(XY_SZ); */
    /* epswr_tics(eps, HOR, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    /* epswr_tics(eps, VER, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    epswr_rectangle(eps, -0.25, +0.25+wx, -0.25, +0.25+wy, FALSE, TRUE);
    
    char *cmt = NULL;
    asprintf(&cmt, "Testing with nf = %d, isolines = %c, ns = %d", nf, "FT"[isolines], ns);
    epswr_comment(eps, cmt);
    free(cmt);
        
    /* Choose domain to plot: */
    interval_t B[2]; /* For a quadrangle plot. */
    B[0] = (interval_t){{ 0.1, 0.9 }};
    B[1] = (interval_t){{ 0.3, 0.8 }};
    double xa[2], xb[2], xc[2]; /* For triangle plot. */
    xa[0] = LO(B[0]); xa[1] = LO(B[1]);
    xb[0] = HI(B[0]); xb[1] = (LO(B[1]) + HI(B[1]))/2;
    xc[0] = LO(B[0]); xc[1] = HI(B[1]);
    
    auto void mess(double x[], int32_t nx, double f[], int32_t nf);
      /* A messy black-box function to plot. */
    
    void mess(double x[], int32_t nx, double f[], int32_t nf)
      { 
        demand(nx == DDIM, "bad {nx}");
        demand(nf >= DDIM, "bad {nf}");
        /* Position: map {x} from polar to Cartesian: */
        double tht = M_PI/2*x[0];
        double rho = fmin(wx,wy)*x[1];
        f[0] = rho*cos(tht);
        f[1] = rho*sin(tht);
        /* Color values are waves in domain: */
        double r = hypot(x[0],x[1]); 
        if (nf >= 3) { f[2] = 0.8*cos(1.5 * 2*M_PI * r); }
        if (nf >= 4) { f[3] = 0.8*cos(2.5 * 2*M_PI * x[0]); }
        if (nf >= 5) { f[4] = 0.8*cos(3.5 * 2*M_PI * x[1]); }
        if (nf >= 6) { f[5] = 0.8*cos(4.5 * 2*M_PI * atan2(x[1],x[0])/(M_PI/2)); }
        demand(nf <= 6, "{nf} too high");
      }
    
    /* Generate corner values for bilinear test function: */
    double *bf = NULL;  /* {bf[(2*i1 + i0)*nf + j]} is value {j} of corner {(i0,i1)}. */ 
    if (! messy)
      { bf = rn_alloc((int32_t)ipow(2,DDIM)*nf);
        fprintf(stderr, "bilinear corner values:\n");
        fprintf(stderr, "\n");
        srandom(4615);
        int32_t k, j;
        for (k = 0; k < 4; k++)
          { int32_t i0 = k % 2;
            int32_t i1 = k / 2;
            fprintf(stderr, "bf%d%d =", i0,i1);
            double *bfk = &(bf[k*nf]);
            bfk[0] = wx*(0.05 + 0.90*i0 - 0.2*(2*i0-1)*drandom());
            bfk[1] = wy*(0.05 + 0.90*i1 - 0.2*(2*i1-1)*drandom());
            for (j = 2; j < nf; j++) { bfk[j] = 0.8*(2*drandom() - 1); }
            for (j = 0; j < nf; j++) { fprintf(stderr, " %12.4f", bfk[j]); }
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
      }
   
    auto void interp(double z, interval_t *B, double bf0[], double bf1[], double f[]);
      /* Affine interpolation betwen {bf0[0..nf-1]} and {bf1[0..nf-1]} 
        using the position of {z} in {B}. */
        
    void interp(double z, interval_t *B, double bf0[], double bf1[], double f[])
      { double s = (z - B->end[0])/(B->end[1] - B->end[0]);
        int32_t j;
        for (j = 0; j < nf; j++) { f[j] = (1 - s)*bf0[j] + s*bf1[j]; }
      }
    
    auto void bill(double x[], int32_t nx, double f[], int32_t nf);
      /* A billinear function to plot, interpolates the values in {bf}. */
    
    void bill(double x[], int32_t nx, double f[], int32_t nf)
      { 
        demand(nx == DDIM, "bad {nx}");
        demand(nf >= DDIM, "bad {nf}");
        double *bf00 = &(bf[0*nf]);
        double *bf01 = &(bf[2*nf]);
        double *bf10 = &(bf[1*nf]);
        double *bf11 = &(bf[3*nf]);
        
        double bfx0[nf], bfx1[nf];
        interp(x[0], &(B[0]), bf00, bf10, bfx0);
        interp(x[0], &(B[0]), bf01, bf11, bfx1);
        interp(x[1], &(B[1]), bfx0, bfx1, f);
      }

    /* Gather the parameters for {epswr_plot_2D_quad}: */
    epswr_plot_2D_func_t *F = (messy ? &mess : &bill);
    epswr_set_fill_color(eps, 0.950, 0.400, 0.000);
    epswr_set_pen(eps, 0.000, 0.000, 0.333,  0.10,  0.0, 0.0);
    epswr_plot_2D_style_t *st = BuildStyle(nf, isolines, bands);
    bool_t fill = TRUE;
    bool_t draw = (nf == 2) || isolines;
    if (quad)
      { epswr_plot_2D_quad(eps, nf, F, B, ns,ns, st, fill, draw);
        if (draw) 
          { epswr_set_pen(eps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
            epswr_plot_2D_quad_outline(eps, nf, F, B, ns,ns);
          }
      }
    else
      { epswr_plot_2D_tri(eps, nf, F, xa,xb,xc, ns, st, fill, draw);
        if (draw) 
          { epswr_set_pen(eps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
            epswr_plot_2D_tri_outline(eps, nf, F, xa,xb,xc, ns);
          }
      }
  }
  
epswr_plot_2D_style_t *BuildStyle
  ( int32_t nf,
    bool_t isolines,
    bool_t bands
  )
  {
    epswr_plot_2D_style_t *st = (epswr_plot_2D_style_t *)notnull(malloc(sizeof(epswr_plot_2D_style_t)), "no mem");
    
    /* Choose the values to plot: */
    st->nc = (bands ? 0 : (isolines ? nf - 3 : nf - 2));
    st->ic = (st->nc == 0 ? -1 : (isolines ? 3 : 2));  
    st->ib = (bands && (nf > 2) ? 2 : -1);
    st->iv = (isolines && (nf > 2) ? 2 : -1);
    
    /* Choose the nominal range for function values: */
    double vMax = 1.0;
    double vMin = -vMax;

    /* Choose the isoline spacing: */
    st->vStep = 0.2;
    st->vStart = st->vStep/2;

    /* Compute the number of isolines: */
    double eps = 1.0e-13 * fmin(st->vStep, fabs(vMax - vMin));
    st->kMin = epswr_inf_isoline(st->vStart, st->vStep, vMin - eps);
    st->kMax = epswr_sup_isoline(st->vStart, st->vStep, vMax + eps);

    /* Build the color band tables: */
    int32_t N = 0;
    st->Rtb = NULL;
    st->Gtb = NULL;
    st->Btb = NULL;
    if (isolines || bands)
      { epswr_make_color_table(
          st->vStart, st->vStep, st->kMin, st->kMax,
          0.000, 0.500, 1.000,
          1.000, 1.000, 1.000,
          1.000, 0.167, 0.000,
          &N, &(st->Rtb), &(st->Gtb), &(st->Btb));
      }
      
    return st;
  }      

void SetPlotWindow(epswr_figure_t *eps, double scale, int32_t ix, int32_t iy, int32_t nx, int32_t ny)
  {
    assert((ix >= 0) && (ix < nx));
    assert((iy >= 0) && (iy < ny));

    double hsz = scale*XSZ;
    double vsz = scale*YSZ;
    double hmin = HMG + ix*(hsz + HMG), hmax = hmin + hsz;
    double vmin = VMG + iy*(vsz + VMG), vmax = vmin + hsz;
    
    /* Set the window: */
    epswr_set_window(eps, hmin, hmax, vmin, vmax, FALSE, 0, XSZ, 0, YSZ);
  }
