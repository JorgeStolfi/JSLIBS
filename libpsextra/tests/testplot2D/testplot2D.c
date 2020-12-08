#define PROG_NAME "testplot2D"
#define PROG_DESC "test of {pswr_plot_2D.h}"
#define PROG_VERS "1.0"

/* Last edited on 2011-06-06 17:51:16 by stolfi */

#define testplot_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>

#include <pswr.h>
#include <pswr_iso.h>
#include <pswr_color.h>
#include <pswr_plot_2D.h>

#define DDIM (2) 
  /* Dimension of domain. */

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */

#define FIG_XSZ (12.0)
#define FIG_YSZ (12.0)
  /* Subfigure dimensions in client units. */

#define MIN_EPS_HMG (4.0)
#define MIN_EPS_VMG (4.0)
  /* Minimum margins for EPS figures (pt). */

#define MIN_PS_HMG (72.0)
#define MIN_PS_VMG (72.0)
  /* Minimum margins for PS documents (pt). */

int main (int argc, char **argv);
void SetPlotWindow(PSStream *ps, double xsz, double ysz, int ix, int iy, int nx, int ny);
void DoTests(int nx, int ny, bool_t epsformat);
void DoPaintings(PSStream *ps, int nx, int ny, bool_t quad);
void PlotFunc2D
  ( PSStream *ps,
    int ix,
    int iy,
    int nx, 
    int ny,
    bool_t messy,
    int nf,
    bool_t isolines,
    bool_t bands,
    bool_t quad,
    int ns
  );
  /* Draws a function plot on file {ps}.  The page
    is divided conceptually into an array of {nx × ny} sub-figures.
    The plot is done within the sub-figure in column
    {ix} and row {iy} (counting from 0 at bottom left).
    
    If {quad} is TRUE, tests the quadrangular region plots.
    If {quad} is FALSE, tests the triangular region plots.
    
    The function will compute {nf} values (must be 2 or more).
    It is plotted with the given{isolines} and {bands} prameters,
    and {ns} steps along each side. */ 

pswr_plot_2D_style_t *BuildStyle
  ( int nf,
    bool_t isolines,
    bool_t bands
  );
  /* Builds a style parameter for {}. */

int main (int argc, char **argv)
  { DoTests(4,5,FALSE);
    DoTests(4,5,TRUE);
    return 0;
  }

void DoTests(int nx, int ny, bool_t epsformat)
  { 
    PSStream *ps = NULL;
    if (epsformat)
      { /* Choose scale: */
        double hscale = 6.0*72/(nx*FIG_XSZ);
        double vscale = 6.0*72/(ny*FIG_XSZ);
        double scale = fmin(hscale, vscale);
        /* EPS figure size proportional to total image size. */
        double eps_hsz = ceil(nx*scale*FIG_XSZ + 2*MIN_EPS_HMG);
        double eps_vsz = ceil(ny*scale*FIG_YSZ + 2*MIN_EPS_VMG);
        ps = pswr_new_stream(OUT_PREFIX, NULL, TRUE, "doc", NULL, FALSE, eps_hsz, eps_vsz);
      }
    else
      { ps = pswr_new_stream(OUT_PREFIX, NULL, FALSE, "doc", "letter", FALSE, 0, 0); }
    DoPaintings(ps, nx, ny, FALSE);
    DoPaintings(ps, nx, ny, TRUE);
    pswr_close_stream(ps);
  }

double FA(double x, double y);
double FB(double x, double y);

void DoPaintings(PSStream *ps, int nx, int ny, bool_t quad)
  { 
    int nsBig = 25; /* Number of plot steps for messy functions. */
    int nsSma = 3;  /* Number of plot steps for smooth functions. */
    
    pswr_new_canvas(ps, (quad ? "qua" : "tri"));

    bool_t messy, isolines, bands;
    
    messy = TRUE;
    isolines = FALSE;
    bands = FALSE;
    PlotFunc2D(ps, 0, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 0, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 0, 2, nx, ny, messy, 4, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 0, 3, nx, ny, messy, 5, isolines, bands, quad, nsBig);
                                                  
    messy = TRUE;
    isolines = TRUE;
    bands = FALSE;
    PlotFunc2D(ps, 1, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 1, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 1, 2, nx, ny, messy, 4, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 1, 3, nx, ny, messy, 5, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 1, 4, nx, ny, messy, 6, isolines, bands, quad, nsBig);
                                                  
    messy = TRUE;
    isolines = TRUE;
    bands = TRUE;
    PlotFunc2D(ps, 2, 0, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 2, 1, nx, ny, messy, 3, isolines, bands, quad, nsBig);
    
    messy = TRUE;
    isolines = FALSE;
    bands = TRUE;
    PlotFunc2D(ps, 2, 3, nx, ny, messy, 2, isolines, bands, quad, nsBig);
    PlotFunc2D(ps, 2, 4, nx, ny, messy, 3, isolines, bands, quad, nsBig);
                                                  
    messy = FALSE;
    isolines = TRUE;
    bands = FALSE;
    PlotFunc2D(ps, 3, 0, nx, ny, messy, 2, isolines, bands, quad, nsSma);
    PlotFunc2D(ps, 3, 1, nx, ny, messy, 3, isolines, bands, quad, nsSma);
    PlotFunc2D(ps, 3, 2, nx, ny, messy, 4, isolines, bands, quad, nsSma);
    PlotFunc2D(ps, 3, 3, nx, ny, messy, 5, isolines, bands, quad, nsSma);
    PlotFunc2D(ps, 3, 4, nx, ny, messy, 6, isolines, bands, quad, nsSma);
  }

void PlotFunc2D
  ( PSStream *ps,
    int ix,
    int iy,
    int nx, 
    int ny,
    bool_t messy,
    int nf,
    bool_t isolines,
    bool_t bands,
    bool_t quad,
    int ns
  )
  { assert((ix >= 0) && (ix < nx));
    assert((iy >= 0) && (iy < ny));
    SetPlotWindow(ps, FIG_XSZ, FIG_YSZ, ix, iy, nx, ny);
    
    double wx = FIG_XSZ;
    double wy = FIG_YSZ;
    
    /* Usable area {[0 _ wx]×[0 _ wy]} */
    
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    /* int ticlo = 0; */
    /* int tichi = (int)floor(FIG_XY_SZ); */
    /* pswr_tics(ps, HOR, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    /* pswr_tics(ps, VER, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    pswr_rectangle(ps, -0.25, +0.25+wx, -0.25, +0.25+wy, FALSE, TRUE);
    
    char *cmt = NULL;
    asprintf(&cmt, "Testing with nf = %d, isolines = %c, ns = %d", nf, "FT"[isolines], ns);
    pswr_comment(ps, cmt);
    free(cmt);
        
    /* Choose domain to plot: */
    interval_t B[2]; /* For a quadrangle plot. */
    B[0] = (interval_t){{ 0.1, 0.9 }};
    B[1] = (interval_t){{ 0.3, 0.8 }};
    double xa[2], xb[2], xc[2]; /* For triangle plot. */
    xa[0] = LO(B[0]); xa[1] = LO(B[1]);
    xb[0] = HI(B[0]); xb[1] = (LO(B[1]) + HI(B[1]))/2;
    xc[0] = LO(B[0]); xc[1] = HI(B[1]);
    
    auto void mess(double x[], int nx, double f[], int nf);
      /* A messy black-box function to plot. */
    
    void mess(double x[], int nx, double f[], int nf)
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
      { bf = (double *)notnull(malloc(ipow(2,DDIM)*nf*sizeof(double)), "no mem");
        fprintf(stderr, "bilinear corner values:\n");
        fprintf(stderr, "\n");
        srandom(4615);
        int k, j;
        for (k = 0; k < 4; k++)
          { int i0 = k % 2;
            int i1 = k / 2;
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
        int j;
        for (j = 0; j < nf; j++) { f[j] = (1 - s)*bf0[j] + s*bf1[j]; }
      }
    
    auto void bill(double x[], int nx, double f[], int nf);
      /* A billinear function to plot, interpolates the values in {bf}. */
    
    void bill(double x[], int nx, double f[], int nf)
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

    /* Gather the parameters for {pswr_plot_2D_quad}: */
    pswr_plot_2D_func_t *F = (messy ? &mess : &bill);
    pswr_set_fill_color(ps, 0.950, 0.400, 0.000);
    pswr_set_pen(ps, 0.000, 0.000, 0.333,  0.10,  0.0, 0.0);
    pswr_plot_2D_style_t *st = BuildStyle(nf, isolines, bands);
    bool_t fill = TRUE;
    bool_t draw = (nf == 2) || isolines;
    if (quad)
      { pswr_plot_2D_quad(ps, nf, F, B, ns,ns, st, fill, draw);
        if (draw) 
          { pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
            pswr_plot_2D_quad_outline(ps, nf, F, B, ns,ns);
          }
      }
    else
      { pswr_plot_2D_tri(ps, nf, F, xa,xb,xc, ns, st, fill, draw);
        if (draw) 
          { pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
            pswr_plot_2D_tri_outline(ps, nf, F, xa,xb,xc, ns);
          }
      }
  }
  
pswr_plot_2D_style_t *BuildStyle
  ( int nf,
    bool_t isolines,
    bool_t bands
  )
  {
    pswr_plot_2D_style_t *st = (pswr_plot_2D_style_t *)notnull(malloc(sizeof(pswr_plot_2D_style_t)), "no mem");
    
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
    st->kMin = pswr_inf_isoline(st->vStart, st->vStep, vMin - eps);
    st->kMax = pswr_sup_isoline(st->vStart, st->vStep, vMax + eps);

    /* Build the color band tables: */
    int N = 0;
    st->Rtb = NULL;
    st->Gtb = NULL;
    st->Btb = NULL;
    if (isolines || bands)
      { pswr_make_color_table(
          st->vStart, st->vStep, st->kMin, st->kMax,
          0.000, 0.500, 1.000,
          1.000, 1.000, 1.000,
          1.000, 0.167, 0.000,
          &N, &(st->Rtb), &(st->Gtb), &(st->Btb));
      }
      
    return st;
  }      

void SetPlotWindow(PSStream *ps, double xsz, double ysz, int ix, int iy, int nx, int ny)
  {
    /* Compute minimum margins in pt: */
    double min_hmg = (ps->eps ? MIN_EPS_HMG : MIN_PS_HMG);
    double min_vmg = (ps->eps ? MIN_EPS_VMG : MIN_PS_VMG);
    /* Get the total page/figure size in pt: */
    double tot_hsz, tot_vsz;
    pswr_get_canvas_size(ps, &tot_hsz, &tot_vsz);
    /* Compute max width and height available for all subfigures in pt: */
    double max_hsz = tot_hsz - 2*min_hmg;
    double max_vsz = tot_vsz - 2*min_vmg;
    /* Compute scale factor (pt per client unit): */
    double hscale = max_hsz/(nx*xsz);
    double vscale = max_vsz/(ny*ysz);
    double scale = fmin(hscale, vscale);
    /* Subfigure dimensions in points: */
    double fig_hsz = scale*xsz;
    double fig_vsz = scale*ysz;
    /* Recompute true margins in points: */
    double hmg = (ps->hCanvasSize - nx*fig_hsz)/2;
    double vmg = (ps->vCanvasSize - ny*fig_vsz)/2;
    /* Compute subfigure coords in canvas: */
    double fig_hlo = hmg + ix*fig_hsz;
    double fig_hhi = fig_hlo + fig_hsz;
    double fig_vlo = vmg + iy*fig_vsz;
    double fig_vhi = fig_vlo + fig_vsz;
    /* Set the window: */
    pswr_set_window(ps, 0, xsz, 0, ysz, fig_hlo, fig_hhi, fig_vlo, fig_vhi);
  }
