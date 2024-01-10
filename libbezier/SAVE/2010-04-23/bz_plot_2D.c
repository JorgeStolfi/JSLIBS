/* See {bz_plot_2D.h}. */
/* Last edited on 2009-08-30 18:01:59 by stolfi */

#define bz_plot_2D_C_COPYRIGHT "Copyright © 2007 by the State University of Campinas (UNICAMP)."

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <rn.h>
#include <pswr.h>
#include <pswr_color.h>
#include <pswr_iso.h>
#include <pswr_shade_quad.h>

#include <bz_patch.h>

#include <bz_plot_2D.h>
     
void bz_plot_2D_patch_and_func
  ( PSStream *ps, 
    bz_patch_t *bz,
    bz_patch_rdim_t nf,
    pswr_plot_2D_func_t *func,
    int minRank,
    int maxRank,
    double tol,
    interval_t B[],
    int nx,
    int ny,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  )
  { 
    demand(bz->m == 2, "wrong domain dimension");
    demand(bz->n >= 2, "bad range dimension");
    demand(bz->n <= bz_patch_MAX_RDIM, "bad range dimension");
    int BSTEPS = 10;
    double Rf = 1.0, Gf = 0.9, Bf = 0.6;
    int nfaces = ipow(3, bz->m);
    int fix;
    fprintf(stderr, "+bz_plot_2D dim = %d\n", bz->m);
    /* Plot faces */
    for (fix = 0; fix < nfaces; fix++)
      { if (box_face_dimension(fix, bz->m) == 2)
        { bz_patch_t *t = bz_plot_2D_extract_face(bz, fix);
          bz_plot_2D_fill_face(ps, t, BSTEPS, Rf, Gf, Bf);
          free(t->c); free(t);
        }
      }
    /* Plot edges */
    pswr_set_pen(ps, 0.0, 0.0, 0.0,   0.15, 0.0, 0.0);
    for (fix = 0; fix < nfaces; fix++)
      { if (box_face_dimension(fix, bz->m) == 1)
        { bz_patch_t *t = bz_plot_2D_extract_face(bz, fix);
          bz_plot_2D_edge(ps, t, BSTEPS);
          free(t->c); free(t);
        }
      }
    /* Plot vertices */
    pswr_set_pen(ps, 0.0, 0.0, 0.0,   0.15, 0.0, 0.0);
    for (fix = 0; fix < nfaces; fix++)
      if (box_face_dimension(fix, bz->m) == 0) 
        { bz_patch_t *t = bz_plot_2D_extract_face(bz, fix);
          bz_plot_2D_vertex(ps, t, BSTEPS);
          free(t->c); free(t);
        }
    fprintf(stderr, "-bz_plot_2D\n");
  }

PSStream *bz_plot_2D_init
  ( char *prefix,
    bool_t epsformat,
    char *paperSize,
    char *caption,
    interval_t bbox[]
  )
  { double hpt, vpt;
    if (epsformat)
      { hpt = 8.0*72.0; vpt = hpt * 4.0/3.0; }
    else
      { pswr_get_paper_dimensions(paperSize, &hpt, &vpt); }
    PSStream *ps = pswr_new_stream(prefix, NULL, epsformat, "doc", paperSize, hpt+6, vpt+8);
    pswr_set_canvas_layout(ps, 3.0, 4.0, TRUE, 0.0, 0.0, (caption == NULL ? 0 : 1), 1,1);
    double dx = 0.1*(HI(bbox[0]) - LO(bbox[0]))/2;
    double dy = 0.1*(HI(bbox[1]) - LO(bbox[1]))/2;
    pswr_new_picture
      ( ps, 
        LO(bbox[0])-dx, HI(bbox[0])+dx, 
        LO(bbox[1])-dy, HI(bbox[1])+dy
      );
    pswr_set_pen(ps, 0,0,0, 0.15, 0,0);
    if (caption != NULL)
      { pswr_add_caption(ps, caption, 0.0); }
    pswr_set_pen(ps, 1,0,0, 0.25, 0,0);
    return ps;
  }

void bz_plot_2D_finish ( PSStream *ps )
  { pswr_close_stream(ps); }

void bz_plot_2D_patch_and_func
  ( PSStream *ps, 
    bz_patch_t *bz,
    bz_patch_rdim_t nf,
    bz_plot_2D_func_t *func,
    int nx,
    bool_t isolines
  )
  {
    /* The domain dimension must be 2: */
    bs_ddim_t ddim = 2;
    
    /* Get and check {nb}: */
    int nb; /* Number of values from {bz}. */
    if (bz == NULL)
      { nb = 0; }
    else
      { demand(bz->m == ddim, "Bezier patch has wrong domain dimension");
        nb = bz->n;
      }
    
    /* Check {nf}: */
    if (func == NULL) 
      { demand(nf == 0, "null function should have zero range dimension"); }
    else
      { demand(nf >= 0, "function has bad range dimension"); }
    
    /* Check the total range dimension {nb + nf} of the concatenated function {T}: */
    demand(nb + nf > ddim, "total range dimension too low");
    
    /* Number of function values {nv}: */
    int nv = nb + nf - ddim;
    if (isolines)
      { demand(nv == 1, "isolines/bands require exactly one function value"); }
    else
      { demand(nv <= 3, "too many function values for shading"); }

    auto void concfunc(int nx, double [x], int nt, double t[]);
    /* Evaluates the concatenated function {T(x)}, puts the result in {t[0..nt-1]} */
    
    void concfunc(int nx, double [x], int nt, double t[])
      {
        assert(nt == nb + nf);
        assert(nx == ddim);
        if ((bz != NULL) && (nb != 0))
          { /* Evaluate {bz(x)}, save in {t[0..nb-1]}: */
            double *b = &(t[0]);
            bz_patch_eval(bz, x, b);
          }
        if ((func != NULL) && (nf != 0))
          { /* Evaluate {func(x)}, save in {t[nb..nb + nf-1]}: */
            double *f = &(t[nb]);
            func(ddim, x, nf, f);
          }
      }

    /* Choose the nominal range for function values: */
    double fMax = 1.0;
    double fMin = -fMax;

    /* Choose the isoline spacing: */
    double fStep = 0.2;
    double fSync = fStep/2;

    /* Compute the number of isolines: */
    int kMin = pswr_inf_isoline(fSync, fStep, -fMax);
    int kMax = pswr_sup_isoline(fSync, fStep, +fMax);

    /* Build the color band tables: */
    int N = 0;
    double *Rtb = NULL, *Gtb = NULL, *Btb = NULL;
    if (isolines)
      { pswr_make_color_table(
          fSync, fStep, kMin, kMax,
          0.000, 0.500, 1.000,
          1.000, 1.000, 1.000,
          1.000, 0.167, 0.000,
          &N, &Rtb, &Gtb, &Btb);
      }
   
    int ix0, ix1;
    inerval_t B[DDIM];
    for (ix0 = 0; ix0 < nx; ix0++) 
      { LO(B[0]) = ((double)ix0)/((double)nx);
        HI(B[0]) = ((double)ix0+1)/((double)nx);
        for (ix1 = 0; ix1 < nx; ix1++) 
          { LO(B[1]) = ((double)ix1)/((double)nx);
            HI(B[1]) = ((double)ix1+1)/((double)nx);
            pswr_plot_2D_face(ps, nt, concfunc, B, nx, isolines, Rtb, Gtb, Btb);
          }
      }
      
    if (isolines)
      { free(Rtb); free(Gtb); free(Btb); }
  }
  
     
/* INTERNAL DEFINITIONS */

#define MAX_VALUES (bz_plot_2D_func_MAX_VALUES)
  /* Shorter name for the max function values (besides shape coords). */

#define MIN_RDIM (bz_plot_2D_func_DDIM)
  /* Minimum dimension for the range of {F}. */

#define MAX_RDIM (bz_plot_2D_func_DDIM + bz_plot_2D_func_MAX_VALUES)
  /* Maximum dimension for the range of {F}. */

#define bz_plot_2D_func_NF (9)
  /* Number of faces in a cell of the domain grid ({3^DDIM}). */

#define bz_plot_2D_func_TOL_MM (0.5)
  /* Flatness tolerance (mm). */

#define bz_plot_2D_func_DEBUG_SHAPE FALSE
  /* TRUE to debug {make_shape} instead of {dg_enum_faces}. */

/* IMPLEMENTATIONS */

void bz_plot_2D_fill_face_shade
  ( PSStream *ps,
    double c00[], 
    double c01[],
    double c10[], 
    double c11[],
    bz_patch_rdim_t n,
    int plotDepth,
    double fMin,
    double fMax
  )
  {
    int i0, i1, j;
    int steps = ipow(2, plotDepth);
    double p[bz_patch_MAX_RDIM], q[bz_patch_MAX_RDIM], r[bz_patch_MAX_RDIM], s[bz_patch_MAX_RDIM];
    for (i0 = 1; i0 <= steps; i0++)
      { for (i1 = 0; i1 <= steps; i1++)
          { double u0 = ((double)i0)/((double)steps);
            double u1 = ((double)i1)/((double)steps);
            double v0 = ((double)i0-1)/((double)steps);
            double v1 = ((double)i1)/((double)steps);
            for (j = 0; j < n; j++)
              { p[j] = dg_bilinear(u0,u1,c00[j],c01[j],c10[j],c11[j]);
                r[j] = dg_bilinear(v0,v1,c00[j],c01[j],c10[j],c11[j]);

              }
            if (i1 > 0) { dg_shade_quadrilateral(ps, p, q, r, s, n, fMin, fMax); }
            for (j = 0; j < n; j++) { q[j] = p[j]; s[j] = r[j]; }
          }
      }
  }

void bz_plot_2D_fill_face_bands
  ( PSStream *ps,
    double c00[], 
    double c01[],
    double c10[], 
    double c11[],
    bz_patch_rdim_t n,
    int plotDepth,
    double fStart,  /* Synchronize isolines with this level. */
    double fStep,   /* isoline spacing. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  )
  {
    demand(n == 3, "can't plot bands with more than one function value"); 
    int i0, i1, j;
    int steps = ipow(2, plotDepth);
    double p[bz_patch_MAX_RDIM], q[bz_patch_MAX_RDIM], r[bz_patch_MAX_RDIM], s[bz_patch_MAX_RDIM];
    for (i0 = 1; i0 <= steps; i0++)
      { for (i1 = 0; i1 <= steps; i1++)
          { double u0 = ((double)i0)/((double)steps);
            double u1 = ((double)i1)/((double)steps);
            double v0 = ((double)i0-1)/((double)steps);
            double v1 = ((double)i1)/((double)steps);
            for (j = 0; j < n; j++)
              { p[j] = dg_bilinear(u0,u1,c00[j],c01[j],c10[j],c11[j]);
                r[j] = dg_bilinear(v0,v1,c00[j],c01[j],c10[j],c11[j]);

              }
            if (i1 > 0) 
              { dg_bands_in_quadrilateral
                  (ps, p, q, r, s, n, fStart, fStep, kMin, kMax, R, G, B);
              }
            for (j = 0; j < n; j++) { q[j] = p[j]; s[j] = r[j]; }
          }
      }
  }

void bz_plot_2D_fill_face_isolines
  ( PSStream *ps,
    double c00[], 
    double c01[],
    double c10[], 
    double c11[],
    bz_patch_rdim_t n,
    int plotDepth,
    double fStart,  /* Synchronize isolines with this level. */
    double fStep,   /* Isoline spacing. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  )
  {
    demand(n == 3, "can't plot isolines with more than one function value"); 
    int i0, i1, j;
    int steps = ipow(2, plotDepth);
    double p[bz_patch_MAX_RDIM], q[bz_patch_MAX_RDIM], r[bz_patch_MAX_RDIM], s[bz_patch_MAX_RDIM];
    for (i0 = 1; i0 <= steps; i0++)
      { for (i1 = 0; i1 <= steps; i1++)
          { double u0 = ((double)i0)/((double)steps);
            double u1 = ((double)i1)/((double)steps);
            double v0 = ((double)i0-1)/((double)steps);
            double v1 = ((double)i1)/((double)steps);
            for (j = 0; j < n; j++)
              { p[j] = dg_bilinear(u0,u1,c00[j],c01[j],c10[j],c11[j]);
                r[j] = dg_bilinear(v0,v1,c00[j],c01[j],c10[j],c11[j]);

              }
            if (i1 > 0) 
              { dg_isolines_in_quadrilateral(ps, p, q, r, s, n, fStart, fStep, kMin, kMax); }
            for (j = 0; j < n; j++) { q[j] = p[j]; s[j] = r[j]; }
          }
      }
  }

void dg_shade_quadrilateral
  ( PSStream *ps,
    double p00[],   /* Coords and values at corner [0,0]. */
    double p01[],   /* Coords and values at corner [0,1]. */
    double p10[],   /* Coords and values at corner [1,0]. */
    double p11[],   /* Coords and values at corner [1,1]. */
    bz_patch_rdim_t n,      /* Dimension of range space. */
    double fMin,    /* Minimum function value. */
    double fMax     /* Maximum function value. */
  )
  {
    double ctr[bz_patch_MAX_RDIM];
    int j;
    for (j = 0; j < n; j++) 
      { ctr[j] = (p00[j]+p01[j]+p10[j]+p11[j])/4.0; }
    /* Paint triangles: */
    dg_shade_triangle(ps, p00, p01, ctr, n, fMin, fMax);
    dg_shade_triangle(ps, p01, p11, ctr, n, fMin, fMax);
    dg_shade_triangle(ps, p11, p10, ctr, n, fMin, fMax);
    dg_shade_triangle(ps, p10, p00, ctr, n, fMin, fMax);
  }

void dg_bands_in_quadrilateral
  ( PSStream *ps,
    double p00[],   /* Coords and values at corner [0,0]. */
    double p01[],   /* Coords and values at corner [0,1]. */
    double p10[],   /* Coords and values at corner [1,0]. */
    double p11[],   /* Coords and values at corner [1,1]. */
    bz_patch_rdim_t n,      /* Dim of point vectors (incl. function values). */
    double fStart,  /* Synchronize isolines with this level. */
    double fStep,   /* Isoline spacing. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  )
  {
    demand(n == 3, "can't plot bands with more than one function value");
    pswr_bands_in_quadrilateral
      ( ps,
        p00[0], p00[1], p00[2],
        p01[0], p01[1], p01[2],
        p10[0], p10[1], p10[2],
        p11[0], p11[1], p11[2],
        fStart, fStep, kMin, kMax,
        R, G, B
      );
  }

void dg_isolines_in_quadrilateral
  ( PSStream *ps,
    double p00[],   /* Coords and values at corner [0,0]. */
    double p01[],   /* Coords and values at corner [0,1]. */
    double p10[],   /* Coords and values at corner [1,0]. */
    double p11[],   /* Coords and values at corner [1,1]. */
    bz_patch_rdim_t n,      /* Dimension of range space. */
    double fStart,   /* Synchronize isolines with this level. */
    double fStep,   /* Isoline spacing. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  )
  {
    demand(n == 3, "can't plot isolines with more than one function value");
    pswr_isolines_in_quadrilateral
      ( ps,
        p00[0], p00[1], p00[2],
        p01[0], p01[1], p01[2],
        p10[0], p10[1], p10[2],
        p11[0], p11[1], p11[2],
        fStart, fStep, kMin, kMax
      );
  }

void dg_shade_triangle
  ( PSStream *ps,
    double a[],     /* First vertex. */
    double b[],     /* Second vertex. */
    double c[],     /* Third vertex. */
    bz_patch_rdim_t n,      /* Dimension of range space. */
    double fMin,    /* Minimum function value. */
    double fMax     /* Maximum function value. */
  )
  { double R, G, B;
    double f[bz_patch_MAX_RDIM];
    int j;
    affirm (n >= 2, "invalid point coordinates");
    /* Compute mean function values normalized to [-1 _ +1]: */
    for (j = 0; j < n-2; j++)
      { /* Compute mean function value: */
        f[j] = (a[j+2]+b[j+2]+c[j+2])/3.0;
        /* Clip to {[fMin _ fMax]}: */
        if (f[j] < fMin) { f[j] = fMin; }
        if (f[j] > fMax) { f[j] = fMax; }
        /* Normalize to [-1 _ +1]: */
        if (f[j] > 0.0) 
          { f[j] /= fMax; }
        else if (f[j] < 0.0) 
          { f[j] /= -fMin; }
      }
    /* Compute color {R,G,B} for white-background plotting: */
    double Rs = 1.0, Gs = 0.342, Bs = 0.000;
    double Y0 = 0.95;
    double Y1 = 0.30;
    if (n == 2)
      { R = 1.0; G = 0.9; B = 0.8; }
    else if (n == 3)
      { pswr_color_scale_1(f[0], Rs, Gs, Bs, Y0, &R, &G, &B); }
    else if (n == 4)
      { pswr_color_scale_2(f[0], f[1], Rs, Gs, Bs, Y0, &R, &G, &B); }
    else if (n >= 5)
      { pswr_color_scale_3(f[0], f[1], f[2], Y0, Y1, &R, &G, &B); }
    pswr_set_fill_color(ps, R,G,B);
    pswr_triangle(ps, a[0], a[1], b[0], b[1], c[0], c[1], TRUE, FALSE);
  }

void bz_plot_2D_vertex(PSStream *ps, bz_patch_t *b, int steps)
  {
    fprintf(stderr, "+bz_plot_2D_vertex dim = %d\n", b->m);
    affirm(b->m == 0, "invalid vertex patch");
    pswr_set_fill_color(ps, 0.0,0.0,0.0);
    pswr_dot(ps, b->c[0], b->c[1], 0.25, TRUE, FALSE);
    fprintf(stderr, "-bz_plot_2D_vertex\n");
  }

void bz_plot_2D_edge(PSStream *ps, bz_patch_t *b, int steps)
  {
    fprintf(stderr, "+bz_plot_2D_edge dim = %d\n", b->m);
    affirm(b->m == 1, "invalid edge patch");
    double p[2], q[2];
    int i;
    fprintf(stderr, "shape =\n");
    bz_patch_print(stderr, b, "%6.2f");
    for (i = 0; i <= steps; i++)
      { double ui = ((double)i)/((double)steps);
        bz_patch_eval(b, &ui, p);
        if (i > 0) { pswr_segment(ps, q[0], q[1], p[0], p[1]); }
        q[0] = p[0]; q[1] = p[1];
      }
    fprintf(stderr, "-bz_plot_2D_edge\n");
  }

void bz_plot_2D_fill_face(PSStream *ps, bz_patch_t *b, int steps, double Rf, double Gf, double Bf)
  {
    int i, j;
    fprintf(stderr, "+bz_plot_2D_fill_face dim = %d\n", b->m);
    affirm(b->m == 2, "invalid face patch");
    fprintf(stderr, "shape =\n");
    bz_patch_print(stderr, b, "%6.2f");
    double u[2], v[2], p[MAX_RDIM], q[MAX_RDIM], r[MAX_RDIM], s[MAX_RDIM];
    q[0] =  q[1] = s[1] = s[0] = 0; /* To pacify the compiler. */
    for (i = 1; i <= steps; i++)
      { for (j = 0; j <= steps; j++)
          { u[0] = ((double)i)/((double)steps);
            u[1] = ((double)j)/((double)steps);
            bz_patch_eval(b, u, p);
            v[0] = ((double)i-1)/((double)steps);
            v[1] = ((double)j)/((double)steps);
            bz_patch_eval(b, v, r);
            if (j > 0) 
              { double xm = (p[0]+q[0]+r[0]+s[0])/4.0;
                double ym = (p[1]+q[1]+r[1]+s[1])/4.0;
                pswr_set_fill_color(ps, Rf,Gf,Bf);
                pswr_triangle(ps, p[0], p[1], q[0], q[1], xm, ym, TRUE, FALSE);
                pswr_triangle(ps, p[0], p[1], r[0], r[1], xm, ym, TRUE, FALSE);
                pswr_triangle(ps, s[0], s[1], q[0], q[1], xm, ym, TRUE, FALSE);
                pswr_triangle(ps, s[0], s[1], r[0], r[1], xm, ym, TRUE, FALSE);
              }
            q[0] = p[0]; q[1] = p[1];
            s[0] = r[0]; s[1] = r[1];
          }
      }
    fprintf(stderr, "-bz_plot_2D_fill_face\n");
  }

bz_patch_t *bz_plot_2D_extract_face(bz_patch_t *b, box_face_index_t fix)
  {
    mbz_dim_t bm = b->m, tm = 0;
    box_signed_dir_t loc[DDIM];
    box_face_signature(fix, bm, loc);
    bz_patch_t *t = (bz_patch_t *)notnull(malloc(sizeof(bz_patch_t)), "no mem");
    bz_degree_t tg[DDIM];
    int i;
    /* Extract degree sequence: */
    for(i = 0; i < bm; i++)
      { if (loc[i] == box_ZER) { tg[tm] = b->g[i]; tm++; } }
    (*t) = bz_patch_new(tm, b->n, tg);
    bz_patch_get_face(b, loc, t);
    return t;
  }

