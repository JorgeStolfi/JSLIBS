/* See {pswr_plot_2D.h}. */
/* Last edited on 2020-10-01 20:33:55 by jstolfi */

#define pswr_plot_2D_C_COPYRIGHT "Copyright © 2007 by the State University of Campinas (UNICAMP)."

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <rn.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <pswr.h>
#include <pswr_color.h>
#include <pswr_iso.h>
#include <pswr_shade_tri.h>

#include <pswr_plot_2D.h>

void pswr_plot_2D_quad
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t B[],
    int n0,
    int n1,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  )
  {
    /* The domain dimension must be 2: */
    int ddim = 2;
    
    /* We need at least 2 values for positional coordinates: */
    demand(nf >= 2, "too few function values to plot.");
    demand(F != NULL, "cannot plot a null function.");
    
    /* Number of tile corners along each scanline: */
    double f[(n0+1)*nf]; /* Value vectors at corners of previous scanline. */  
    /* The values at tile corner {k0} are {f[k0*nf+j]}, for {j=0..nf-1}. */
    
    double fp[nf], fq[nf]; /* Work vectors for {F}-vals at corners of a 4-sided tile. */
    double fmm[nf];        /* Work vector for the {F}-vals at the center of the tile. */
    double x11[ddim];      /* Upper right corner of current tile. */
    double xmm[ddim];      /* Center of current tile. */
    
    int pass; /* Pass 0 is fill, pass 1 is draw: */
    for (pass = 0; pass < 2; pass++)
      { if ((pass == 0) && (! fill)) { continue; }
        if ((pass == 1) && (! draw)) { continue; }

        bool_t fillP = (pass == 0);
        bool_t drawP = (pass == 1);

        /* Scan the tile corners: */
        int k0, k1;
        for (k1 = 0; k1 <= n1; k1++)
          { /* Now, if {k1 > 0}, {f} contains the {F}-vals for the vertices in scanline {k1-1}. */

            /* Compute domain coordinate {x11[1]} of this scanline: */
            double u1 = ((double)k1)/((double)n1);
            x11[1] = (1 - u1)*LO(B[1]) + u1*HI(B[1]);

            /* Compute domain coordinate {xmm[1]} of tile centers in this scanline: */
            double v1 = ((double)k1 - 0.5)/((double)n1);
            xmm[1] = (1 - v1)*LO(B[1]) + v1*HI(B[1]);

            /* Initialize the pointers to the top corners of current tile (will swap): */
            double *f01 = fp;
            double *f11 = fq;

            /* Initialize the pointers to the bottom corners of current tile: */
            double *f00 = NULL;
            double *f10 = NULL;

            for (k0 = 0; k0 <= n0; k0++)
              { /* Now {f} contains the {F}-values for the vertices {0..k0-1} of scanline {k1}, */
                /* and vertices {k0..n0} of scanline {k1-1}. */
                /* Also, if {k0 > 0}, */
                /*   {f10} points to a segment of {f}. */
                /*   {f01} and {f11} are scratch arrays. */
                /*   {f10} has the {F}-vals for row {k1-1}, column {k0-1}. */
                /*   {f11} has the {F}-vals of row {k1}, column {k0-1}. */

                /* We will plot the domain tile between rows {k1-1} and {k1}, cols {k0-1} and {k0}. */
                /* If {k1 == 0} the tile does not exist but we need to compute its top corners. */

                /* Swap the pointers to the top corners of the current tile: */
                { double *tmp = f01; f01 = f11; f11 = tmp; }
                /* Now {f01} is top left corner of curr tile, and {f11} is scratch area. */

                /* Update the pointers to the bottom corners of current tile: */
                if (k0 > 0) { f00 = f10; } else { f00 = NULL; }
                f10 = &(f[k0*nf]);
                /* Now {f00} and {f10} point to segments of {f[..]}; */
                /* {f00} is the {F}-vals of {k1-1,k0-1} (if any) and {f10} of {k1-1,k0}. */

                /* Compute the domain coordinate {x11[0]} of upper right tile corner: */
                double u0 = ((double)k0)/((double)n0);
                x11[0] = (1 - u0)*LO(B[0]) + u0*HI(B[0]);

                /* Compute domain coordinate {xmm[0]} of tile center: */
                double v0 = ((double)k0 - 0.5)/((double)n0);
                xmm[0] = (1 - v0)*LO(B[0]) + v0*HI(B[0]);

                /* Compute the function values for the top right corner: */
                F(x11, ddim, f11, nf);

                /* Compute the function values for the center of the tile: */
                F(xmm, ddim, fmm, nf);

                if (k0 > 0)
                  { if(k1 > 0) 
                      { /* Paint the four triangles: */
                        pswr_plot_2D_tri_atom
                          ( ps, f00, f01, fmm, nf, st, fillP, drawP );

                        pswr_plot_2D_tri_atom                     
                          ( ps, f01, f11, fmm, nf, st, fillP, drawP );

                        pswr_plot_2D_tri_atom                     
                          ( ps, f11, f10, fmm, nf, st, fillP, drawP );

                        pswr_plot_2D_tri_atom                     
                          ( ps, f10, f00, fmm, nf, st, fillP, drawP );
                      }
                  }

                if (k0 > 0)
                  { /* Save top left corner as bottom left corner for next row: */
                    int j;
                    for (j = 0; j < nf; j++) { f00[j] = f01[j]; }
                  }
              }
            /* Save top right corner of last tile as last bottom right corner: */
            int j;
            for (j = 0; j < nf; j++) { f10[j] = f11[j]; }
          }
      }
  }
   
void pswr_plot_2D_tri
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    double xa[],
    double xb[],
    double xc[],
    int ns,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  )
  {
    /* The domain dimension must be 2: */
    int ddim = 2;
    
    /* We need at least 2 values for positional coordinates: */
    demand(nf >= 2, "too few function values to plot.");
    demand(F != NULL, "cannot plot a null function.");
    
    /* Number of tile corners along each scanline: */
    int nv = ns+1;   /* Number of corners. */
    double f[nv*nf]; /* Value vectors at corners of previous scanline. */  
    /* The values at tile corner {ib} are {f[ib*nf+j]}, for {j=0..nf-1}. */
    
    double fp[nf], fq[nf]; /* Work vectors for {F}-vals at corners of a 4-sided tile. */
    double x11[ddim];      /* Upper right corner of current chip pair. */
    
    int pass; /* Pass 0 is fill, pass 1 is draw: */
    for (pass = 0; pass < 2; pass++)
      { if ((pass == 0) && (! fill)) { continue; }
        if ((pass == 1) && (! draw)) { continue; }
        
        bool_t fillP = (pass == 0);
        bool_t drawP = (pass == 1);

        /* Scan the vertices (with overshoot of 1 on {kb}): */
        int ka, kb, kc; /* Indices of vertices (sum is {ns}, all non-negative when inside). */
        for (kc = 0; kc <= ns; kc++)
          { /* Now, if {kc > 0}, {f} contains the {F}-vals for the vertices in scanline {kc-1}. */

            /* Initialize the pointers to the {F}-vals on vertex row {kc} (will swap): */
            double *f01 = fp;
            double *f11 = fq;

            /* Initialize the pointers to the {F}-vals on row {kc-1}: */
            double *f00 = NULL;
            double *f10 = NULL;

            /* Scan the domain chips upper-bounded by scanline {kc}, in pairs: */
            for (kb = 0; kb <= ns-kc+1; kb++)
              { ka = ns - kb - kc;

                /* Now {f} contains the {F}-values for the vertices {0..kb-1} of scanline {kc}, */
                /* and vertices {kb..ns} of scanline {kc-1}. */
                /* Also, if {kb > 0}, */
                /*   {f10} points to a segment of {f}. */
                /*   {f01} and {f11} are scratch arrays. */
                /*   {f10} has the {F}-vals for row {kc-1}, column {kb-1}. */
                /*   {f11} has the {F}-vals of row {kc}, column {kb-1}. */

                /* We will plot the two chips between rows {kc-1} and {kc}, cols {kb-1} and {kb}. */
                /* On the last iteration only the first of these two chips is inside {T}. */
                /* If {kc == 0} the chips do not exist but we need to compute their top corners. */

                /* Swap the pointers to the top corners of the current chip pair: */
                { double *tmp = f01; f01 = f11; f11 = tmp; }
                /* Now {f01} is {F}-vals of vertex {kc,kb-1}, and {f11} is scratch area. */

                /* Update the pointers to the bottom corners of current tile: */
                if (kb > 0) { f00 = f10; } else { f00 = NULL; }
                if (kb <= ns) { f10 = &(f[kb*nf]); } else { f10 = NULL; }
                /* Now {f00} and {f10} point to segments of {f[..]}; */
                /* {f00} is the {F}-vals of {kc-1,kb-1} (if any) and {f10} of {kc-1,kb} (ditto). */

                if ((kb > 0) && (kc > 0)) 
                  { /* Paint the upward-pointing chip in chip pair: */
                    pswr_plot_2D_tri_atom
                      ( ps, f00, f10, f01, nf, st, fillP, drawP );
                  }

                if (ka >= 0)
                  { /* Compute domain coordinates {x11[0..1]} of vertex {kc,kb}: */
                    double ua = ((double)ka)/((double)ns);
                    double ub = ((double)kb)/((double)ns);
                    double uc = ((double)kc)/((double)ns);
                    x11[0] = ua*xa[0] + ub*xb[0] + uc*xc[0];
                    x11[1] = ua*xa[1] + ub*xb[1] + uc*xc[1];

                    /* Evaluate {F} at vertex {kc,kb}: */
                    F(x11, ddim, f11, nf);

                    if ((kb > 0) && (kc > 0))
                      { /* Paint the downward-pointing chip in chip pair: */
                        pswr_plot_2D_tri_atom
                          ( ps, f10, f11, f01, nf, st, fillP, drawP );
                      }
                  }

                if (kb > 0)
                  { /* Save top left corner as bottom left corner for next row: */
                    int j;
                    for (j = 0; j < nf; j++) { f00[j] = f01[j]; }
                  }
              }
          }
      }
  }
 
void pswr_plot_2D_quad_outline
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t B[],
    int n0,
    int n1
  )
  {
    /* The domain dimension must be 2: */
    int ddim = 2;
    
    double *x0; /* Either {xa} or {xb}. */
    double *x1; /* Either {xb} or {xa}. */
    
    auto void L(double z[], int nz, double g[], int ng);
      /* Maps a real {z[0]} from the interval {[0_1]} to 
        a point {x[0..1]} on the domain segment 
        from {x0[0..1]} to {x1[0..1]}, affinely; 
        then calls {F} on {xp} and returns the result
        in {g[0..ng-1]}. Requires {ng == nf}. */
        
    void L(double z[], int nz, double g[], int ng)
      { assert(nz == 1);
        assert(ng == nf);
        double x[ddim];
        x[0] = (1 - z[0])*x0[0] + z[0]*x1[0];
        x[1] = (1 - z[0])*x0[1] + z[0]*x1[1];
        F(x, ddim, g, ng);
      }
    
    double xa[ddim], xb[ddim]; /* Two adjacent corners of {B[0]×B[1]}. */
    interval_t U = (interval_t){{ 0.0, 1.0 }};  /* A unit interval for the parameter. */
    
    x0 = xa; x1 = xb;
    x0[0] = LO(B[0]); x0[1] = LO(B[1]); 
    x1[0] = HI(B[0]); x1[1] = LO(B[1]);
    pswr_plot_2D_line(ps, nf, L, &U, n0);
    
    { double *tmp = x0; x0 = x1; x1 = tmp; }
    x1[0] = HI(B[0]); x1[1] = HI(B[1]); 
    pswr_plot_2D_line(ps, nf, L, &U, n1);
    
    { double *tmp = x0; x0 = x1; x1 = tmp; }
    x1[0] = LO(B[0]); x1[1] = HI(B[1]); 
    pswr_plot_2D_line(ps, nf, L, &U, n0);
    
    { double *tmp = x0; x0 = x1; x1 = tmp; }
    x1[0] = LO(B[0]); x1[1] = LO(B[1]); 
    pswr_plot_2D_line(ps, nf, L, &U, n1);
  }

void pswr_plot_2D_tri_outline
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    double xa[],
    double xb[],
    double xc[],
    int ns
  )
  {
    /* The domain dimension must be 2: */
    int ddim = 2;
    
    double *x0; /* Either {xa}, {xb}, or {xc}. */
    double *x1; /* Either {xb}, {xc}, or {xa}. */
    
    auto void L(double z[], int nz, double g[], int ng);
      /* Maps a real {z[0]} from the interval {[0_1]} to 
        a point {x[0..1]} on the domain segment 
        from {x0[0..1]} to {x1[0..1]}, affinely; 
        then calls {F} on {xp} and returns the result
        in {g[0..ng-1]}. Requires {ng == nf}. */
        
    void L(double z[], int nz, double g[], int ng)
      { assert(nz == 1);
        assert(ng == nf);
        double x[ddim];
        x[0] = (1 - z[0])*x0[0] + z[0]*x1[0];
        x[1] = (1 - z[0])*x0[1] + z[0]*x1[1];
        F(x, ddim, g, ng);
      }
    
    interval_t U = (interval_t){{ 0.0, 1.0 }};  /* A unit interval for the parameter. */
    
    x0 = xa; x1 = xb;
    pswr_plot_2D_line(ps, nf, L, &U, ns);
    
    x0 = xb; x1 = xc;
    pswr_plot_2D_line(ps, nf, L, &U, ns);
    
    x0 = xc; x1 = xa;
    pswr_plot_2D_line(ps, nf, L, &U, ns);
  }

void pswr_plot_2D_line
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t *B,
    int ns
  )
  {
    /* The domain dimension must be 2: */
    int ddim = 1;
    
    /* We need at least 2 values for positional coordinates: */
    demand(nf >= 2, "too few function values to plot.");
    demand(F != NULL, "cannot plot a null function.");
    
    double fp[nf], fq[nf]; /* Function values at previous and next point. */

    /* Initialize the pointers to the {F}-vals at ends of the curr seg (will swap): */
    double *f0 = fp;
    double *f1 = fq;

    /* Scan the sample points: */
    double xp[ddim];  /* Endpoint of current segment. */
    int i;
    for (i = 0; i <= ns; i++)
      { /* Compute coordinate {xp[1]} of this point: */
        double u = ((double)i)/((double)ns);
        xp[0] = (1 - u)*LO(*B) + u*HI(*B);
        
        /* Swap the pointers to the {F}-vals at ends of the curr seg: */
        { double *tmp = f0; f0 = f1; f1 = tmp; }
        /* Now {f0} is the previous point, and {f1} is a scratch area. */
            
        /* Compute the function values for the current point: */
        F(xp, ddim, f1, nf);
            
        /* Plot segment, if we have both endpoints: */
        if (i > 0) { pswr_segment(ps, f0[0], f0[1], f1[0], f1[1]); }
      }
  }

void pswr_plot_2D_tri_atom
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  )
  { 
    if (fill)
      {
        if ((st->nc > 0) && (st->ic >= 0) && (st->ic + st->nc <= nf))
          { pswr_plot_2D_tri_atom_shade
              ( ps, fa, fb, fc, nf, st->ic, st->nc );
          }
        else if ((st->ib >= 0) && (st->ib < nf))
          { pswr_plot_2D_tri_atom_bands
              ( ps, fa, fb, fc, nf, st->ib, 
                st->vStart,st->vStep, st->kMin,st->kMax, 
                st->Rtb,st->Gtb,st->Btb
              );
          }
        else
          { pswr_plot_2D_tri_atom_solid
              ( ps, fa, fb, fc, nf );
          }
      }
    
    if (draw)
      { if ((st->iv >= 0) && (st->iv < nf))
          { pswr_plot_2D_tri_atom_isolines
              ( ps, fa, fb, fc, nf, st->iv, 
                st->vStart,st->vStep, st->kMin,st->kMax
              );
          }
      }
  }


void pswr_plot_2D_tri_atom_solid
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf 
  )  
  { demand(nf >= 2, "can't plot without coordinates"); 
    pswr_triangle
      ( ps,
        fa[0], fa[1],
        fb[0], fb[1],
        fc[0], fc[1],
        TRUE, FALSE
      );
  }

void pswr_plot_2D_tri_atom_shade
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int ic,
    int nc
  )  
  { demand(nf >= 2, "can't plot without coordinates"); 
    demand(nc > 0, "can't smooth-shade without function values"); 
    demand((ic >= 0) && (ic + nc <= nf), "invalid function value indices"); 
    /* Map corner function values to corner colors: */
    double *f[3] = { fa, fb, fc };
    double R[3], G[3], B[3]; /* Vertex colors. */
    int k;
    for (k = 0; k < 3; k++)
      { 
        /* Get function values for vertex {k} of part. */
        double *fk = f[k];
        /* Choose the color for vertex {k}. */
        double Rs = 1.000, Gs = 0.342, Bs = 0.000;
        double Y0 = 0.95;
        double Y1 = 0.30;
        double Rc, Gc, Bc;
        if (nc == 1)
          { pswr_color_scale_1(fk[ic], Rs, Gs, Bs, Y0, &(Rc), &(Gc), &(Bc)); }
        else if (nc == 2)
          { pswr_color_scale_2(fk[ic], fk[ic+1], Rs, Gs, Bs, Y0, &(Rc), &(Gc), &(Bc)); }
        else if (nc >= 3)
          { pswr_color_scale_3(fk[ic], fk[ic+1], fk[ic+2], Y0, Y1, &(Rc), &(Gc), &(Bc)); }
        frgb_t RGB = (frgb_t){{ (float)Rc, (float)Gc, (float)Bc }};
        frgb_clip_rgb(&RGB);
        R[k] = RGB.c[0];
        G[k] = RGB.c[1];
        B[k] = RGB.c[2];
      }
    /* Paint quadrilateral with bilinear shading: */
    pswr_shade_triangle
      ( ps,
        fa[0], fa[1], R[0], G[0], B[0],
        fb[0], fb[1], R[1], G[1], B[1],
        fc[0], fc[1], R[2], G[2], B[2],
        -1
      );
  }
  
void pswr_plot_2D_tri_atom_bands
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int ib,
    double vStart,
    double vStep, 
    int kMin,
    int kMax,
    double *Rtb, 
    double *Gtb, 
    double *Btb 
  )  
  { demand(nf >= 2, "can't plot without coordinates"); 
    demand((ib >= 0) && (ib < nf), "invalid function value index"); 
    if ((Rtb != NULL) && (Gtb != NULL) && (Btb != NULL))
      { pswr_bands_in_triangle
          ( ps,
            fa[0], fa[1], fa[ib],
            fb[0], fb[1], fb[ib],
            fc[0], fc[1], fc[ib],
            vStart, vStep, kMin, kMax,
            Rtb, Gtb, Btb
          );
      }
  }
  
void pswr_plot_2D_tri_atom_isolines
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int iv,
    double vStart,
    double vStep, 
    int kMin,
    int kMax 
  )  
  { demand(nf >= 2, "can't plot without coordinates"); 
    demand((iv >= 0) && (iv < nf), "invalid function value index"); 
    pswr_isolines_in_triangle
      ( ps,
        fa[0], fa[1], fa[iv],
        fb[0], fb[1], fb[iv],
        fc[0], fc[1], fc[iv],
        vStart, vStep, kMin, kMax
      );
  }

