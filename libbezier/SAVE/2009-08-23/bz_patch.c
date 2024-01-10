/* See bezier.h */
/* Last edited on 2009-08-23 13:24:59 by stolfi */

#define _GNU_SOURCE
#include <bz_basic.h>
#include <bz_patch.h>

#include <affirm.h>
#include <interval.h>
#include <box.h>
#include <bool.h>
#include <jsmath.h>

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

void bz_patch_debug(char *pref, bz_patch_t *b, char *suff); 
void bz_patch_point_debug(char *pref, double *x, bz_patch_ddim_t n, char *suff);

double bz_patch_test_func(int j, double p[]);
  /* Component {j} of the dependent function of the domain point {p}. */

bz_patch_t bz_patch_new(bz_patch_ddim_t m, bz_patch_rdim_t n, const bz_degree_t g[])
  { bz_patch_t b;
    int tsz, i;
    /* Compute total size {tsz} of control array: */
    /* fprintf(stderr, "+ bz_patch_new\n");  */
    tsz = n; for (i = 0; i < m; i++) { tsz *= (g[i]+1); }
    b.m = m;
    b.n = n;
    for (i = 0; i < bz_patch_MAX_DDIM; i++) { b.g[i] = ( i < m ? g[i] : 0); }
    b.c = (double *)notnull(malloc(tsz*sizeof(double)), "out of mem");
    /* fprintf(stderr, "- bz_patch_new %d[%d]\n", (int)(b.c), tsz);  */ 
    return b;
  }

bz_patch_t bz_patch_uniform_new(bz_patch_ddim_t m, bz_patch_rdim_t n, bz_degree_t g)
  { int i;
    bz_degree_t ug[bz_patch_MAX_DDIM];
    for (i = 0; i < m; i++) { ug[i] = g; }
    return bz_patch_new(m, n, ug);
  }
  
bz_patch_t bz_patch_facet_new
  ( bz_patch_ddim_t m, 
    bz_patch_rdim_t n, 
    const bz_degree_t g[], 
    box_axis_index_t j
  )
  { bz_degree_t fg[bz_patch_MAX_DDIM];
    int bi, fi;
    for (bi = 0, fi = 0; bi < m; bi++)
      { if (bi != j) { fg[fi] = g[bi]; fi++; } }
    return bz_patch_new(m-1, n, fg);
  }

void bz_patch_free(bz_patch_t b)
  {
    free(b.c);
  }

double *bz_patch_control_point(bz_patch_t *b, bz_patch_cpindex_t e[])
  { int ix, i;
    /* Compute index {ix} of point in point array, by the formula
      {ix = SUM{ PROD{(g[j]+1) : j=i+1..m-1}*e[i] :  i = 0..m-1}} */
    ix = 0; for (i = 0; i < b->m; i++) { ix = ix*(b->g[i]+1) + e[i]; }
    /* Return address of first coordinate; */
    return &(b->c[ix * b->n]);
  }

void bz_patch_eval(bz_patch_t *b, double x[], double sx[])
  { bz_patch_eval_m(b->m, b->n, b->g, b->c, x, sx); }
  
void bz_patch_eval_box(bz_patch_t *b, interval_t box[], bz_patch_t *f)
  { affirm(FALSE, "not implemented"); }
  
void bz_patch_compute_steps(
    bz_patch_t *b,  /* A Bézier patch. */
    int step[],   /* (OUT) control point index increment for each axis. */
    int *sz       /* (OUT) total number of control points. */
  )
  {
    int tsz = 1, i;
    bz_degree_t *g = &(b->g[0]);
    for (i = b->m-1; i >= 0; i--) { step[i] = tsz; tsz *= (g[i] + 1); }
    (*sz) = tsz;
  }

void bz_patch_eval_m
  ( bz_patch_ddim_t m,           /* Dimension of parameter space. */
    bz_patch_rdim_t n,           /* Dimension of image space. */
    const bz_degree_t g[], /* Degrees on each coordinate. */
    double *c,           /* Bezier control points. */
    double x[],          /* Argument vector (size {m}). */
    double fx[]          /* OUT: Image vector (size {n}). */
  )
  { /* Recursion on the argument dimension {m}: */
    /* fprintf(stderr, "+ bz_patch_eval m = %d", m); */
    /* bz_point_debug(" x = (", x, m, ")\n"); */
    if (m == 0)
      { /* A Bezier patch from {R^0} to {R^n} is just a point of {R}: */
        int i; for (i = 0; i < n; i++) { fx[i] = c[i]; }
      }
    else
      { double rc[(g[0]+1)*n];
        int tsz1, i; 
        /* Compute size of control array for an {(m-1)}-dimensional patch: */
        tsz1 = n; for (i = 1; i < m; i++) { tsz1 *= (g[i]+1); }
        /* Recursively evaluate {g+1} patches of dimension {m-1}: */
        for (i = 0; i <= g[0]; i++)
          { bz_patch_eval_m(m-1, n, &(g[1]), &(c[tsz1*i]), &(x[1]), &(rc[n*i])); }
        /* Now do a 1-dimensional DeCasteljau interpolation based on {x[0]}: */
        bz_patch_eval_1(n, g[0], rc, x[0], fx);
      }
    /* fprintf(stderr, "- bz_patch_eval_m  m = %d", m); */
    /* bz_point_debug(" fx = (", fx, n, ")\n"); */
  }
  
void bz_patch_eval_1
  ( bz_patch_rdim_t n,     /* Dimension of image space. */
    bz_degree_t g,   /* Degree of curve. */
    double *c,    /* Bezier control points. */
    double x,      /* Argument value. */
    double fx[]    /* OUT: Image vector (size {n}). */
  )
  { /* The DeCasteljau algorithm: */
    int i, s, r;
    for (i = 0; i < n; i++)
      { double y = 1.0 - x;
        for (r = g; r > 0; r--)
          { for (s = 0; s < r; s++)
              { c[n*s + i] = (y*c[n*s + i] + x*c[n*(s+1) + i]); }
          }
        fx[i] = c[i];
      }
  }

void bz_patch_get_facet(bz_patch_t *b, box_axis_t ax, interval_side_t dir, bz_patch_t *t)
  { int tm = t->m;
    int bm = b->m;
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    affirm(t->n == b->n, "range dimension mismatch");
    affirm(tm == bm-1, "domain dimension mismatch");
    { int it, ib;
      for (it = 0, ib = 0; it < tm; it++,ib++)
        { if (ib == ax) { ib++; }
          affirm(tg[it] == bg[ib], "degree mismatch");
        }
    }
    { int n = t->n;
      static int bstep[bz_patch_MAX_DDIM];
      static int tstep[bz_patch_MAX_DDIM];
      int bsz, tsz, tix;
      /* Computes the index increments and total sizes for the control point arrays: */ 
      bz_patch_compute_steps(b, bstep, &bsz);
      bz_patch_compute_steps(t, tstep, &tsz);
      /* Copy control points: */
      for (tix = 0; tix < tsz; tix++)
        { /* Compute index {bix} such that {b->c[bix]} corresponds to {t->c[tix]}: */
          int bix = 0, tt = tix, it, ib;
          for (ib = bm-1, it = tm-1; ib >= 0; ib--) 
            { if (ib == ax) 
                { bix += (dir == 1 ? bg[ax]*bstep[ax] : 0); }
              else
                { int tnumi = tg[it]+1, tei = tt % tnumi;
                  bix += tei*bstep[ib];
                  tt /= tnumi; it--;
                }
            }
          /* Copy the coefficient vector: */
          { double *tc = &(t->c[n*tix]), *bc = &(b->c[n*bix]); 
            int j;
            for (j = 0; j < n; j++,tc++,bc++) { (*tc) = (*bc); }
          }
        }
    }
  }

void bz_patch_get_face
  ( bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  )
  { int tm = t->m;
    int bm = b->m;
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    affirm(t->n == b->n, "range dimension mismatch");
    { int it, ib;
      for (ib = 0,it = 0; ib < b->m; ib++)
        { if (dir[ib] == 00) 
            { affirm(it < tm, "dimension too low for location");
              affirm(tg[it] == bg[ib], "degree mismatch");
              it++; 
            }
        }
      affirm(tm == it, "dimension inconsistent with location");
    }
    { int n = b->n;
      static int bstep[bz_patch_MAX_DDIM];
      static int tstep[bz_patch_MAX_DDIM];
      int bsz, tsz, tix;
      
      bz_patch_compute_steps(b, bstep, &bsz);
      bz_patch_compute_steps(t, tstep, &tsz);
      /* Copy control points: */
      for (tix = 0; tix < tsz; tix++)
        { /* Compute index {bix} such that {b->c[bix]} corresponds to {t->c[tix]}: */
          int bix = 0, tt = tix, it, ib;
          for (ib = bm-1, it = tm-1; ib >= 0; ib--) 
            { if (dir[ib] == -1) 
                { /* bix += 0*bstep[ib]; */ }
              else if (dir[ib] == +1) 
                { bix += bg[ib]*bstep[ib]; }
              else
                { int tnumi = tg[it]+1, tei = tt % tnumi;
                  bix += tei*bstep[ib];
                  tt /= tnumi; it--;
                }
            }
          /* Copy the coefficient vector: */
          { double *tc = &(t->c[n*tix]), *bc = &(b->c[n*bix]); 
            int j;
            for (j = 0; j < n; j++,tc++,bc++) { (*tc) = (*bc); }
          }
        }
    }
  }

void bz_patch_set_face
  ( bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  )
  { int tm = t->m;
    int bm = b->m;
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    affirm(t->n == b->n, "range dimension mismatch");
    { int it, ib;
      for (ib = 0,it = 0; ib < b->m; ib++)
        { if (dir[ib] == 00) 
            { affirm(tg[it] == bg[ib], "degree mismatch");
              it++; 
            }
        }
      affirm(tm == it, "domain dimension mismatch");
    }
    { int n = b->n;
      static int bstep[bz_patch_MAX_DDIM];
      static int tstep[bz_patch_MAX_DDIM];
      int bsz, tsz, bix;
      
      bz_patch_compute_steps(b, bstep, &bsz);
      bz_patch_compute_steps(t, tstep, &tsz);
      /* Modify control points: */
      for (bix = 0; bix < bsz; bix++)
        { /* Compute {tix} such that {t->c[tix]} is  the control point of
            {t} closest to {b->c[ix]}; compute the index {btix} such that
            {b->c[btix]} is the corresponding control point in {b};
            and compute the influence factor {dfac} of that control 
            point on {b->c[bix]}: */
          int bb = bix, tix = 0, btix = 0, it, ib;
          double dfac = 1.0;
          for (ib = bm-1, it = tm-1; ib >= 0; ib--) 
            { int bnumi = bg[ib] + 1, bei = bb % bnumi;
              double xi = ((double)bei)/((double)bg[ib]);
              bb /= bnumi;
              if (dir[ib] == -1)
                { dfac *= (1.0 - xi); /* btix += 0*bstep[ib]; */ }
              else if (dir[ib] == +1)
                { dfac *= xi; btix += bg[ib]*bstep[ib]; }
              else
                { btix += bei*bstep[ib]; tix += bei*tstep[it]; it--; }
            }
          if (dfac > 0.0)
            { /* Modify {b->c[bix]} by {t->c[tix]-b->c[btix]} times {dfac}: */
              double *bc = &(b->c[n*bix]);
              double *tc = &(t->c[n*tix]);
              int j;
              if (dfac == 1.0)
                { affirm(bix == btix, "indexing bug");
                  for (j = 0; j < n; j++,tc++,bc++) { (*bc) = (*tc); }
                }
              else
                { double *btc = &(b->c[n*btix]);
                  affirm(bix != btix, "indexing bug");
                  for (j = 0; j < n; j++,tc++,bc++,btc++)
                    { double delta = dfac*((*tc) - (*btc));
                      (*bc) += delta;
                    } 
                }
              
            }
        }
    }
  }
  
bz_degree_t bz_patch_max_degree(bz_patch_t *b)
  { bz_degree_t gmax = -1;
    bz_degree_t *bg = &(b->g[0]);
    bz_patch_ddim_t bm = b->m;
    while(bm > 0) { bm--; if (bg[bm] > gmax) { gmax = bg[bm]; } }
    return gmax;
  }

void bz_patch_raise_degree(bz_patch_t *b, bz_patch_t *t)
  { /* Assumption: Linear interpolation of the indices? */
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    affirm(t->m == b->m, "domain dimension mismatch");
    affirm(t->n == b->n, "range dimension mismatch");
    { int i;
      for (i = 0; i < b->m; i++)
        { affirm(tg[i] >= bg[i], "degree cannot be lowered");
          affirm(bg[i] == 1, "raising implemented only from degree 1");
        }
    }
    { int m = b->m;
      int n = b->n;
      static double x[bz_patch_MAX_DDIM];
      static int bstep[bz_patch_MAX_DDIM];
      static int tstep[bz_patch_MAX_DDIM];
      int bsz, tsz, tix;
      bz_patch_compute_steps(b, bstep, &bsz);
      bz_patch_compute_steps(t, tstep, &tsz);
      /* Interpolate control points of {t} from those of {b}: */
      for (tix = 0; tix < tsz; tix++)
        { /* Compute nominal position {x[0..m-1]} of control point {t->c[tix]}: */
          int tt = tix, i;
          for (i = m - 1; i >= 0; i--) 
            { int tnumi = tg[i] + 1, tei = tt % tnumi;
              x[i] = ((double)tei)/((double)tg[i]);
              tt /= tnumi;
            }
          /* Evaluate {b} at {x} to obtain the control point: */
          bz_patch_eval(b, x, &(t->c[n*tix]));
        }
    }
  }

void bz_patch_invert
  ( double p[],
    bz_patch_t *b,
    double tol,
    double x[]
  )
  { /* Subdivide the shape by DeCasteljau (possibly unequal), until */
    /* it can be approximated by an affine function, then invert that. */
    affirm(FALSE, "not implemented");
  }

void bz_patch_compute_bbox(bz_patch_t *b, interval_t box[])
  { /* Uses the fact that the cell is contained in the hull of the control pts. */
    int tsz, i, j, t;
    int m = b->m, n = b->n;
    double *c = b->c;
    /* Compute total size {tsz} of the control point array: */
    tsz = n; for (i = 0; i < m; i++) { tsz *= (b->g[i]+1); }
    /* Now find the control point range in each range coordinate {j}: */
    for (j = 0; j < n; j++)
      { double lo = INFINITY;
        double hi = -INFINITY;
        for (t = j; t < tsz; t += n)
          { double ci = c[t];
            if (ci < lo) { lo = ci; }
            if (ci > hi) { hi = ci; }
          }
        box[j] = (interval_t){{lo, hi}};
      }
  }
  
void bz_patch_multiaffine_approx(bz_patch_t *b, bz_patch_t *t)
  { /* Just get the corner control points. */
    /* Then subtract from original to bound the error. */
    int m = b->m;
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    bool_t trivial = TRUE;
    affirm(t->m == b->m, "domain dimension mismatch");
    affirm(t->n == b->n, "range dimension mismatch");
    { int i;
      for (i = 0; i < m; i++)
        { affirm(tg[i] == 1, "approximant is not multiaffine"); 
          if (bg[i] != 1) { trivial = FALSE; }
        }
    }
    { int n = b->n;
      int tsz = ipow(2, m);
      if (trivial)
        { int ij;
          double *tc = &(t->c[0]), *bc = &(b->c[0]); 
          /* Copy the control points literally: */
          for (ij = 0; ij < n*tsz; ij++,tc++,bc++) { (*tc) = (*bc); }
        }
      else
        { static int bstep[bz_patch_MAX_DDIM];
          int bsz, tix;
          bz_patch_compute_steps(b, bstep, &bsz);
          /* Copy the corner control points: */
          /* The indices {tix} and {bix} do not include the factor {n}. */
          for (tix = 0; tix < tsz; tix++)
            { /* Compute index {bix} s.t. {b->c[bix]} corresponds to {t->c[tix]}: */
              int bix = 0, tt = tix, i;
              for (i = m-1; i >= 0; i--) 
                { int tei = tt & 1;
                  bix += tei*bg[i]*bstep[i];
                  tt >>= 1;
                }
              { double *tc = &(t->c[n*tix]), *bc = &(b->c[n*bix]); 
                int j;
                for (j = 0; j < n; j++,tc++,bc++) { (*tc) = (*bc); }
              }
            }
        }
    }
  }
  
double bz_patch_multiaffine_error(bz_patch_t *b, bz_patch_t *t)
  { /* Just get the corner control points. */
    /* Then subtract from original to bound the error. */
    const bz_degree_t *bg = &(b->g[0]);
    const bz_degree_t *tg = &(t->g[0]);
    bool_t trivial = TRUE;
    affirm(t->m == b->m, "domain dimension mismatch");
    affirm(t->n == b->n, "range dimension mismatch");
    { int i;
      for (i = 0; i < b->m; i++)
        { affirm(tg[i] == 1, "approximant is not multiaffine"); 
          if (bg[i] != 1) { trivial = FALSE; }
        }
    }
    if (trivial)
      { return 0.0; }
    else
      { int m = b->m;
        int n = b->n;
        static double x[bz_patch_MAX_DDIM];
        static double tx[bz_patch_MAX_DDIM];
        static int bstep[bz_patch_MAX_DDIM];
        int bsz, bix, i;
        double dmax = 0.0;
        bz_patch_compute_steps(b, bstep, &bsz);
        for (bix = 0; bix < bsz; bix++)
          { /* Compute nominal coordinates {x} of control point {b->c[bix]}: */
            int bb = bix;
            for (i = m-1; i >= 0; i--) 
              { int bnumi = bg[i] + 1, bei = bb % bnumi;
                x[i] = ((double)bei)/((double)bg[i]);
                bb /= bnumi;
              }
            /* Evaluate the nultiaffine spline {t} at {x}: */
            bz_patch_eval(t, x, tx);
            /* Compare control point {b->c[bix]} with interpolated point {tx}: */
            { double *tc = tx, *bc = &(b->c[n*bix]); 
              int j;
              for (j = 0; j < n; j++,tc++,bc++) 
                { double dj = fabs((*tc) - (*bc));
                  if (dj > dmax) { dmax = dj; }
                }
            }
          }
        return dmax;
      }
  }

void bz_patch_affine_approx(bz_patch_t *b, double c[], double M[], double *err)
  { /* Compute approximate center, mean derivatives. */
    /* Then re-convert to Bézier and subtract from original to bound the error. */
    affirm(FALSE, "not implemented");
  }
  
void bz_patch_grad(bz_patch_t *b, box_axis_t ax, bz_patch_t *t)
  { affirm(FALSE, "not implemented"); }

void bz_patch_split
  ( bz_patch_t *b, 
    box_axis_t a, 
    double ratio,
    bz_patch_t *bLO,
    bz_patch_t *bHI
  )
  { int m = b->m, n = b->n;
    const bz_degree_t *g = &(b->g[0]);
    double *c = b->c, *cLO = bLO->c, *cHI = bHI->c;
    
    /* The DeCasteljau interpolation must act on the sequence of {g[a]+1}
      points {c'[t] = c[e[0],..e[m-1]]} where all {e[i]} stay fixed
      except {t = e[a]} which varies from {0} to {g[a]}. 
      
      Recall that the index {r} of control point {c[e[0],..e[m-1]]} 
      in the linearized control point array is
      
        {r = SUM{ G[k]*e[k] :  k = 0..m-1}}
        
      where {G[k] = PROD{g[j]+1 : j = k+1..m-1}}.
      
      Breaking this sum at the term {k = a} we get {r = u + v + w}, where
      
        {u = SUM{ G[k]*e[k] :  k = 0..a-1}}
        {v = G[a]*e[a]}
        {w = SUM{ G[k]*e[k] :  k = a+1..m-1}}
        
      Extracting out the common factors, we get
      
        {u = G[a-1] * SUM{ PROD{g[j]+1 : j = k+1..a-1}*e[k] :  k = 0..a-1}}
        
      As we vary the {e[i]} from {0} to {g[i]}, the terms {u,v,w} vary 
      from 0 (inclusive) to limits {limu,limv,limw} (exclusive)
      with steps {du,dv,dw}, where 
      
        {limw = G[a]}           {dw = 1}
        {limv = G[a-1]}         {dv = G[a] = limw}
        {limu = G[m-1]}         {du = G[a-1] = limv}
        
      Keep in mind that {r = u + v + w} is still a *point* index, and
      each point has {n} coordinates. So {r} must be multiplied by
      {n}, and added to the component index {i} which ranges
      over {0..n-1}. */
    
    int u, v, v1, w, i, k, limu, limv, limw, du, dv, dw;
    double s = 1.0 - ratio, t = ratio;
    
    /* fprintf(stderr, "+ pz_patch_split\n"); */ 
    /* Compute {limu,limv,limw} (already times {n}): */
    limw = n; for(k = a+1; k < m; k++) { limw*= (g[k]+1); }
    limv = limw * (g[a]+1);
    limu = limv; for(k = 0; k < a; k++) { limu*= (g[k]+1); }
    /* Compute the increments {du,dv,dw} (already times {n}): */
    du = limv; dv = limw; dw = n;
    
    /* Now do it: */
    /* fprintf(stderr, "= pz_patch_split cLO = %d cHI = %d\n", (int)cLO, (int)cHI); */ 
    /* fprintf(stderr, "limu = %6d  limv = %6d  limw = %6d\n", limu, limv, limw); */
    /* fprintf(stderr, "du =   %6d  dv =   %6d  dw =   %6d\n", du, dv, dw);  */
    for (u = 0; u < limu; u += du)
      { for (w = 0; w < limw; w += dw)
          { for (i = 0; i < n; i++)
              { int uwi = u + w + i;
                /* DeCasteljau algorithm along {v} index: */
                for (v = 0; v < limv; v += dv) 
                  { cHI[uwi + v] = c[uwi + v]; }
                for (v = 0; v < limv; v += dv)
                  { cLO[uwi + v] = cHI[uwi];
                    for (v1 = 0; v1 < limv - v - dv; v1 += dv)
                      { cHI[uwi + v1] = s*cHI[uwi + v1] + t*cHI[uwi + v1 + dv]; }
                  }
              }
          }
      }
    /* fprintf(stderr, "- pz_patch_split\n"); */ 
  }

double bz_patch_try_flatten_face(
    bz_patch_t *b, 
    box_signed_dir_t dir[], 
    double tol
  )
  { /* Data for Bézier patch {b}: */
    int bm = b->m;
    bz_degree_t *bg = &(b->g[0]);
    /* Data for its face {F}: */
    int tm;
    bz_degree_t tg[bz_patch_MAX_DDIM];
    double err;
    int trivial = TRUE;

    /* Get signature, dimension, and degree vector of face: */
    { int ib;
      for (ib = 0, tm = 0; ib < bm; ib++)
        { if (dir[ib] == 00) 
            { tg[tm] = bg[ib]; tm++;
              if (bg[ib] != 1) { trivial = FALSE; } 
            } 
        }
    }
    if (trivial)
      { /* Face is multi-affine, hence already flat: */
        return 0.0;
      } 
    else
      { /* Could be done more efficiently... */
        bz_patch_t t = bz_patch_new(tm, b->n, tg);
        bz_patch_t a = bz_patch_uniform_new(tm, b->n, 1);
        fprintf(stderr, "+ bz_try_flatten_face\n");
        bz_patch_debug("b = \n", b, "\n");
        bz_patch_get_face(b, dir, &t);
        bz_patch_debug("orig t = \n", &t, "\n");
        bz_patch_multiaffine_approx(&t, &a);
        bz_patch_debug("flat a = \n", &a, "\n");
        err = bz_patch_multiaffine_error(&t, &a);
        if (err <= tol)
          { /* Flatten face {f} of {b}, propagate correction to superfaces: */
            fprintf(stderr, "flattening the face:\n");
            bz_patch_raise_degree(&a, &t);
            bz_patch_debug("flat t = \n", &t, "\n");
            bz_patch_set_face(b, dir, &t);
            bz_patch_debug("modf b = \n", b, "\n");
            err = 0.0; 
          }
        free(t.c);
        free(a.c);
        fprintf(stderr, "- bz_patch_try_flatten_face\n");
        return err;
      }
  }

void bz_patch_print(FILE *wr, bz_patch_t *b, char *fmt)
  {
    int ix, j;
    int rowsz, planesz, bsz;
    double *bc = &(b->c[0]);
    int n = b->n;
    { int i; bsz = 1; for (i = 0; i < b->m; i++) { bsz *= b->g[i]+1; } }
    rowsz = ( b->m >= 1 ? b->g[b->m-1]+1 : bsz + 1 );
    planesz = ( b->m >= 2 ? b->g[b->m-2]+1 : bsz + 1 );
    for (ix = 0; ix < bsz; ix++)
      { for (j = 0; j < n; j++,bc++)
          { if (j != 0) { fprintf(wr, " "); } fprintf(wr, fmt, (*bc)); }
        fprintf(wr, ((ix+1) % rowsz == 0 ? "\n" : "  "));
        if ((ix+1) % planesz == 0) { fprintf(wr, "\n");  }
      }
  }

bz_patch_t bz_patch_from_box(bz_patch_ddim_t d, interval_t box[])
  { 
    /* fprintf(stderr, "+ dg_patch_from_box\n"); */
    bz_patch_t b = bz_patch_uniform_new(d, d, 1);
    int ix; 
    int bsz = ipow(2, d);
    for (ix = 0; ix < bsz; ix++)
      { /* Get labels {e[0..d-1]} of Bezier coeff {b.c[ix]} */
        int t = ix, i;
        for (i = 0; i < d; i++) 
          { int ei = (t % 2); t /= 2;
            b.c[d*ix + i] = box[i].end[ei];
          }
      }
    /* fprintf(stderr, "- dg_patch_from_box\n"); */
    return b;
  }

void bz_patch_debug(char *pref, bz_patch_t *b, char *suff)
  {
    fprintf(stderr, "%s", pref);
    bz_patch_print(stderr, b, "%6.2f");
    fprintf(stderr, "%s", suff);
  }

void bz_patch_point_debug(char *pref, double x[], bz_patch_ddim_t n, char *suff)
  { int j;
    fprintf(stderr, "%s", pref);
    for (j = 0; j < n; j++,x++)
      { if (j != 0) { fprintf(stderr, " "); } fprintf(stderr, "%6.2f", (*x)); }
    fprintf(stderr, "%s", suff);
  }

bz_patch_t bz_patch_make_test(bz_patch_ddim_t d, bz_patch_rdim_t n, bz_degree_t g)
  { 
    bz_patch_t b = bz_patch_uniform_new(d, n, g);
    int ix; 
    int g1 = g+1;
    int bsz = ipow(g1, d);
    
    fprintf(stderr, "+ bz_patch_make_test\n");
    if (g == 1)
      { 
        /* An irregular quadrilateral: */
        affirm(d == 2, "implemented only for 2D domains");
        for (ix = 0; ix < bsz; ix++)
          { int e[d];
            double p[d], w;
            /* Get labels {e[0..d-1]} of Bezier coeff {b.c[ix]} */
            int t = ix, i, j;
            fprintf(stderr, "e = (");
            for (i = 0; i < d; i++) 
              { e[i] = (t % g1); t /= g1; fprintf(stderr, "%d", e[i]); }
            fprintf(stderr, ")\n");
            w = 1.0 + 0.3*e[0];
            p[0] = w*(2.0*e[0] - 3.0*e[1] + 3.0);
            p[1] = w*(3.0*e[0] + 2.0*e[1]);
            /* Compute function values {p[d..n-1] at {p[0..d-1]}: */
            for (j = d; j < n; j++) { p[j] = bz_patch_test_func(j-d, p); }
            fprintf(stderr, "p = (");
            for (j = 0; j < n; j++) { fprintf(stderr, " %.6f", p[j]); }
            fprintf(stderr, " )\n");
            /* Store in {b.c}: */
            for (j = 0; j < n; j++) { b.c[n*ix + j] = p[j]; }
          }
      }
    else if (g == 2)
      { 
        affirm(d >= 2, "domain must be at least 2D");
        /* Build Bézier patch as rectangle with parabolic bend: */
        for (ix = 0; ix < bsz; ix++)
          { int e[d];
            double p[n];
            /* Get labels {e[0..d-1]} of Bezier coeff {b.c[ix]} */
            int t = ix, i, j;
            for (i = 0; i < d; i++) { e[i] = (t % g1); t /= g1; }
            if ((e[0] == 1) && (e[1] == 2))
              { p[0] = 3.0; p[1] = 3.0; }
            else
              { p[0] = 2.0*e[0]; p[1] = e[1]; }
            /* Compute function values {p[d..n-1] at {p[0..d-1]}: */
            for (j = d; j < n; j++) { p[j] = bz_patch_test_func(j-d, p); }
            /* Store in {b.c}: */
            for (j = 0; j < n; j++) { b.c[n*ix + j] = p[j]; }
          }
      }
    else if (g == 3)
      {
        double fd = (double)d;
        double ALPHA = 0.75*(sqrt(2)-1); /* Bézier constant for circle appr. */
        affirm(d >= 2, "domain must be at least 2D");
        /* Build Bézier patch approximating one quadrant of a thick sphere: */
        for (ix = 0; ix < bsz; ix++)
          { int e[d];
            double p[n];
            /* Get labels {e[0..d-1]} of Bezier coeff {b.c[ix]} */
            int t = ix, i, j;
            for (i = 0; i < d; i++) { e[i] = (t % g1); t /= g1; }
            { /* Bézier coeff for {y = x^d}, for {x} in {[(1/2)^{1/d} _ 1]}: */
              double rmin = 0.5, smin = pow(rmin,1.0/fd), vmin = fd*pow(smin,fd-1);
              double rmax = 1.0, smax = pow(rmax,1.0/fd), vmax = fd*pow(smax,fd-1);
              switch(e[0])
                { case 0: p[0] = rmin; break;
                  case 1: p[0] = rmin + vmin*(smax - smin)/3.0; break;
                  case 2: p[0] = rmax - vmax*(smax - smin)/3.0; break;
                  case 3: p[0] = rmax; break;
                  default: affirm(FALSE, "bad Bezier label");
                }
              fprintf(stderr, "  layer = %d radius = %6.4f\n", e[0], p[0]);
              /* Bézier coeffs for circle */
              for (i = 1; i < d; i++) 
                { double x = p[i-1], nx, ny;
                  switch (e[i])
                    { case 0: nx = x; ny = 0.0; break;
                      case 1: nx = x; ny = ALPHA*x; break;
                      case 2: nx = ALPHA*x; ny = x; break;
                      case 3: nx = 0.0; ny = x; break;
                      default: /* GCC pacifier: */ nx = ny = 0;
                    }
                  p[i-1] = nx; p[i] = ny;
                }
            }
            /* Compute function values {p[d..n-1] at {p[0..d-1]}: */
            for (j = d; j < n; j++) { p[j] = bz_patch_test_func(j-d, p); }
            /* Store in {b.c}: */
            for (j = 0; j < n; j++) { b.c[n*ix + j] = p[j]; }
          }
      }
    else
      { affirm(FALSE, "not implemented for this degree"); }
    /* Phew! */
    fprintf(stderr, "- bz_patch_make_test\n");
    return b;
  }


/* The Archimedean perimeter parameter: */
#define PI (M_PI)

double bz_patch_test_func(int j, double p[])
  { switch(j % 3)
    { case 0:  return sin(3.0*M_PI*(p[0]+p[1]+0.25));
      case 1:  return cos(3.0*M_PI*(p[0]-p[1]+0.25));
      case 2:  return cos(3.0*M_PI*(p[0]-p[1]+0.25));
      default: affirm(FALSE, "duh??");
        return 0.0;
    }
  }
