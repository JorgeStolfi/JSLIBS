/* See cpk_build.h */
/* Last edited on 2025-01-01 03:14:29 by stolfi */ 

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <r2x2.h>
#include <r2.h>
#include <affirm.h>
#include <interval.h>
#include <bool.h>

#include <hxg_paint.h>
#include <hxg_canvas.h>
#include <cpk_debug.h>
#include <cpk_io.h>
#include <cpk_graph.h>
#include <cpk_weight.h>
#include <cpk_build.h>


/* INTERNAL PROTOTYPES */

r2x2_t cpk_get_rotation(r2_t *dir);
  /* Returns the 2x2 rotation matrix {M} that takes the 
     vector {dir} to the X-axis. */

cpk_domain_t cpk_map_domain(cpk_domain_t *C, r2_t *o, r2x2_t *M);
  /* Modifies the domain specs {C} so that the set of {C}-valid points
    is translated by {-o}, and rotated by {M}. Assumes that {M} is
    orthonormal (i.e. distance-preserving). */

r2_vec_t cpk_map_points(r2_vec_t P, r2_t *o, r2x2_t *M);
  /* Maps each point of {P} by subtracting {o} from it,
    and multiplying the result (as a row vector) by
    the matrix {M}.  Infinite points are just copied. */
    
r2_vec_t cpk_unmap_points(r2_vec_t UrbM, r2x2_t *M, r2_t *o);
  /* Undoes the effect of {cpk_map_points(P,o,M)}. */

hxg_canvas_t *cpk_make_canvas_canon(double step, interval_t B[]);
  /* Returns a canvas {cvs} with an uniform hexagonal pixel array,
    with pixel spacing {step}, that covers the rectangle {B}. */
     
void cpk_set_urban_gridpoints_canon(hxg_canvas_t *cvs, bool_t valid[], r2_vec_t *R);
  /* Sets all pixels of {cvs} that lie inside the polygon {R} to TRUE.
    The polygon {R} must be closed and may have multiple islands
    and/or holes, serparated by infinite points. */
     
void cpk_erase_invalid_gridpoints_canon
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    cpk_domain_t *C, 
    cpk_policy_t *P
  );
  /* Eliminates every candidate to auction center in {cvs} that lies
    too close to an existing or planned station, to an intermunicipal
    border, or to an international border, as determined by {C} and {P}.
    Such pixels are set to zero; other pixels remain unchanged. */
  
void cpk_compute_URB_weights
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    float URB[], 
    double rAuc, 
    double rRad
  );
  /* Assumes {valid[k]} is true iff pixel {k} of {cvs} is within the
    urbanized area. For every pixel {k} such that {valid[k]} is true,
    computes the weight term {URB[k]} for that pixel's position. */

void cpk_compute_DEM_CLS_weights
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    uint32_t DEM[], 
    float CLS[], 
    r2_vec_t Dem, 
    double rAuc
  );
  /* For very pixel {k} such that {valid[k]} is TRUE, 
    computes the weight terms {DEM[k]} and {CLS[k]} for 
    that pixel's position. */

void cpk_proximity_graph
  ( hxg_canvas_t *cvs, /* Canvas. */
    bool_t valid[],     /* Which pixels to take. */
    double pixW[],      /* Weight of each pixel. */
    r2_vec_t *V,        /* OUT: Vertex coordinates. */
    double_vec_t *W,    /* OUT: Value of each vertex. */
    double dLim,        /* Distance threshold. */
    ui2_vec_t *E         /* OUT: Proximity edges. */
  );
  /* Takes a canvas {cvs} and two vectors {valid} and {pixW}, with one
    element per pixel, stored in the standard way. Returns an undirected
    graph {(V,E)} whose vertices are the positions of all pixels {k} of
    {cvs} with {valid[k] = TRUE}, with an edge {(i,j)} whenever {i<j}
    and {dist(V[i],V[j]) < dLim}. Also sets {W[i]} to the pixel's value
    {pixW[k]}. The vectors {V}, {W}, and {E} are allocated by the
    procedure. */

/* IMPLEMENTATIONS */

#define ND_NV_MAX 30000
/* Max number of vertices in no-demand mode */

void cpk_build_graph
  ( cpk_domain_t *C,      /* Valid points specs. */
    cpk_policy_t *P,      /* RadCom allocation parameters. */
    r2_vec_t *V,          /* OUT: List of candidate points. */
    double_vec_t *W,      /* OUT: Demand coverage of each candidate auction point. */
    ui2_vec_t *E,         /* OUT: Proximity (= incompatibility) edges. */
    bool_t verbose        /* TRUE prints diagnostic messages. */
  )
  {
    if (verbose) { cpk_trace_entry(); }
    
    /* Min dist between auction centers. */
    double dAAMin = P->dMin + 2*P->rAuc; 
    /* Min dist between auction center and any installed station. */
    double dASMin = P->dMin + P->rAuc;   
    if (verbose) 
      { 
        fprintf(stderr, "Derived parameters:\n");
        fprintf(stderr, "  min station-station distance = " XY_FFMT "\n", P->dMin);
        fprintf(stderr, "  min auction-station distance = " XY_FFMT "\n", dASMin);
        fprintf(stderr, "  min auction-auction distance = " XY_FFMT "\n", dAAMin);
        fprintf(stderr, "  min dist auction-municipalities = " XY_FFMT "\n", P->dMun);
        fprintf(stderr, "  min dist auction-nations = " XY_FFMT "\n", P->dNat);
      }

    /* Normalize and print weights: */
    cpk_normalize_weight_coeffs(&(P->wcoeffs));
    if (verbose)
      { fprintf(stderr, "Coefficients for weight terms\n");
        cpk_print_weights(stderr, &(P->wcoeffs));
      }
    
    /* Choose a reference point: */
    r2_t o = (r2_t){{0,0}};
    
    /* Choose the grid orientation. */
    double tilt = M_PI/7; /* An arbitrary angle. */
    r2_t dir = (r2_t){{cos(tilt),sin(tilt)}};
    if (verbose)
      { fprintf(stderr, "dir = (" XY_FFMT "," XY_FFMT ")\n", X(dir), Y(dir)); }

    /* Map the domain {C} so that {G->dir} is horizontal and {G->o} at (0,0): */
    r2x2_t M = cpk_get_rotation(&dir);
    cpk_domain_t CM = cpk_map_domain(C, &o, &M);
        
    /* Get the mapper domain's bounding box: */
    interval_t BM[2];
    r2_bbox(CM.Urb.ne, CM.Urb.e, BM, TRUE);
    if (verbose || (! verbose))
      { fprintf
          ( stderr, 
            "Mapped bounding box = [" XY_FFMT " _ " XY_FFMT "] × [" XY_FFMT " _ " XY_FFMT "]\n",
            LO(BM[0]), HI(BM[0]),   LO(BM[1]), HI(BM[1])
          );
      }        
    
    /* Choose the grid step: */
    double step;
    if (P->for_demand) 
      { /* Goal is about 12--18 candidates within {P->rAuc} from any point: */
        /* Guess: use {P->rAuc/2} minus a little bit. */
        step = 0.98 * P->rAuc/2;
      }
    else
      { /* Goal is to make a tight regular packing possible. */
        /* Compute area of bounding box: */
        double BMArea = (HI(BM[0])-LO(BM[0]))*(HI(BM[1])-LO(BM[1]));
        /* Find the minimum step so that the number of points is reasonable: */
        double min_step = sqrt(BMArea/ND_NV_MAX/(sqrt(3)/4));
        /* Use some exact divisor of {dAAMin}, plus a tiny slosh. */
        double n = floor(1.01*dAAMin/min_step);
        if (n < 1) { n = 1; }
        if (n > 8) { n = 8; }
        step = 1.01 * dAAMin/n;
      }
    if (verbose || (! verbose))
      { fprintf(stderr, "step = " XY_FFMT, step);
        fprintf(stderr, " ( = rAuc/%6.3f = dAAMin/%6.3f)\n", P->rAuc/step, dAAMin/step);
      }
    
    /* Get a canvas with all the points within the urbanized area set to 1: */
    hxg_canvas_t *cvs = cpk_make_canvas_canon(step, BM);
    if (verbose)
      { fprintf
          ( stderr, 
            "Allocated canvas o = (%f,%f) n = (%d,%d) s = (%f,%f)\n",
            X(cvs->org), Y(cvs->org), 
            cvs->size[0], cvs->size[1],
            X(cvs->step), Y(cvs->step)
          );
      }        
      
    uint32_t nPix = (uint32_t)(cvs->size[0]*cvs->size[1]);
    
    /* Set {valid[k] = TRUE} iff inside the urban area: */
    if (verbose) { fprintf(stderr, "Painting the polygon...\n"); }
    bool_vec_t valid = bool_vec_new(nPix);
    for (int32_t i = 0; i < nPix; i++) { valid.e[i] = FALSE; }
    cpk_set_urban_gridpoints_canon(cvs, valid.e, &(CM.Urb));
    
    /* Compute the {URB} term for all pixels within the urban area: */
    float_vec_t URB = float_vec_new(nPix);
    cpk_compute_URB_weights(cvs, valid.e, URB.e, P->rAuc, P->dMin/2);

    /* Invalidate all pixels too close to borders or other stations: */
    cpk_erase_invalid_gridpoints_canon(cvs, valid.e, &CM, P);
    
    /* Compute the {DEM} and {CLS} weight terms for all valid pixels: */
    uint32_vec_t DEM = uint32_vec_new(nPix);
    float_vec_t CLS = float_vec_new(nPix);
    cpk_compute_DEM_CLS_weights(cvs, valid.e, DEM.e, CLS.e, CM.Dem, P->rAuc);
    
    /* If in demand-driven mode, invalidate pixels that cover no DIs: */
    if (P->for_demand)
      { for (int32_t i = 0; i < nPix; i++) { if (DEM.e[i] == 0) { valid.e[i] = FALSE; } } }
    
    /* Compute pixel weights {pixW} for valid pixels: */
    if (verbose) { fprintf(stderr, "Computing vertex weights...\n"); }
    double_vec_t pixW = double_vec_new(nPix);
    { int32_t iV = 0; 
      for (int32_t i = 0; i < nPix; i++)
        { if (valid.e[i])
            { pixW.e[i] = cpk_compute_weight
                ( &(P->wcoeffs), 
                  (double)URB.e[i], 
                  (double)DEM.e[i], 
                  (double)CLS.e[i]
                );
              if (verbose) 
                { if ((iV < 20) || (iV >= nPix - 20))
                    { fprintf 
                        ( stderr, 
                          "  %7d  URB = %6.4f  DEM = %d  CLS = %6.4f\n", 
                          i, URB.e[i], DEM.e[i], CLS.e[i]
                        );
                    }
                  else if (iV == 20)
                    { fprintf(stderr, "  ...\n"); }
                }
              iV++;
            }
          else 
            { pixW.e[i] = 0.0; }
        }
    }
    
    /* Housecleaning: */
    free(URB.e);
    free(DEM.e);
    free(CLS.e);
    
    /* Now extract the positive vertices and edges: */
    r2_vec_t VM;
    cpk_proximity_graph(cvs, valid.e, pixW.e, &VM, W, dAAMin, E);

    /* Unmap the vertices: */
    (*V) = cpk_unmap_points(VM, &M, &o);
    
    /* Housecleaning: */
    cpk_domain_free(&CM);
    free(valid.e); 
    free(pixW.e);
    free(VM.e);
    
    if (verbose) { cpk_trace_exit(); }
  }
    
void cpk_proximity_graph
  ( hxg_canvas_t *cvs, /* Canvas. */
    bool_t valid[],    /* Which pixels to take. */
    double pixW[],     /* Weight of each pixel. */
    r2_vec_t *V,       /* OUT: Vertex coordinates. */
    double_vec_t *W,   /* OUT: Value of each vertex. */
    double dLim,       /* Distance threshold. */
    ui2_vec_t *E        /* OUT: Proximity edges. */
  )
  {
    uint32_t nx = cvs->size[0], ny = cvs->size[1]; 
    uint32_t nPix = nx*ny; uint32_t k; 
    /* Count nonzero pixels, build index {I}: */
    uint32_t nV = 0;
    int32_vec_t I = int32_vec_new(nPix); /* Linear pixel index {k} to vertex number {iV}. */
    { for (k = 0; k < nPix; k++) 
        { if (valid[k]) { I.e[k] = (int32_t)nV; nV++; } else { I.e[k] = -1; } }
    }
    /* fprintf(stderr, "Counted %d nonzero pixels.\n", nV); */

    /* Gather their indices and positions, in sweep order (Y, then X): */
    /* Also replace the pixel value by {i+1} where {i} is the index in {V}. */
    i2_vec_t K = i2_vec_new(nV); /* Vertex number {iV} to pixel indices {ix,iy}. */
    (*V) = r2_vec_new(nV);
    (*W) = double_vec_new(nV);
    { for (int32_t iy = 0; iy < ny; iy++)
        { uint32_t k = (uint32_t)iy*nx;
          for (int32_t ix = 0; ix < nx; ix++,k++)
            { if (valid[k]) 
                { i2_t u = (i2_t){{ix,iy}};
                  int32_t iV = I.e[k];
                  assert(iV >= 0);
                  K.e[iV] = u;
                  V->e[iV] = hxg_canvas_pixel_pos(cvs, u);
                  W->e[iV] = pixW[k];
                }
            }
        }
      /* fprintf(stderr, "Collected %d nonzero pixels.\n", nV); */
    }
    
    /* Get neighborhood templates {D[0],D[1]} for even and odd scanlines: */
    i2_vec_t D[2];
    { for (int32_t r = 0; r < 2; r++)
        { D[r] = hxg_canvas_half_disk_arcs(cvs, (i2_t){{0,r}}, dLim);
          /* fprintf(stderr, "Template for (iy%%2)=%d has %d arcs.\n", r, D[r].ne); */
        }
    }
    /* Get the proximity edges. */
    { /* Allocate the edge vector: */
      (*E) = ui2_vec_new(V->ne*D[0].ne); /* Estimated size. */
      uint32_t nE = 0;
      for (int32_t iV = 0; iV < nV; iV++)
        { /* Get pixel indices {ix,iy} of {V[iV]}: */
          i2_t *Ki = &(K.e[iV]);
          int32_t ix = X(*Ki), iy = Y(*Ki); 
          /* Get template edges appropriate for scanline {iy}: */
          i2_vec_t *Di = &(D[iy & 1]);
          uint32_t nD = Di->ne;
          /* Apply template to {ix,iy}: */
          for (int32_t s = 0; s < nD; s++)
            { /* Get edge {(ex,ey)} from template: */
              i2_t *e = &(Di->e[s]);
              int32_t ex = e->c[0], ey = e->c[1]; 
              /* Compute destination pixel {jx,jy}: */
              int32_t jx = ix + ex, jy = iy + ey;
              if ((jx >= 0) && (jx < nx) && (jy >= 0) && (jy < ny))
                { /* Get the pixel's linearized index {k}: */
                  uint32_t k = (uint32_t)jy*nx + (uint32_t)jx;
                  if (valid[k])
                    { /* Destination is valid, get its index {j} in {V}: */
                      int32_t jV = I.e[k];
                      /* Sanity checks: */
                      assert((iV >= 0) && (iV < nV));
                      assert((jV >= 0) && (jV < nV));
                      assert(jV > iV);
                      assert(X(K.e[jV]) == jx);
                      assert(Y(K.e[jV]) == jy);
                      /* Add edge to {E}: */
                      ui2_vec_expand(E, (vec_index_t)nE); /* Just to be sure... */
                      E->e[nE] = (ui2_t){{ (uint32_t)iV, (uint32_t)jV }};
                      nE++;
                    }
                }
            }
        }
      ui2_vec_trim(E, nE);
    }
    
    /* Housecleaning: */
    free(D[0].e);
    free(D[1].e);
    free(K.e);
    free(I.e);
  }



r2x2_t cpk_get_rotation(r2_t *dir)
  { 
    /* Get the grid's scanline direction: */
    r2_t u = *dir;
    
    /* Normalize {u} to unit length: */
    r2_dir(&u, &u);
    
    /* Get the direction {v} orthogonal to {u}, ccw from it: */
    double ux = X(u), uy = Y(u);
    double vx = -uy, vy = ux;
    return (r2x2_t){{{ux,vx},{uy,vy}}};
  }

cpk_domain_t cpk_map_domain(cpk_domain_t *C, r2_t *o, r2x2_t *M)
  { 
    cpk_domain_t CM;
    /* Rewrite the polygon and the fixed sites in the {o,u,v} coord system: */
    CM.Urb = cpk_map_points(C->Urb, o, M);
    CM.Exs = cpk_map_points(C->Exs, o, M);
    CM.Auc = cpk_map_points(C->Auc, o, M);
    CM.Mun = cpk_map_points(C->Mun, o, M);
    CM.Nat = cpk_map_points(C->Nat, o, M);
    CM.Dem = cpk_map_points(C->Dem, o, M);
    return CM;
  }

r2_vec_t cpk_map_points(r2_vec_t P, r2_t *o, r2x2_t *M) 
  { r2_vec_t Q = r2_vec_new(P.ne);
    for (int32_t k = 0; k < P.ne; k++) 
      { r2_t *Pk = &(P.e[k]);
        r2_t *Qk = &(Q.e[k]);
        if (r2_is_finite(Pk))
          { r2_t d; 
            r2_sub(Pk, o, &d); 
            r2x2_map_row(&d, M, Qk); 
          }
        else
          { *Qk = *Pk; }
      }
    return Q;
  }    
  
r2_vec_t cpk_unmap_points(r2_vec_t P, r2x2_t *M, r2_t *o) 
  { r2_vec_t Q = r2_vec_new(P.ne);
    for (int32_t k = 0; k < P.ne; k++) 
      { r2_t *Pk = &(P.e[k]);
        r2_t *Qk = &(Q.e[k]);
        if (r2_is_finite(Pk))
          { r2_t d; 
            /* The inverse of a rotation matrix is its transpose: */ 
            r2x2_map_col(M, Pk, &d); 
            r2_add(o, &d, Qk);
          }
        else
          { *Qk = *Pk; }
      }
    return Q;
  }    
  
hxg_canvas_t *cpk_make_canvas_canon(double step, interval_t B[])
  {     
    demand(step >= 1.0e-120, "invalid {step}");
    demand(LO(B[0]) < HI(B[0]), "invalid empty box");
    demand(LO(B[1]) < HI(B[1]), "invalid empty box");
    
    /* Grid steps in {x} and {y}: */
    double sx = step; 
    double sy = step*sqrt(3)/2;
        
    /* Round the corners to grid points ignoring odd-even shifts: */
    int32_t ixlo = 2*(int32_t)floor(LO(B[0])/sx/2 - 1.0e-6);
    int32_t ixhi = 2*(int32_t)ceil(HI(B[0])/sx/2 + 1.0e-6);
    int32_t iylo = 2*(int32_t)floor(LO(B[1])/sy/2 - 1.0e-6);
    int32_t iyhi = 2*(int32_t)ceil(HI(B[1])/sy/2 + 1.0e-6);
    
    assert(ixhi > ixlo);
    assert(iyhi > iylo);
        
    /* Compute the grid dimensions: */
    uint32_t nx = (uint32_t)(ixhi - ixlo); assert(nx >= 2);
    uint32_t ny = (uint32_t)(iyhi - iylo); assert(ny >= 2);
        
    /* Allocate the grid: */
    hxg_canvas_t *cvs = hxg_canvas_new(nx, ny, step, FALSE);
    cvs->org.c[0] += ixlo*cvs->step.c[0];
    cvs->org.c[1] += iylo*cvs->step.c[0];
    return cvs;
  }

void cpk_set_urban_gridpoints_canon(hxg_canvas_t *cvs, bool_t valid[], r2_vec_t *R)
  {
    /* Paint the polygon {R} over the grid: */
    auto void set_true(uint32_t ix, uint32_t iy, uint32_t k);
    
    hxg_paint_polygon(cvs, R, set_true);

    void set_true(uint32_t ix, uint32_t iy, uint32_t k) { valid[k] = TRUE; }
  }
 
void cpk_erase_invalid_gridpoints_canon
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    cpk_domain_t *C, 
    cpk_policy_t *P
  )
  {
    /* Erase all pixels near existing or planned stations: */
    /* fprintf(stderr, "Erasing areas near existing/planned stations...\n"); */
    auto void set_false(uint32_t ix, uint32_t iy, uint32_t k);

    uint32_t nP =  C->Exs.ne;
    r2_t *c = C->Exs.e;
    for (int32_t i = 0; i < nP; i++)
      { hxg_paint_circle(cvs, &(c[i]), P->dMin + P->rAuc, set_false); }
        
    /* fprintf(stderr, "Erasing areas near predetermined auction points...\n"); */
    auto void set_false(uint32_t ix, uint32_t iy, uint32_t k);
      { uint32_t nP =  C->Auc.ne;
        r2_t *c = C->Auc.e;
        for (int32_t i = 0; i < nP; i++)
            { hxg_paint_circle(cvs, &(c[i]), P->rAuc + P->dMin + P->rAuc, set_false); }
      }
        
    /* Erase all pixels near intermunicipal or international borders: */
    /* fprintf(stderr, "Erasing border zones...\n"); */
    hxg_paint_sausage(cvs, &(C->Mun), P->dMun, set_false);
    hxg_paint_sausage(cvs, &(C->Nat), P->dNat, set_false);

    void set_false(uint32_t ix, uint32_t iy, uint32_t k)
      { valid[k] = FALSE; }
  }

void cpk_compute_URB_weights
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    float URB[], 
    double rAuc, 
    double rRad
  )
  {
    /* For each valid (i.e. nonzero) pixel at a point {p}, paints a
    disk with center {p} and radius {rAuc+rRad}. To each nonzero pixel
    {q} within this circle, adds {f(dist(p,q))/fsum}, where {f(r)} is 1 for
    {r < rRad-rAuc}, 0 for {r > rRad+rAuc}, and some smooth
    interpolator between these two values; and {fsum} is a nomrmalization
    factor so that the maximum final value is 1. */

    uint32_t nx = cvs->size[0], ny = cvs->size[1];
    uint32_t nPix = nx*ny; 
    for (int32_t k = 0; k < nPix; k++) { URB[k] = 0.0; }
    double dMax = rRad + rAuc;
    double dMin = fabs(rRad - rAuc);
    double dW = dMax - dMin;
    
    auto double cprob(int32_t ix, int32_t iy, r2_t *p);
      /* Approx probability that a disk of radius {rRad} centered
        at a random point within the auction disk of pixel {ix,iy}
        will cover the point {p}. */
    
    /* Compute the normalization factor {fsum}: */
    r2_t *o = &(cvs->org);
    double fsum = hxg_canvas_sum_over_circle(cvs, *o, dMax, cprob, o);

    /* Now distribute the audience over the auction points: */
    for (int32_t iy = 0; iy < ny; iy++)
      { uint32_t k = (uint32_t)iy*nx; 
        for (int32_t ix = 0; ix < nx; ix++,k++)
          { if (valid[k]) 
              { /* Pixel {k} is inside the urban area -- a potential listener. */
                r2_t pM = hxg_canvas_pixel_pos(cvs, (i2_t){{ix,iy}}); /* Listene's position. */

                /* Count it as potential audience for by nearby auction points: */
                auto void add_audience(uint32_t jx, uint32_t jy, uint32_t r);
                hxg_paint_circle(cvs, &pM, dMax, add_audience);
                
                void add_audience(uint32_t jx, uint32_t jy, uint32_t r)
                  { if ((valid[r]) && (jx >= 0) && (jx < nx) && (jy >= 0) && (jy < ny)) 
                      { /* Get prob {f} that an auction at {jx,jy} will serve {pM}: */
                        double f = cprob((int32_t)jx, (int32_t)jy, &pM);
                        /* Add to pixel, normalized: */
                        URB[r] += (float)(f/fsum);
                      }
                  }
              }
          }
      }


    double cprob(int32_t ix, int32_t iy, r2_t *p)
      { /* Get the pixel's coordinates: */
        r2_t q = hxg_canvas_pixel_pos(cvs, (i2_t){{ ix, iy }});
        /* Compute distance from potential audience {p}: */
        double d = r2_dist(&q, p);
        /* Approximate prob that an auction at {q} will serve {p}: */
        double f;
        if (d <= dMin)
          { f = 1; }
        else if (d >= dMax)
          { f = 0; }
        else 
          { double s = (d-dMin)/dW;
            f = 1 - (3 - 2*s)*s*s;
          }
        return f;
      }
  }
 
void cpk_compute_DEM_CLS_weights
  ( hxg_canvas_t *cvs, 
    bool_t valid[], 
    uint32_t DEM[], 
    float CLS[], 
    r2_vec_t Dem, 
    double rAuc
  )
  {        
    uint32_t nx = cvs->size[0], ny = cvs->size[1];
    uint32_t nPix = nx*ny; 
    for (int32_t k = 0; k < nPix; k++) { DEM[k] = 0.0; CLS[k] = 0; }
    
    uint32_t nA = Dem.ne;
    double r2 = rAuc*rAuc;
    for (int32_t i = 0; i < nA; i++)
      { r2_t *Ck = &(Dem.e[i]);
        /* Increment all valid pixels inside the circle: */

        auto void bump_valid(uint32_t ix, uint32_t iy, uint32_t k);

        hxg_paint_circle(cvs, Ck, rAuc, bump_valid);

        void bump_valid(uint32_t ix, uint32_t iy, uint32_t k)
          { if (valid[k])
              { /* Count one more eligible DI for this auction point: */
                DEM[k]++;
                /* Compute the closeness function {h}: */
                r2_t pM = hxg_canvas_pixel_pos(cvs, (i2_t){{ (int32_t)ix, (int32_t)iy }});
                double d2 = r2_dist_sqr(&pM, Ck);
                double h = (d2 >= r2 ? 0 : 1 - d2/r2);
                /* Add it to the {CLS} score of this auction point: */
                CLS[k] += (float)h;
              }
          }
      }
  }

