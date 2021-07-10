/* See mroc_stl.h */
/* Last edited on 2021-07-08 00:35:48 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <ppv_types.h>
#include <ppv_array.h>
#include <r3.h>
#include <jsmath.h>
#include <affirm.h>

#include <mroc.h>
#include <mroc_stl.h>

/* INTERNAL PROTOTYPES */
    

/* IMPLEMENTATIONS */

void mroc_stl_classify_tetra_vertices(double f[], int32_t *niP, int32_t ki[], int32_t *noP, int32_t ko[])
  {
    int32_t ni = 0; /* Number of corners inside the object. */
    int32_t no = 0; /* Number of corners outside the object. */
    int32_t k;
    for (k = 0; k < 4; k++)
      { assert(f[k] != 0.5); 
        if (f[k] > 0.5) 
          { ki[ni] = k; ni++; }
        else if (f[k] < 0.5) 
          { ko[no] = k; no++; }
        else
          { affirm(FALSE, "tetrahedron corner is ambiguous"); }
      }
    (*niP) = ni;
    (*noP) = no;
  }

void mroc_stl_write_surface_in_tetra
  ( FILE *wr, 
    r3_t p[], 
    double f[], 
    float eps,
    int32_t *ntP,
    int32_t *neP
  )
  {
    /* Check for intersection: */
    int32_t ni; /* Number of corners inside the object. */
    int32_t ki[4]; /* Indices of inside corners are {ki[0..ni-1]}. */
    int32_t no; /* Number of corners outside the object. */
    int32_t ko[4]; /* Indices of outside corners are {ko[0..no-1]}. */
    
    mroc_stl_classify_tetra_vertices(f, &ni, ki, &no, ko);
    assert(ni + no == 4);
    
    if ((ni == 0) || (ni == 4)) { return; }
    /* Tetrahedron intersects surface: */
    int32_t nq = 0;   /* Number of edges that cross the surface. */
    i3_t q[4];        /* Points where edges cross surface are {q[0..nq-1]} in cyclic order, quantized. */
    if (ni == 1)
      { /* Intersection is a single triangle: */ 
        assert(no == 3);
        int32_t k;
        for (k = 0; k < 3; k++)
          { q[nq] = mroc_stl_edge_crossing(&(p[ki[0]]), f[ki[0]], &(p[ko[k]]), f[ko[k]], eps); nq++; }
        assert(nq == 3);
        mroc_stl_write_i3_triangle(wr, &(q[0]), &(q[1]), &(q[2]), eps, ntP, neP);
      }
     else if (ni == 3)
      { /* Intersection is a single triangle: */ 
        assert(no == 1);
        int32_t k;
        for (k = 0; k < 3; k++)
          { q[nq] = mroc_stl_edge_crossing(&(p[ko[0]]), f[ko[0]], &(p[ki[k]]), f[ki[k]], eps); nq++; }
        assert(nq == 3);
        mroc_stl_write_i3_triangle(wr, &(q[0]), &(q[1]), &(q[2]), eps, ntP, neP);
      }
    else if (ni == 2)
      { /* Intersection is a quadrilateral. */
        /* Compute the 4 corners in cyclic order: */
        q[nq] = mroc_stl_edge_crossing(&(p[ki[0]]), f[ki[0]], &(p[ko[0]]), f[ko[0]], eps); nq++;
        q[nq] = mroc_stl_edge_crossing(&(p[ki[0]]), f[ki[0]], &(p[ko[1]]), f[ko[1]], eps); nq++;
        q[nq] = mroc_stl_edge_crossing(&(p[ki[1]]), f[ki[1]], &(p[ko[1]]), f[ko[1]], eps); nq++;
        q[nq] = mroc_stl_edge_crossing(&(p[ki[1]]), f[ki[1]], &(p[ko[0]]), f[ko[0]], eps); nq++;
        assert(nq == 4);
        /* Choose diagonal to split: */
        if (i3_dist_sqr(&(q[0]), &(q[2])) < i3_dist_sqr(&(q[1]), &(q[3])))
          { /* Split by 0-2: */
            mroc_stl_write_i3_triangle(wr, &(q[0]), &(q[1]), &(q[2]), eps, ntP, neP);
            mroc_stl_write_i3_triangle(wr, &(q[0]), &(q[2]), &(q[3]), eps, ntP, neP);
          }
        else
          { /* Split by 1-3: */
            mroc_stl_write_i3_triangle(wr, &(q[0]), &(q[1]), &(q[3]), eps, ntP, neP);
            mroc_stl_write_i3_triangle(wr, &(q[2]), &(q[1]), &(q[3]), eps, ntP, neP);
          }
      }
    else
      { fatalerror("invalid {ni,no}"); }
  }
 
void mroc_stl_write_tetra_faces
  ( FILE *wr, 
    r3_t p[], 
    double f[], 
    int32_t *ntP
  )
  {
    /* Check for intersection: */
    int32_t ni; /* Number of corners inside the object. */
    int32_t ki[4]; /* Indices of inside corners are {ki[0..ni-1]}. */
    int32_t no; /* Number of corners outside the object. */
    int32_t ko[4]; /* Indices of outside corners are {ko[0..no-1]}. */
    
    mroc_stl_classify_tetra_vertices(f, &ni, ki, &no, ko);
    assert(ni + no == 4);

    if ((ni == 0) || (ni == 4)) { return; }
    /* Tetrahedron intersects surface: */

    /* Compute barycenter: */
    r3_t ctr;        /* Barycenter of tetrahedron. */
    r3_t p01, p23;
    r3_mix(0.5, &(p[0]), 0.5, &(p[1]), &p01);
    r3_mix(0.5, &(p[2]), 0.5, &(p[3]), &p23);
    r3_mix(0.5, &p01, 0.5, &p23, &ctr);

    /* Srink the tetrahedron towards the center: */
    r3_t q[4];       /* Corners of shrunk tetrahedron. */
    int32_t k;
    double vfac = t2s_TETRA_SHRINK;
    double cfac = 1 - vfac; 
    for (k = 0; k < 4; k++) { r3_mix(cfac, &ctr, vfac, &(p[k]), &(q[k]));  }
    
    /* Write the faces of the srunk tetrahedron: */ 
    for (k = 0; k < 4; k++)
      { int32_t k1 = (k+1) % 4;
        int32_t k2 = (k+2) % 4;
        int32_t k3 = (k+3) % 4;
        mroc_stl_write_r3_triangle(wr, &(q[k1]), &(q[k2]), &(q[k3])); (*ntP)++;
      }
  }
    
#define mroc_stl_MAX_QUANTIZED (10000000)
  /* Safety parameter: max abs value of quantized coordinate, to avoid overflows. */

i3_t mroc_stl_edge_crossing(r3_t *p0, double f0, r3_t *p1, double f1, float eps)
  {
    /* Compute the relative position {r} of crossing point along segment: */
    double df = f1 - f0; 
    demand(df != 0, "same value");
    double r = (0.5 - f0)/df;
    demand ((r > 0) && (r < 1), "does not cross 0.5");

    i3_t q; /* Crossing point. quantized. */
    int32_t j;
    for (j = 0; j < 3; j++)
      { /* Compute crossing point coordinate {qj}: */
        float qj = (float)((1-r)*p0->c[j] + r*p1->c[j]);
        /* Round it to even integer multiple of {eps}: */
        q.c[j] = (int32_t)iroundfrac(qj, eps, 2, 0, INT32_MAX); 
        if (abs(q.c[j]) > mroc_stl_MAX_QUANTIZED)
          { fprintf(stderr, "** quantized coordinate %+11.7f --> %d too big; use larger {eps}\n", qj, q.c[j]); 
            assert(FALSE);
          }
      }
    return q;
  }

void mroc_stl_write_i3_triangle
  ( FILE *wr, 
    i3_t *p0, 
    i3_t *p1, 
    i3_t *p2,
    float eps,
    int32_t *ntP,
    int32_t *neP
  )
  {
    /* Check for repeated vertices: */
    if (i3_eq(p0, p1) || i3_eq(p1, p2) || i3_eq(p2, p0)) { (*neP)++; return; }
    
    /* Unround and write out: */
    r3_t pf0 = mroc_stl_unround_point(p0, eps);
    r3_t pf1 = mroc_stl_unround_point(p1, eps);
    r3_t pf2 = mroc_stl_unround_point(p2, eps);
    mroc_stl_write_r3_triangle(wr, &pf0, &pf1, &pf2); (*ntP)++;      
  }
     
r3_t mroc_stl_unround_point(i3_t *p, float eps)
  { 
    r3_t r;
    r.c[0] = ((double)eps)*((double)p->c[0]);
    r.c[1] = ((double)eps)*((double)p->c[1]);
    r.c[2] = ((double)eps)*((double)p->c[2]);
    return r;
  }

void mroc_stl_write_r3_triangle
  ( FILE *wr, 
    r3_t *p0, 
    r3_t *p1, 
    r3_t *p2
  )
  {
    /* Compute normal: */
    r3_t d = mroc_stl_compute_normal(p0, p1, p2);
    
    /* Write face: */
    fprintf(wr, "\n");
    fprintf(wr, "facet\n");
    r3_gen_print(wr, &d, "%+7.4f", "normal ", " ", "\n");
    fprintf(wr, "outer loop\n");
    r3_gen_print(wr, p0, "%8.3f", "  vertex ", " ", "\n");
    r3_gen_print(wr, p1, "%8.3f", "  vertex ", " ", "\n");
    r3_gen_print(wr, p2, "%8.3f", "  vertex ", " ", "\n");
    fprintf(wr, "endloop\n");
    fprintf(wr, "endfacet\n");
  }
  
r3_t mroc_stl_compute_normal(r3_t *p0, r3_t *p1, r3_t *p2)
  {
    r3_t a, b, d;
    r3_sub(p1, p0, &a);
    r3_sub(p2, p0, &b);
    r3_cross(&a, &b, &d);
    double len = r3_dir(&d, &d);
    if (len < 1.0e-6)
      { fprintf(stderr, "!! degenerate triangle -- defaulting normal\n");
        r3_gen_print(stderr, p0, "%.4f", "  p0 = ( ", " ", " )\n");
        r3_gen_print(stderr, p1, "%.4f", "  p1 = ( ", " ", " )\n");
        r3_gen_print(stderr, p2, "%.4f", "  p2 = ( ", " ", " )\n");
        r3_gen_print(stderr, &a, "%.4f", "  a  = ( ", " ", " )\n");
        r3_gen_print(stderr, &b, "%.4f", "  b  = ( ", " ", " )\n");
        d = (r3_t){{ 0.0, 0.0, 1.0 }};
      }
    return d;
  }
