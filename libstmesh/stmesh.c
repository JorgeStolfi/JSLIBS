/* See {stmesh.h} */
/* Last edited on 2016-04-21 18:56:18 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <i3.h>
#include <r3.h>
#include <bvtable.h>
#include <bvhash.h>

/* #include <stmesh_STL.h> */
  
#include <stmesh.h>

void stmesh_print_edge_degrees(FILE *wr, stmesh_t mesh)
  {
    uint32_t ne = stmesh_edge_count(mesh);
    uint32_t nf = stmesh_face_count(mesh);
    
    uint32_t deg_max = nf;
    uint32_t d;
    
    /* Scan edges and count edges of the same degree: */
    uint32_t *ne_deg = notnull(malloc((deg_max+1)*sizeof(uint32_t)), "no mem"); 
    for (d = 0; d <= deg_max; d++) { ne_deg[d] = 0; }
    stmesh_edge_unx_t uxe;
    for (uxe = 0; uxe < ne; uxe++)
      { stmesh_edge_t e = stmesh_get_edge(mesh, uxe);
        d = stmesh_edge_degree(e);
        if (d == 0) 
          { fprintf(stderr, "** edge %u has zero degree\n", uxe);
            assert(d > 0);
          }
        assert(d <= deg_max);
        ne_deg[d]++;
      }
    
    /* Print degree statistics: */
    bool_t border = FALSE;
    int32_t ne_tot = 0;
    int32_t nf_tot = 0;
    for (d = 0; d <= deg_max; d++)
      { if (ne_deg[d] > 0)
          { fprintf(wr, "there are %u edges shared by %u faces\n", ne_deg[d], d);
            ne_tot += ne_deg[d];
            nf_tot += ne_deg[d]*d;
            if ((d % 2) != 0) { border = TRUE; }
          }
      }
    assert((nf_tot % 3) == 0);
    nf_tot /= 3;
    assert(ne_tot == ne);
    assert(nf_tot == nf);
    
    /* if there are edges with odd degree, the mesg is not a closed surface: */
    if (border) { fprintf(stderr, "!! warning: mesh is not closed\n"); }
    
    free(ne_deg);
  }
     
r3_t stmesh_unround_point(i3_t *p, float eps)
  { 
    r3_t r;
    r.c[0] = ((double)eps)*((double)p->c[0]);
    r.c[1] = ((double)eps)*((double)p->c[1]);
    r.c[2] = ((double)eps)*((double)p->c[2]);
    return r;
  }

void stmesh_print_bounding_box(FILE *wr, stmesh_t mesh)
  { 
    float eps = stmesh_get_eps(mesh);
    i3_t minQ, maxQ;
    stmesh_get_bounding_box(mesh, &minQ, &maxQ);
    fprintf(wr, "bounding box:\n");
    int k;
    for (k = 0; k < 3; k++)
      { int32_t minQk = minQ.c[k]; float minFk = eps*(float)minQk;
        int32_t maxQk = maxQ.c[k]; float maxFk = eps*(float)maxQk;
        fprintf(wr, "  %c [ %d _ %d ] = [ %.8f mm _ %.8f mm ]\n", "XYZ"[k], minQk, maxQk, minFk, maxFk); 
      }
  }
  
bool_t stmesh_edge_crosses_plane(stmesh_edge_t e, int32_t pZ)
  {
    bool_t verbose = FALSE;
    stmesh_vert_t v[2]; /* Endpoints of {e}. */
    stmesh_edge_get_endpoints(e, v);
    int32_t vZ[2]; /* Quantized {Z}-coords of endpoints of {e}. */
    int r;
    for (r = 0; r < 2; r++) 
      { i3_t p = stmesh_vert_get_pos(v[r]); vZ[r] = p.c[2];
        if (verbose) { fprintf(stderr, "endpoint %d at Z = %+d\n", r, vZ[r]); }
      }
    assert((vZ[0] != pZ) && (vZ[1] != pZ));
    if ((vZ[0] < pZ) && (vZ[1] > pZ)) { return TRUE; }
    if ((vZ[0] > pZ) && (vZ[1] < pZ)) { return TRUE; }
    return FALSE;
  }

r2_t stmesh_edge_plane_intersection(stmesh_edge_t e, int32_t pZ, double eps)
  {
    stmesh_vert_t v[2]; /* Endpoints of {e}. */
    stmesh_edge_get_endpoints(e, v);
    i3_t p[2]; /* Quantized coordinates of the endpoints. */
    int r;
    for (r = 0; r < 2; r++) { p[r] = stmesh_vert_get_pos(v[r]); }
    int32_t nZ = pZ - p[0].c[2];
    int32_t dZ = p[1].c[2] - p[0].c[2];
    assert(dZ != 0);
    double t = ((double)nZ)/((double)dZ);
    assert((t > 0.0) && (t < 1.0));
    r2_t u;
    int k;
    for (k = 0; k < 2; k++)
      { int32_t dk = p[1].c[k] - p[0].c[k];
        u.c[k] = eps*(p[0].c[k] + t*dk);
      }
    return u;
  }

void stmesh_face_get_sliced_sides(stmesh_t mesh, stmesh_face_t f, int32_t pZ, stmesh_edge_t e[])
  { 
    bool_t verbose = FALSE;
    
    stmesh_edge_t side[3]; /* The sides of {f} */
    stmesh_face_get_sides(f, side);
    if (verbose) { fprintf(stderr, "mesh face %u and plane at Z = %+d\n", stmesh_face_get_unx(mesh, f), pZ); }
    int k;
    int m = 0;
    for (k = 0; k < 3; k++)
      { stmesh_edge_t ek = side[k];
        if (stmesh_edge_crosses_plane(ek, pZ))
          { if (verbose) 
              { stmesh_edge_unx_t uxek = stmesh_edge_get_unx(mesh, ek);
                fprintf(stderr, "  side %d = edge %u crosses\n", k, uxek);
              }
            assert(m < 2);
            e[m] = stmesh_edge_natural(ek);
            m++;
          }
      }
    demand(m == 2, "face does not cross plane");
  }

stmesh_t stmesh_build
  ( float eps, 
    uint32_t nv, 
    i3_t vpos[], 
    uint32_t ne, 
    stmesh_vert_unx_pair_t endv[],
    uint32_t nf, 
    stmesh_edge_unx_triple_t side[],
    bool_t checkSorted
  )
  {  
    bool_t verify = TRUE;
    
    demand(nv < stmesh_nv_MAX, "too many vertices");
    demand(ne < stmesh_ne_MAX, "too many edges");
    demand(nf < stmesh_nf_MAX, "too many triangles");
   
    stmesh_t mesh = stmesh_new_desc(eps, nv,ne,nf);

    /* Build the vertex table: */
    stmesh_vert_unx_t uxv;
    for (uxv = 0; uxv < nv; uxv++)
      { i3_t *vp = &(vpos[uxv]);
        stmesh_vert_unx_t uxv2 = stmesh_add_vert(mesh, vp);
        assert(uxv2 == uxv);
      }

    /* Build the edge table: */
    stmesh_edge_unx_t uxe;
    for (uxe = 0; uxe < ne; uxe++)
      { /* Get the indices of the endpoints from the {endv} list: */
        stmesh_vert_unx_pair_t *ev = &(endv[uxe]);
        stmesh_edge_unx_t uxe2 = stmesh_add_edge(mesh, &(ev->c[0]));
        assert(uxe2 == uxe);
      }
        
    /* Fill the face table: */
    stmesh_face_unx_t uxf;
    for (uxf = 0; uxf < nf; uxf++)
      { /* Get the indices of its (unoriented) side edges, from the {side} parameter: */
        stmesh_edge_unx_triple_t *fe = &(side[uxf]);
        stmesh_face_unx_t uxf2 = stmesh_add_face(mesh, &(fe->c[0]));
        assert(uxf2 == uxf);
      }

    stmesh_print_bounding_box(stderr, mesh);
    stmesh_print_edge_degrees(stderr, mesh);
    
    /* Paranoia: */
    if (verify) { stmesh_check(mesh); }
    
    return mesh;
  }
    
void stmesh_check(stmesh_t mesh)
  {
    
    /* Check vertices */
    uint32_t nv = stmesh_vert_count(mesh);
    stmesh_vert_unx_t uxv;
    for (uxv = 0; uxv < nv; uxv++)
      { stmesh_vert_t v = stmesh_get_vert(mesh, uxv);
        stmesh_vert_check(mesh, v);
      }

    /* Check edges */
    uint32_t ne = stmesh_edge_count(mesh);
    stmesh_edge_unx_t uxe;
    for (uxe = 0; uxe < ne; uxe++)
      { stmesh_edge_t e = stmesh_get_edge(mesh, uxe);
        stmesh_edge_check(mesh, e);
      }

    uint32_t nf = stmesh_face_count(mesh);
    stmesh_face_unx_t uxf;
    for (uxf = 0; uxf < nf; uxf++)
      { stmesh_face_t f = stmesh_get_face(mesh, uxf);
        stmesh_face_check(mesh, f);
      }
  }

void stmesh_vert_check(stmesh_t mesh, stmesh_vert_t v)
  {
    /* Check the index: */
    stmesh_vert_unx_t uxv = stmesh_vert_get_unx(mesh, v);
    assert(v == stmesh_get_vert(mesh, uxv));
  }

void stmesh_edge_check(stmesh_t mesh, stmesh_edge_t e)
  {
    int k;

    /* Check the index: */
    stmesh_edge_unx_t uxe = stmesh_edge_get_unx(mesh, e);
    stmesh_edge_t e0 = stmesh_get_edge(mesh, uxe); /* Natural vesion of {e} */
    
    /* Get the endpoints of {e}: */
    stmesh_vert_t ve[2];
    stmesh_edge_get_endpoints(e, ve);
    
    /* Compare with reversed versions: */
    bool_t found = FALSE; /* Found a version of {e} that matches {e0}. */
    stmesh_edge_t c = e;
    int t;
    for (t = 0; t < 2; t++)
      { /* At this point, {c = reverse^t(e)}. */
        stmesh_edge_t c1 = stmesh_edge_reverse(e, t);
        assert(c1 == c);
        
        /* Compare with natural version: */
        if (c == e0) { found = TRUE; }

        /* Check the endpoints of {c}: */
        stmesh_vert_t vc[2];
        stmesh_edge_get_endpoints(c, vc);
        for (k = 0; k < 2; k++)
          { 
            /* Compare endpoint  {vc[k]} with the endpoint of {e}: */
            if (t == 0)
              { assert(vc[k] == ve[k]); }
            else
              { assert(vc[k] == ve[1 - k]); }
            
            /* Check single endpoint: */
            stmesh_vert_t uk = stmesh_edge_get_endpoint(c, k);
            assert(vc[k] == uk);
          }
        c = stmesh_edge_reverse(c, 1);
      }
    assert(c == e);
    assert(found);
  }
   
void stmesh_face_check(stmesh_t mesh, stmesh_face_t f)
  {
    int k;

    /* Check the index: */
    stmesh_face_unx_t uxf = stmesh_face_get_unx(mesh, f);
    stmesh_face_t f0 = stmesh_get_face(mesh, uxf); /* Natural vesion of {f} */

    /* Get the sides and corners of {f}: */
    stmesh_edge_t ef[3];
    stmesh_face_get_sides(f, ef);

    stmesh_vert_t vf[3];
    stmesh_face_get_corners(f, vf);
    
    /* Compare with flipped and shifted versions: */
    bool_t found = FALSE; /* Found a version of {f} that matches {f0}. */
    stmesh_face_t g = f;
    int s, t;
    for (s = 0; s < 3; s++)
      { /* At this point, {g = shift^s(f)}. */
        stmesh_face_t g1 = stmesh_face_shift(f, s);
        assert(g1 == g);
        
        for (t = 0; t < 2; t++)
          { /* At this point, {g = flip^t(shift^s(f))}. */
            stmesh_face_t g2 = stmesh_face_flip(g1, t);
            assert(g2 == g);
            
            /* Compare with natural version: */
            if (g == f0) { found = TRUE; }
        
            /* Get the sides and corners of {g}: */
            stmesh_edge_t eg[3];
            stmesh_face_get_sides(g, eg);

            stmesh_vert_t vg[3];
            stmesh_face_get_corners(g, vg);
            
            /* Compare the edges of {f,g}: */
            for (k = 0; k < 3; k++)
              { if (t == 0)
                  { assert(eg[k] == ef[(k+s)%3]); }
                else
                  { assert(eg[k] == stmesh_edge_reverse(ef[(3-k+s)%3], 1)); }
                stmesh_edge_check(mesh, eg[k]);
              }
            
            /* Compare the corners of {f,g}: */
            for (k = 0; k < 3; k++)
              { if (t == 0)
                  { assert(vg[k] == vf[(k+s)%3]); }
                else
                  { assert(vg[k] == vf[(3-k+s)%3]); }
              }
            g = stmesh_face_flip(g, 1);
          }
        g = stmesh_face_shift(g, 1);
      }
    assert(g == f);
    assert(found);

    /* Check consistency of elements of {f}: */
    for (k = 0; k < 3; k++)
      { int k1 = (k+1)%3;
        int k2 = (k+2)%3;
        stmesh_vert_t uk = stmesh_edge_get_endpoint(ef[k1], 1);
        stmesh_vert_t wk = stmesh_edge_get_endpoint(ef[k2], 0);
        assert(vf[k] == uk);
        assert(vf[k] == wk);
      }

  }
   
