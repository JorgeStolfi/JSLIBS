/* See {stpoly.h} */
/* Last edited on 2016-04-21 18:55:47 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <i2.h>
#include <r2.h>
#include <bvtable.h>
#include <bvhash.h>

#include <stpoly_STP.h>
  
#include <stpoly.h>

/* INTERNAL PROTOTYPES */

typedef struct stpoly_vert_unx_pair_t { stpoly_vert_unx_t v[2]; } stpoly_vert_unx_pair_t;
  /* A pair of vertex indices. */

stpoly_t stpoly_build
  ( float eps, 
    uint32_t nv, 
    i2_t vpos[], 
    uint32_t ne, 
    stpoly_vert_unx_pair_t endv[]
  );
  /* Builds a polygonal figure data structure with fundamental length unit {eps},
    {nv} vertices, and {ne} unoriented edges,
    given the vertex coordinates and the basic topological information.
    
    The quantized coordinates of the vertex with index {uxv} will be {vpos[uxv]}. 
    
    The endpoints of the unoriented edge with index {uxe} will be
    the vertices with indices {endv[uxe].c[0..1]}.  
    
    The procedure also computes derived informetion
    such as vertex degrees, the fields {.minY,.maxY} of each edge,
    and the bounding box of the polygonal figure. */

/* IMPLEMENTATIONS */
  
stpoly_t stpoly_read_STP(char *fileName, bool_t binary, float eps, uint32_t neGuess, bool_t even)
  {
    bool_t debug = TRUE;
    
    char *format = ((char *[2]){ "ascii", "binary" })[binary];
    fprintf(stderr, "reading polygonal figure from file %s (%s)\n", fileName, format);
    fprintf(stderr, "expecting about %u edges\n", neGuess);
    fprintf(stderr, "quantizing vertex coords to%s multiples of %.8f mm\n", (even ? " even" : ""), eps);

    /* While reading the STP file and building the structure, the
      entries must be linked by *indices* of elements, rather
      than pointers, because the hash tables may be reallocated. Only
      after the reading is finished can we translate the indices into
      pointers.
      
      During the reading, a vertex is represented by its
      integer position (as an {i2.t}); and an (unoriented) edge is represented
      by indices of the two endpoints, in increasing order (as an
      {stpoly_vert_unx_pair_t}).
      
      The fields {ue.minY} and {ue.maxY}, that are not needed during reading,
      are set later, when these pairs of integers are converted to
      records of types {stpoly_vert_rep_t} and {stpoly_edge_rep_t}. */

    /* Lookup table to uniquify edges: */
    size_t sz_edge = sizeof(stpoly_vert_unx_pair_t);
    uint32_t ng_edge = neGuess; /* Expected number of edges. */
    bvtable_t *tb_edge = bvtable_new(sz_edge, ng_edge);
    auto uint64_t hash_edge(void *ep, size_t sz);
    auto int cmp_edge(void *xp, void *yp, size_t sz);

    /* Lookup table to uniquify vertices: */
    size_t sz_vert = sizeof(i2_t);
    uint32_t ng_vert = ng_edge;
    bvtable_t *tb_vert = bvtable_new(sz_vert, ng_vert);
    auto uint64_t hash_vert(void *vp, size_t sz);
    auto int cmp_vert(void *xp, void *yp, size_t sz);
    
    auto void process_STP_edge(int line, stpoly_STP_edge_t *edge);
      /* Procedure that quantizes an STP face {face} and stores it in
         the polygonal figure, if not degenerate. */

    uint32_t ne_read = 0;  /* Number of edges read from the STP file. */
    uint32_t ne_keep = 0;  /* Number of edges retained in the polygonal figure. */
    
    stpoly_STP_read(fileName, binary, &process_STP_edge);
    
    fprintf(stderr, "read %u edges, kept %u\n", ne_read, ne_keep);

    /* Close the tables and get the basic data: */
    uint32_t nv, ne; /* Number of vertices and edges. */
    i2_t *vpos;                     /* Quantized coordinates of vertices. */
    stpoly_vert_unx_pair_t *endv;   /* Edges, as pairs of vertex indices. */
    
    bvtable_close(tb_vert, &nv, (void**)&(vpos)); 
    bvtable_close(tb_edge, &ne, (void**)&(endv)); 
    
    assert(ne == ne_keep);
    fprintf(stderr, "found %u distinct vertices and %u distinct edges\n", nv, ne);

    /* Build the polygonal figure data structure. */
    stpoly_t poly = stpoly_build(eps, nv, vpos, ne, endv);
    
    return poly;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void process_STP_edge(int line, stpoly_STP_edge_t *stp_edge)
      { 
        if (debug) 
          { fprintf(stderr, "edge = "); stpoly_STP_print_edge(stderr, stp_edge); fprintf(stderr, " --> (");  }
        
        ne_read++;
        /* Quantize vertices, assign indices: */
        stpoly_vert_unx_t uxv[2]; /* Indices of endpoint vertices, assigned or recovered. */
        int k;
        for (k = 0; k < 2; k++)
          { /* Quantize the coordinates of endpoint {k}: */
            i2_t vposk = stpoly_STP_round_point(&(stp_edge->v[k]), eps, even);
            if (debug) { fprintf(stderr, " ( %d %d )", vposk.c[0], vposk.c[1]); }
            /* Get the unique vertex index: */
            uxv[k] = bvtable_add(tb_vert, (void*)(&vposk), &hash_vert, &cmp_vert);
            if (debug) { fprintf(stderr, " = v[%d]", uxv[k]); }
            demand(uxv[k] <= stpoly_nv_MAX, "too many vertices in polygonal figure");
          }
        if (debug) { fprintf(stderr, " )"); }
        /* Check for repeated vertices: */
        if (uxv[0] == uxv[1])
          { fprintf(stderr, "%s:%d: !! warning: edge endpoints coincide\n", fileName, line);
            stpoly_STP_print_edge(stderr, stp_edge);
            fprintf(stderr, "\n");
            return;
          }
        
        /* Make sure that {uxv[0]} is the vertex with lowest index: */
        if (uxv[0] > uxv[1]) { stpoly_vert_unx_t t = uxv[0]; uxv[0] = uxv[1]; uxv[1] = t; }

        /* Assign an index to the edge: */
        stpoly_vert_unx_pair_t uxendv;     /* Temporary record for edge {k}. */
        uxendv.v[0] = uxv[0];
        uxendv.v[1] = uxv[1];
        stpoly_edge_unx_t uxe = bvtable_add(tb_edge, (void *)(&uxendv), &hash_edge, &cmp_edge);
        if (debug) { fprintf(stderr, " = e[%d]\n", uxe); }
        if (uxe < ne_keep) 
          { /* Repeated edge: */
            fprintf(stderr, "%s:%d: !! repeated edge, ignored\n", fileName, line);
          }
        else
          { demand(ne_keep < stpoly_ne_MAX, "too many edges in polygonal figure");
            ne_keep++;
          }
      }

    /* Hashing and comparison procedures for faces, edges, and vertices: */
    
    uint64_t hash_edge(void *ep, size_t sz)
      { assert(sz == sizeof(stpoly_vert_unx_pair_t));
        stpoly_vert_unx_pair_t *e = (stpoly_vert_unx_pair_t *)ep;
        /* Requires endpoint indices in increasing order: */
        assert(e->v[0] < e->v[1]);
        /* Hash the endpoint indices: */
        uint64_t h = bvhash_bytes(ep, sizeof(stpoly_vert_unx_pair_t));
        return h;
      }
        
    auto int cmp_edge(void *xp, void *yp, size_t sz)
      { assert(sz == sizeof(stpoly_vert_unx_pair_t));
        stpoly_vert_unx_pair_t *x = (stpoly_vert_unx_pair_t *)xp;
        stpoly_vert_unx_pair_t *y = (stpoly_vert_unx_pair_t *)yp;
        /* Compare the endpoint indices lexicographically: */
        int k;
        for (k = 0; k < 2; k++)
          { stpoly_vert_unx_t uxvx = x->v[k];
            stpoly_vert_unx_t uxvy = y->v[k];
            if (uxvx < uxvy)
              { return -1; }
            else if (uxvx > uxvy)
              { return +1; }
          }
        return 0;
      }
        
    uint64_t hash_vert(void *vp, size_t sz)
      { assert(sz == sizeof(i2_t));
        /* Hash the quantized coords: */
        uint64_t h = bvhash_bytes(vp, sizeof(i2_t));
        return h;
      }
        
    auto int cmp_vert(void *xp, void *yp, size_t sz)
      { assert(sz == sizeof(i2_t));
        i2_t *x = (i2_t *)xp;
        i2_t *y = (i2_t *)yp;
        /* Compare quantized coords lexicographically in order {Y,X}: */
        int k;
        for (k = 0; k < 2; k++)
          { int32_t xk = x->c[1-k];
            int32_t yk = y->c[1-k];
            if (xk < yk)
              { return -1; }
            else if (xk > yk)
              { return +1; }
          }
        return 0;
      }
  }
  
void stpoly_print_vert_degrees(FILE *wr, stpoly_t poly)
  {
    uint32_t nv = stpoly_vert_count(poly);
    uint32_t ne = stpoly_edge_count(poly);
    
    uint32_t deg_max = ne;
    uint32_t d;
     
    /* Scan vertices and count vertices of the same degree: */
    uint32_t *nv_deg = notnull(malloc((deg_max+1)*sizeof(uint32_t)), "no mem"); 
    for (d = 0; d <= deg_max; d++) { nv_deg[d] = 0; }
    stpoly_edge_unx_t uxv;
    for (uxv = 0; uxv < nv; uxv++)
      { stpoly_vert_t v = stpoly_get_vert(poly, uxv);
        d = stpoly_vert_degree(v);
        if (d == 0) 
          { fprintf(stderr, "** vertex %u has zero degree\n", uxv);
            assert(d > 0);
          }
        assert(d <= deg_max);
        nv_deg[d]++;
      }
    
    /* Print degree statistics: */
    bool_t border = FALSE;
    int32_t nv_tot = 0;
    int32_t ne_tot = 0;
    for (d = 0; d <= deg_max; d++)
      { if (nv_deg[d] > 0)
          { fprintf(wr, "there are %u vertices shared by %u edges\n", nv_deg[d], d);
            nv_tot += nv_deg[d];
            ne_tot += nv_deg[d]*d;
            if ((d % 2) != 0) { border = TRUE; }
          }
      }
    assert((ne_tot % 2) == 0);
    ne_tot /= 2;
    assert(nv_tot == nv);
    assert(ne_tot == ne);
    
    /* if there are edges with odd degree, the figure is not a closed surface: */
    if (border) { fprintf(stderr, "!! warning: polygonal figure is not closed\n"); }
    
    free(nv_deg);
  }
     
r2_t stpoly_unround_point(i2_t *p, float eps)
  { 
    r2_t r;
    r.c[0] = ((double)eps)*((double)p->c[0]);
    r.c[1] = ((double)eps)*((double)p->c[1]);
    return r;
  }

void stpoly_print_bounding_box(FILE *wr, stpoly_t poly)
  { 
    float eps = stpoly_get_eps(poly);
    i2_t minQ, maxQ;
    stpoly_get_bounding_box(poly, &minQ, &maxQ);
    fprintf(wr, "bounding box:\n");
    int k;
    for (k = 0; k < 2; k++)
      { int32_t minQk = minQ.c[k]; float minFk = eps*(float)minQk;
        int32_t maxQk = maxQ.c[k]; float maxFk = eps*(float)maxQk;
        fprintf(wr, "  %c [ %d _ %d ] = [ %.8f mm _ %.8f mm ]\n", "XY"[k], minQk, maxQk, minFk, maxFk); 
      }
  }
  
bool_t stpoly_edge_crosses_line(stpoly_edge_t e, int32_t pY)
  {
    bool_t verbose = FALSE;
    stpoly_vert_t v[2]; /* Endpoints of {e}. */
    stpoly_edge_get_endpoints(e, v);
    int32_t vY[2]; /* Quantized {Y}-coords of endpoints of {e}. */
    int r;
    for (r = 0; r < 2; r++) 
      { i2_t p = stpoly_vert_get_pos(v[r]); vY[r] = p.c[1];
        if (verbose) { fprintf(stderr, "endpoint %d at Y = %+d\n", r, vY[r]); }
      }
    assert((vY[0] != pY) && (vY[1] != pY));
    if ((vY[0] < pY) && (vY[1] > pY)) { return TRUE; }
    if ((vY[0] > pY) && (vY[1] < pY)) { return TRUE; }
    return FALSE;
  }

double stpoly_edge_line_intersection(stpoly_edge_t e, int32_t pY, double eps)
  {
    stpoly_vert_t v[2]; /* Endpoints of {e}. */
    stpoly_edge_get_endpoints(e, v);
    i2_t p[2]; /* Quantized coordinates of the endpoints. */
    int r;
    for (r = 0; r < 2; r++) { p[r] = stpoly_vert_get_pos(v[r]); }
    int32_t nY = pY - p[0].c[1];
    int32_t dY = p[1].c[1] - p[0].c[1];
    assert(dY != 0);
    double t = ((double)nY)/((double)dY);
    assert((t > 0.0) && (t < 1.0));
    int32_t dX = p[1].c[0] - p[0].c[0];
    double uX = eps*(p[0].c[0] + t*dX);
    return uX;
  }

stpoly_t stpoly_build
  ( float eps, 
    uint32_t nv, 
    i2_t vpos[], 
    uint32_t ne, 
    stpoly_vert_unx_pair_t endv[]
  )
  {  
    bool_t verify = TRUE;
    
    demand(nv < stpoly_nv_MAX, "too many vertices");
    demand(ne < stpoly_ne_MAX, "too many edges");
   
    stpoly_t poly = stpoly_new_desc(eps, nv, ne);

    /* Build the vertex table: */
    stpoly_vert_unx_t uxv;
    for (uxv = 0; uxv < nv; uxv++)
      { i2_t *vp = &(vpos[uxv]);
        stpoly_vert_unx_t uxv2 = stpoly_add_vert(poly, vp);
        assert(uxv2 == uxv);
      }

    /* Build and fill the edge table: */
    stpoly_edge_unx_t uxe;
    for (uxe = 0; uxe < ne; uxe++)
      { /* Get the indices of the endpoints from the {endv} list: */
        stpoly_vert_unx_pair_t *ev = &(endv[uxe]);
        stpoly_edge_unx_t uxe2 = stpoly_add_edge(poly, &(ev->v[0]));
        assert(uxe2 == uxe);
      }

    stpoly_print_bounding_box(stderr, poly);
    stpoly_print_vert_degrees(stderr, poly);
    
    /* Paranoia: */
    if (verify) { stpoly_check(poly); }
    
    return poly;
  }
    
void stpoly_check(stpoly_t poly)
  {
    /* Check vertices */
    uint32_t nv = stpoly_vert_count(poly);
    stpoly_vert_unx_t uxv;
    for (uxv = 0; uxv < nv; uxv++)
      { stpoly_vert_t v = stpoly_get_vert(poly, uxv);
        stpoly_vert_check(poly, v);
      }

    /* Check edges */
    uint32_t ne = stpoly_edge_count(poly);
    stpoly_edge_unx_t uxe;
    for (uxe = 0; uxe < ne; uxe++)
      { stpoly_edge_t e = stpoly_get_edge(poly, uxe);
        stpoly_edge_check(poly, e);
      }
  }

void stpoly_vert_check(stpoly_t poly, stpoly_vert_t v)
  {
    /* Check the index: */
    stpoly_vert_unx_t uxv = stpoly_vert_get_unx(poly, v);
    assert(v == stpoly_get_vert(poly, uxv));
  }

void stpoly_edge_check(stpoly_t poly, stpoly_edge_t e)
  {
    int k;

    /* Check the index: */
    stpoly_edge_unx_t uxe = stpoly_edge_get_unx(poly, e);
    stpoly_edge_t e0 = stpoly_get_edge(poly, uxe); /* Natural vesion of {e} */
    
    /* Get the endpoints of {e}: */
    stpoly_vert_t ve[2];
    stpoly_edge_get_endpoints(e, ve);
    
    /* Compare with reversed versions: */
    bool_t found = FALSE; /* Found a version of {e} that matches {e0}. */
    stpoly_edge_t c = e;
    int t;
    for (t = 0; t < 2; t++)
      { /* At this point, {c = reverse^t(e)}. */
        stpoly_edge_t c1 = stpoly_edge_reverse(e, t);
        assert(c1 == c);
        
        /* Compare with natural version: */
        if (c == e0) { found = TRUE; }

        /* Check the endpoints of {c}: */
        stpoly_vert_t vc[2];
        stpoly_edge_get_endpoints(c, vc);
        for (k = 0; k < 2; k++)
          { 
            /* Compare endpoint  {vc[k]} with the endpoint of {e}: */
            if (t == 0)
              { assert(vc[k] == ve[k]); }
            else
              { assert(vc[k] == ve[1 - k]); }
            
            /* Check single endpoint: */
            stpoly_vert_t uk = stpoly_edge_get_endpoint(c, k);
            assert(vc[k] == uk);
          }
        c = stpoly_edge_reverse(c, 1);
      }
    assert(c == e);
    assert(found);
  }
   
