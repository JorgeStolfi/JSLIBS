/* Functions of {stpoly.h} that depend on the representation. See {stpoly_rep.h} */
/* Last edited on 2022-10-20 05:59:33 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <i2.h>
#include <r2.h>
  
#include <stpoly.h>
#include <stpoly_rep.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

#define EDGE_NAT(e) ( (stpoly_edge_t)( ((uint64_t)(e)) & (~1LU) ) )
#define EDGE_D(e) ( (int32_t)( ((uint64_t)(e)) & 1LU ) )
#define EDGE_REV(e) ( (stpoly_edge_t)( ((uint64_t)(e)) ^ 1LU ) )

uint32_t stpoly_vert_count(stpoly_t mesh) 
  { return mesh->nv; }
  
uint32_t stpoly_edge_count(stpoly_t mesh) 
  { return mesh->ne; }
  
float stpoly_get_eps(stpoly_t mesh) 
  { return mesh->eps; }

i2_t stpoly_vert_get_pos(stpoly_vert_t v)
  { return v->pos; }

stpoly_edge_t stpoly_edge_reverse(stpoly_edge_t e, int32_t k)
  { if ((k & 1) == 1)
      { return EDGE_REV(e); }
    else
      { return e; }
  }

stpoly_edge_t stpoly_edge_natural(stpoly_edge_t e)
  { 
    return EDGE_NAT(e);
  }

uint32_t stpoly_vert_degree(stpoly_vert_t v)
  { return v->degree;
  }

stpoly_vert_t stpoly_edge_get_endpoint(stpoly_edge_t e, int32_t k)
  { assert((k == 0) || (k == 1));
    int32_t ed = EDGE_D(e);
    stpoly_edge_t e0 = EDGE_NAT(e);
    return e0->endv[ed ^ k];
  }

void stpoly_edge_get_endpoints(stpoly_edge_t e, stpoly_vert_t v[])
  { int32_t ed = EDGE_D(e);
    stpoly_edge_t e0 = EDGE_NAT(e);
    v[0] = e0->endv[0 ^ ed];
    v[1] = e0->endv[1 ^ ed];
  }

void stpoly_get_bounding_box(stpoly_t mesh, i2_t *minQP, i2_t *maxQP)
  { 
    (*minQP) = mesh->minQ;
    (*maxQP) = mesh->maxQ;
  }

void stpoly_edge_get_yrange(stpoly_edge_t e, int32_t *minYP, int32_t *maxYP)
  { 
    stpoly_edge_t e0 = EDGE_NAT(e);
    (*minYP) = e0->minY;
    (*maxYP) = e0->maxY;
  }

stpoly_t stpoly_new_desc(float eps, uint32_t nv_max, uint32_t ne_max)
  { stpoly_t mesh = notnull(malloc(sizeof(stpoly_rep_t)), "no mem");
    mesh->eps = eps;
    
    mesh->nv_max = nv_max;
    mesh->v = notnull(malloc(nv_max*sizeof(stpoly_vert_rep_t)), "no mem");
    mesh->nv = 0;

    mesh->ne_max = ne_max;
    mesh->e = notnull(malloc(ne_max*sizeof(stpoly_edge_rep_t)), "no mem");
    mesh->ne = 0; 

    mesh->minQ = (i2_t){{ INT32_MAX, INT32_MAX }}; /* For now. */
    mesh->maxQ = (i2_t){{ INT32_MIN, INT32_MIN }}; /* For now. */
    
    return mesh;
  }
        
void stpoly_free(stpoly_t poly)
  {
    free(poly->v);
    free(poly->e);
    free(poly);
  }
  
stpoly_vert_t stpoly_get_vert(stpoly_t mesh, stpoly_vert_unx_t uxv)
  {
    demand(uxv < mesh->nv, "invalid vertex index"); 
    return &(mesh->v[uxv]);
  }

stpoly_edge_t stpoly_get_edge(stpoly_t mesh, stpoly_edge_unx_t uxe)
  {
    demand(uxe < mesh->ne, "invalid edge index");
    return &(mesh->e[uxe]);
  }
 
/* ??? Store the index of the vert/edge in its record? ??? */

stpoly_vert_unx_t stpoly_vert_get_unx(stpoly_t mesh, stpoly_vert_t v)
  { 
    assert((v >= mesh->v) && (v < &(mesh->v[mesh->nv])));
    uint64_t uxv = v - (mesh->v);
    assert(uxv < (uint64_t)(mesh->nv));
    return (stpoly_vert_unx_t)uxv;
  }
  
stpoly_edge_unx_t stpoly_edge_get_unx(stpoly_t mesh, stpoly_edge_t e)
  { 
    assert((e >= mesh->e) && (e < &(mesh->e[mesh->ne])));
    uint64_t uxe = e - (mesh->e);
    assert(uxe < (uint64_t)(mesh->ne));
    return (stpoly_edge_unx_t)uxe;
  }

stpoly_vert_unx_t stpoly_add_vert(stpoly_t mesh, i2_t *pos)
  {
    bool_t debug = FALSE;
    
    /* Allocate the vertex: */
    stpoly_vert_unx_t uxv = mesh->nv;
    demand(uxv < mesh->nv_max, "too many vertices");
    stpoly_vert_rep_t *v = &(mesh->v[uxv]);
    mesh->nv++;

    /* Set the vertex coordinates: */
    v->pos = (*pos);

    /* Update the bounding box: */
    int32_t k;
    for (k = 0; k < 3; k++)
      { int32_t pik = pos->c[k];
        if (pik < mesh->minQ.c[k]) { mesh->minQ.c[k] = pik; }
        if (pik > mesh->maxQ.c[k]) { mesh->maxQ.c[k] = pik; }
      }

    /* Clear the degree: */
    v->degree = 0;

    if (debug) { fprintf(stderr, "added v[%u] = ( %d %d %d)\n", uxv, pos->c[0], pos->c[1], pos->c[2]); }

    return uxv;
  }
    
stpoly_edge_unx_t stpoly_add_edge(stpoly_t mesh, stpoly_vert_unx_t uxv[])
  {
    bool_t debug = FALSE;
    
    demand(uxv[0] < uxv[1], "endpoints must be sorted");
    
    /* Allocate the edge: */
    stpoly_edge_unx_t uxe = mesh->ne;
    demand(uxe < mesh->ne_max, "too many edges");
    stpoly_edge_rep_t *e = &(mesh->e[uxe]);
    mesh->ne++;

    /* Set the edge endpoints: */
    int32_t d;
    for (d = 0; d < 2; d++)
      { assert(uxv[d] < mesh->nv);
        e->endv[d] = &(mesh->v[uxv[d]]);
      }

    /* Get the pointers to the vertex records, and bump degrees, compute {e.minY,e.maxY}: */
    if (debug) { fprintf(stderr, " vertices "); }
    e->minY = INT32_MAX;
    e->maxY = INT32_MIN;
    int32_t k;
    for (k = 0; k < 2; k++) 
      { assert(uxv[k] < mesh->nv);
        stpoly_vert_t vk = &(mesh->v[uxv[k]]);
        if (debug) { fprintf(stderr, " ( %d %d %d)", vk->pos.c[0], vk->pos.c[1], vk->pos.c[2]); }
        int32_t yk = vk->pos.c[1];
        if (yk < e->minY) { e->minY = yk; }
        if (yk > e->maxY) { e->maxY = yk; }

        /* Paranoia check: */
        stpoly_vert_t uk = stpoly_edge_get_endpoint(e, k);
        assert(uk == vk);
        vk->degree++;
      }
    if (debug) { fprintf(stderr, "\n"); }

    if (debug) { fprintf(stderr, "added e[%u] endpoints ( %u %u ) ", uxe, uxv[0], uxv[1]); }

    return uxe;
  }
    
void stpoly_vert_print(FILE *wr, stpoly_t mesh, stpoly_vert_t v)
  { 
    stpoly_vert_unx_t uxv = stpoly_vert_get_unx(mesh, v);
    i2_t vQ = stpoly_vert_get_pos(v);
    fprintf(wr, "vertex %u = %p ( %d %d %d )", uxv, v, vQ.c[0], vQ.c[1], vQ.c[2]);

    fprintf(wr, " degree %d", v->degree);
  }

void stpoly_edge_print(FILE *wr, stpoly_t mesh, stpoly_edge_t e)
  { 
    stpoly_edge_unx_t uxe = stpoly_edge_get_unx(mesh, e);
    
    stpoly_edge_t e0 = EDGE_NAT(e);
    fprintf(wr, "edge %u = %p (%p:%d)", uxe, e, e0, EDGE_D(e));
    fprintf(wr, " Y = [ %d .. %d ]\n", e0->minY, e0->maxY); 

    stpoly_vert_t ve[2];
    stpoly_edge_get_endpoints(e, ve);
    stpoly_vert_unx_t uxv0 = stpoly_vert_get_unx(mesh, ve[0]);
    stpoly_vert_unx_t uxv1 = stpoly_vert_get_unx(mesh, ve[1]);
    fprintf(wr, " endpoints %u %u", uxv0, uxv1);
  }

