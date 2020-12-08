/* Functions of {stmesh.h} that depend on the representation. See {stmesh_rep.h} */
/* Last edited on 2016-04-21 18:44:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <i3.h>
#include <r3.h>
  
#include <stmesh.h>
#include <stmesh_rep.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

#define EDGE_NAT(e) ( (stmesh_edge_t)( ((uint64_t)(e)) & (~1LU) ) )
#define EDGE_D(e) ( (int)( ((uint64_t)(e)) & 1LU ) )
#define EDGE_REV(e) ( (stmesh_edge_t)( ((uint64_t)(e)) ^ 1LU ) )
        
#define FACE_NAT(f) ( (stmesh_face_t)( ((uint64_t)(f)) & (~7LU) ) )
#define FACE_D(f) ( (int)( ((uint64_t)(f)) & 1LU ) )
#define FACE_K(f) ( (int)( ( ((uint64_t)(f)) & 6LU ) >> 1 ) )
#define FACE_FLIP(f) ( (stmesh_face_t)( ((uint64_t)(f)) ^ 1LU ) )

uint32_t stmesh_vert_count(stmesh_t mesh) 
  { return mesh->nv; }
  
uint32_t stmesh_edge_count(stmesh_t mesh) 
  { return mesh->ne; }
  
uint32_t stmesh_face_count(stmesh_t mesh) 
  { return mesh->nf; }
  
float stmesh_get_eps(stmesh_t mesh) 
  { return mesh->eps; }

i3_t stmesh_vert_get_pos(stmesh_vert_t v)
  { return v->pos; }

stmesh_edge_t stmesh_edge_reverse(stmesh_edge_t e, int k)
  { if ((k & 1) == 1)
      { return EDGE_REV(e); }
    else
      { return e; }
  }

stmesh_edge_t stmesh_edge_natural(stmesh_edge_t e)
  { 
    return EDGE_NAT(e);
  }

stmesh_face_t stmesh_face_flip(stmesh_face_t f, int k)
  { if ((k & 1) == 0)
      { return f; }
    else
      { return FACE_FLIP(f); }
  }

stmesh_face_t stmesh_face_shift(stmesh_face_t f, int k)
  { 
    int fd = FACE_D(f);
    int fk = FACE_K(f);
    stmesh_face_t f0 = FACE_NAT(f);
    /* Apply {k} 120 degree rotations in sense {d} to {fk}: */
    k = k % 3;
    fk = (fd == 0 ? fk + k : fk - k);
    /* Reduce {fk} to {0..2} modulo 3: */
    fk = fk % 3;
    if (fk < 0) { fk += 3; }
    /* Reassmeble the address: */
    return (stmesh_face_t)( ((uint64_t)f0) | (uint64_t)((fk << 1) | fd) );
  }

stmesh_face_t stmesh_face_natural(stmesh_face_t f)
  { 
    return FACE_NAT(f);
  }

uint32_t stmesh_edge_degree(stmesh_edge_t e)
  { stmesh_edge_t e0 = EDGE_NAT(e);
    return e0->degree;
  }

stmesh_vert_t stmesh_edge_get_endpoint(stmesh_edge_t e, int k)
  { assert((k == 0) || (k == 1));
    int ed = EDGE_D(e);
    stmesh_edge_t e0 = EDGE_NAT(e);
    return e0->endv[ed ^ k];
  }

void stmesh_edge_get_endpoints(stmesh_edge_t e, stmesh_vert_t v[])
  { int ed = EDGE_D(e);
    stmesh_edge_t e0 = EDGE_NAT(e);
    v[0] = e0->endv[0 ^ ed];
    v[1] = e0->endv[1 ^ ed];
  }


void stmesh_get_bounding_box(stmesh_t mesh, i3_t *minQP, i3_t *maxQP)
  { 
    (*minQP) = mesh->minQ;
    (*maxQP) = mesh->maxQ;
  }


void stmesh_face_get_zrange(stmesh_face_t f, int32_t *minZP, int32_t *maxZP)
  { 
    stmesh_face_t f0 = FACE_NAT(f);
    (*minZP) = f0->minZ;
    (*maxZP) = f0->maxZ;
  }
    
stmesh_edge_t stmesh_face_get_base(stmesh_face_t f)
  { int fk = FACE_K(f);
    int fd = FACE_D(f);
    stmesh_face_t f0 = FACE_NAT(f);
    stmesh_edge_t e = f0->side[fk];
    assert(e != NULL);
    if (fd != 0) { e = EDGE_REV(e); }
    return e;
  }

void stmesh_face_get_sides(stmesh_face_t f, stmesh_edge_t e[])
  { int fk = FACE_K(f);
    int fd = FACE_D(f);
    stmesh_face_t f0 = FACE_NAT(f);
    if (fd == 0)
      { e[0] = f0->side[fk];
        e[1] = f0->side[(fk+1)%3];
        e[2] = f0->side[(fk+2)%3];
      }
    else
      { e[0] = EDGE_REV(f0->side[fk]);
        e[1] = EDGE_REV(f0->side[(fk+2)%3]);
        e[2] = EDGE_REV(f0->side[(fk+1)%3]);
      }
  }

void stmesh_face_get_corners(stmesh_face_t f, stmesh_vert_t v[])
  { stmesh_edge_t e[3];
    stmesh_face_get_sides(f, e);
    int k;
    for (k = 0; k < 3; k++)
      { stmesh_edge_t ek2 = e[(k+2)%3];
        int d = EDGE_D(ek2);
        stmesh_edge_t e0 = EDGE_NAT(ek2);
        v[k] = e0->endv[d];
      }
  }

stmesh_t stmesh_new_desc(float eps, uint32_t nv_max, uint32_t ne_max, uint32_t nf_max)
  { stmesh_t mesh = notnull(malloc(sizeof(stmesh_rep_t)), "no mem");
    mesh->eps = eps;
    
    mesh->nv_max = nv_max;
    mesh->v = notnull(malloc(nv_max*sizeof(stmesh_vert_rep_t)), "no mem");
    mesh->nv = 0;

    mesh->ne_max = ne_max;
    mesh->e = notnull(malloc(ne_max*sizeof(stmesh_edge_rep_t)), "no mem");
    mesh->ne = 0; 

    mesh->nf_max = nf_max;
    mesh->f = notnull(malloc(nf_max*sizeof(stmesh_face_rep_t)), "no mem");
    mesh->nf = 0; 

    mesh->minQ = (i3_t){{ INT32_MAX, INT32_MAX, INT32_MAX }}; /* For now. */
    mesh->maxQ = (i3_t){{ INT32_MIN, INT32_MIN, INT32_MIN }}; /* For now. */
    
    return mesh;
  }

void stmesh_free(stmesh_t mesh)
  {
    free(mesh->v);
    free(mesh->e);
    free(mesh->f);
    free(mesh);
  }

stmesh_vert_t stmesh_get_vert(stmesh_t mesh, stmesh_vert_unx_t uxv)
  {
    demand(uxv < mesh->nv, "invalid vertex index"); 
    return &(mesh->v[uxv]);
  }

stmesh_edge_t stmesh_get_edge(stmesh_t mesh, stmesh_edge_unx_t uxe)
  {
    demand(uxe < mesh->ne, "invalid edge index");
    return &(mesh->e[uxe]);
  }
  
stmesh_face_t stmesh_get_face(stmesh_t mesh, stmesh_face_unx_t uxf)
  {
    demand(uxf < mesh->nf, "invalid face index");
    return &(mesh->f[uxf]);
  }

/* ??? Store the index of the ver, edge,face in its record? ??? */

stmesh_vert_unx_t stmesh_vert_get_unx(stmesh_t mesh, stmesh_vert_t v)
  { 
    assert((v >= mesh->v) && (v < &(mesh->v[mesh->nv])));
    uint64_t uxv = v - (mesh->v);
    assert(uxv < (uint64_t)(mesh->nv));
    return (stmesh_vert_unx_t)uxv;
  }
  
stmesh_edge_unx_t stmesh_edge_get_unx(stmesh_t mesh, stmesh_edge_t e)
  { 
    assert((e >= mesh->e) && (e < &(mesh->e[mesh->ne])));
    uint64_t uxe = e - (mesh->e);
    assert(uxe < (uint64_t)(mesh->ne));
    return (stmesh_edge_unx_t)uxe;
  }
  
stmesh_face_unx_t stmesh_face_get_unx(stmesh_t mesh, stmesh_face_t f)
  { 
    assert((f >= mesh->f) && (f < &(mesh->f[mesh->nf])));
    uint64_t uxf = f - (mesh->f);
    assert(uxf < (uint64_t)(mesh->nf));
    return (stmesh_face_unx_t)uxf;
  }

stmesh_vert_unx_t stmesh_add_vert(stmesh_t mesh, i3_t *pos)
  {
    bool_t debug = FALSE;
    
    /* Allocate the vertex: */
    stmesh_vert_unx_t uxv = mesh->nv;
    demand(uxv < mesh->nv_max, "too many vertices");
    stmesh_vert_rep_t *v = &(mesh->v[uxv]);
    mesh->nv++;

    /* Set the vertex coordinates: */
    v->pos = (*pos);

    /* Update the bounding box: */
    int k;
    for (k = 0; k < 3; k++)
      { int32_t pik = pos->c[k];
        if (pik < mesh->minQ.c[k]) { mesh->minQ.c[k] = pik; }
        if (pik > mesh->maxQ.c[k]) { mesh->maxQ.c[k] = pik; }
      }

    if (debug) { fprintf(stderr, "added v[%u] = ( %d %d %d)\n", uxv, pos->c[0], pos->c[1], pos->c[2]); }

    return uxv;
  }
    
stmesh_edge_unx_t stmesh_add_edge(stmesh_t mesh, stmesh_vert_unx_t uxv[])
  {
    bool_t debug = FALSE;
    
    demand(uxv[0] < uxv[1], "endpoints must be sorted");
    
    /* Allocate the edge: */
    stmesh_edge_unx_t uxe = mesh->ne;
    demand(uxe < mesh->ne_max, "too many edges");
    stmesh_edge_rep_t *e = &(mesh->e[uxe]);
    mesh->ne++;

    /* Set the edge endpoints: */
    int d;
    for (d = 0; d < 2; d++)
      { assert(uxv[d] < mesh->nv);
        e->endv[d] = &(mesh->v[uxv[d]]);
      }

    /* Clear the degree: */
    e->degree = 0;

    if (debug) { fprintf(stderr, "added e[%u] endpoints ( %u %u ) ", uxe, uxv[0], uxv[1]); }

    return uxe;
  }
    
stmesh_face_unx_t stmesh_add_face(stmesh_t mesh, stmesh_edge_unx_t uxe[])
  { 
    bool_t debug = FALSE;
    
    demand((uxe[0] < uxe[1]) && (uxe[1] < uxe[2]), "sides must be sorted");
    
    /* Allocate the edge: */
    stmesh_face_unx_t uxf = mesh->nf;
    demand(uxf < mesh->nf_max, "too many faces");
    stmesh_face_rep_t *f = &(mesh->f[uxf]);
    mesh->nf++;

    if (debug) { fprintf(stderr, "added f[%u] edges ( %u %u %u ) ", uxf, uxe[0], uxe[1], uxe[2]); }

    /* Get the pointers to the edge records, in their natural orientations, and bump degrees: */
    stmesh_edge_rep_t *e[3];
    int k;
    for (k = 0; k < 3; k++) 
      { assert(uxe[k] < mesh->ne);
        e[k] = &(mesh->e[uxe[k]]);
        e[k]->degree++;
      }

    /* The base edge is {fside.c[0]} in natural orientation: */
    f->side[0] = e[0];
    
    /* Get the base edge endpoints: */
    stmesh_vert_t v0 = e[0]->endv[0];
    stmesh_vert_t v1 = e[0]->endv[1];
    
    /* Find {k1} in {1..2} such that {e[k1]} in sense {d1} continues the base edge: */
    int k1 = -1; int d1 = -1;
    for (k = 1; k < 3; k++)
      { if (v1 == e[k]->endv[0]) { k1 = k; d1 = 0; break; }
        else if (v1 == e[k]->endv[1]) { k1 = k; d1 = 1; break; }
      }
    assert(k1 >= 0);
    f->side[1] = stmesh_edge_reverse(e[k1], d1);
    assert(v1 == stmesh_edge_get_endpoint(f->side[1], 0));

    /* Find the proper orientation of the other edge: */
    int k2 = 3 - k1; /* {1-->2, 2-->1}. */
    int d2 = (v0 == e[k2]->endv[0] ? 1 : 0);
    f->side[2] = stmesh_edge_reverse(e[k2], d2);
    assert(v0 == stmesh_edge_get_endpoint(f->side[2], 1));

    /* Compute {f.minZ,f.maxZ}: */
    if (debug) { fprintf(stderr, " vertices "); }

    f->minZ = INT32_MAX;
    f->maxZ = INT32_MIN;
    for (k = 0; k < 3; k++)
      { stmesh_vert_t vk = stmesh_edge_get_endpoint(f->side[k], 0);
        int32_t zk = vk->pos.c[2];
        if (zk < f->minZ) { f->minZ = zk; }
        if (zk > f->maxZ) { f->maxZ = zk; }

        if (debug) { fprintf(stderr, " ( %d %d %d)", vk->pos.c[0], vk->pos.c[1], vk->pos.c[2]); }

        /* Paranoia check: */
        stmesh_vert_t uk = stmesh_edge_get_endpoint(f->side[(k+2)%3], 1);
        assert(uk == vk);
      }
    if (debug) { fprintf(stderr, "\n"); }
     
    return uxf;
  }

void stmesh_vert_print(FILE *wr, stmesh_t mesh, stmesh_vert_t v)
  { 
    stmesh_vert_unx_t uxv = stmesh_vert_get_unx(mesh, v);
    i3_t vQ = stmesh_vert_get_pos(v);
    fprintf(wr, "vertex %u = %p ( %d %d %d )", uxv, v, vQ.c[0], vQ.c[1], vQ.c[2]);
  }

void stmesh_edge_print(FILE *wr, stmesh_t mesh, stmesh_edge_t e)
  { 
    stmesh_edge_unx_t uxe = stmesh_edge_get_unx(mesh, e);
    
    stmesh_edge_t e0 = EDGE_NAT(e);
    fprintf(wr, "edge %u = %p (%p:%d)", uxe, e, e0, EDGE_D(e));

    fprintf(wr, " degree %d", e0->degree);

    stmesh_vert_t ve[2];
    stmesh_edge_get_endpoints(e, ve);
    stmesh_vert_unx_t uxv0 = stmesh_vert_get_unx(mesh, ve[0]);
    stmesh_vert_unx_t uxv1 = stmesh_vert_get_unx(mesh, ve[1]);
    fprintf(wr, " endpoints %u %u", uxv0, uxv1);
  }

void stmesh_face_print(FILE *wr, stmesh_t mesh, stmesh_face_t f)
  { 
    stmesh_face_unx_t uxf = stmesh_face_get_unx(mesh, f);
    
    stmesh_face_t f0 = FACE_NAT(f);
    fprintf(wr, "face %u = %p (%p:%d:%d)", uxf, f, f0, FACE_K(f), FACE_D(f)); 
    fprintf(wr, " Z = [ %d .. %d ]\n", f0->minZ, f0->maxZ); 
    
    /* Print the sides of {f}: */
    stmesh_edge_t ef[3];
    stmesh_face_get_sides(f, ef);
    int k;
    for (k = 0; k < 3; k++)
      { fprintf(wr, "  side %d = ", k);
        stmesh_edge_print(wr, mesh, ef[k]);
        fprintf(wr, " \n");
      }
        
    /* Print the corners of {f}: */
    stmesh_vert_t vf[3];
    stmesh_face_get_corners(f, vf);
    for (k = 0; k < 3; k++)
      { fprintf(wr, "  corner %d = ", k);
        stmesh_vert_print(wr, mesh, vf[k]);
        fprintf(wr, "\n");
      }
  }

