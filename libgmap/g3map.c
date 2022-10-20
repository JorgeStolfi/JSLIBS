/* See {g3map.h}  */
/* Last edited on 2022-10-20 06:29:17 by stolfi */

#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>
#include <affirm.h>

#include <gem.h>

#include <g3map.h>

// #define gmap_INI_NODES 500  
  /* Initial size of node stacks. */

g3map_place_t g3map_step(g3map_place_t a, int32_t i)
  { 
    /* Checks whether {a} lies on an element of dimension {i}: */ 
    int32_t k;
    for (k = 0; k < i; k++) { assert(a != gem_step(a, k)); }
    /* Take the step: */
    return gem_step(a, i);
  }

g3map_place_t g3map_vert_make(void)
  { return gem_node_new(3); }

g3map_place_t g3map_edge_make(g3map_place_t u, g3map_place_t v)
  { 
    /* Make sure that the vertices are unattached:  */
    assert(u != v);
    int32_t k;
    for (k = 0; k < 3; k++)
      { assert(u == gem_step(u, k));
        assert(v == gem_step(v, k));
      }
    gem_splice(u, v, 0);
    return u;
  }
  
void g3map_edge_splice(g3map_place_t a, g3map_place_t b)
  {
    /* Make sure that {a,b} are not the same end of the same edge: */
    demand(a != b, "cannot splice an edge end to itself");

    /* Make sure that the edges are not attached to a cell:  */
    int32_t k;
    for (k = 2; k < 3; k++)
      { assert(a == gem_step(a, k));
        assert(b == gem_step(b, k));
      }
    gem_splice(a, b, 1);
  }

g3map_place_t g3map_face_make(int32_t n, g3map_place_t e[])
  {
    int32_t i;
    for (i = 0; i < n; i++)
      { /* Get the two edge ends to attach: */
        g3map_place_t b0 = gem_step(e[i],0);
        g3map_place_t a1 = e[(i+1)%n];
        /* Make sure that the edges have loose ends: */
        assert(b0 == gem_step(b0, 1));
        assert(a1 == gem_step(a1, 1));
        /* Attach the two ends: */
        g3map_edge_splice(b0, a1);
      }
    return e[0];
  }

void g3map_face_splice(g3map_place_t f, g3map_place_t g)
  {
    /* Make sure that the edges are not attached to a cell:  */
    assert(f == gem_step(f, 3));
    assert(g == gem_step(g, 3));
    
    g3map_place_t fs = gem_step(f, 0); /* Other end of edge {f} */
    g3map_place_t gs = gem_step(g, 0); /* Other end of edge {g} */
    
    /* Make sure that f and g are not the same or opposite ends of the same edge: */
    assert(f != gs);
    assert(f != g);
    
    gem_splice(f, g, 2);
    gem_splice(fs, gs, 2);
  }

#define g3map_INI_NODES 100
  /* Initial size of facial node list. */

void g3map_cell_splice(gem_ref_t a, gem_ref_t b)
  {
    gem_ref_vec_t nodesA = gem_ref_vec_new(g3map_INI_NODES);
    gem_ref_vec_t nodesB = gem_ref_vec_new(g3map_INI_NODES);
    int32_t nnA = 0, nnB = 0;

    gem_traverse(a, 1, &nodesA, &nnA);
    gem_traverse(b, 1, &nodesB, &nnB);
    demand(nnA == nnB, "faces do not match");

    int32_t i;
    for (i = 0; i < nnA; i++) { gem_splice(nodesA.e[i], nodesB.e[i], 3); }
  }
