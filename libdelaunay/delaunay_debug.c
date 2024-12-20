/* See {delaunay_debug.h} */
/* Last edited on 2024-12-05 10:25:01 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <quad.h>
#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <hi2.h>

#include <delaunay.h>
#include <delaunay_debug.h>

#define delaunay_COORD_FMT "6d"
  /* Format (minus the "%") to use when printing a site's homogeneous coord. */

/* The exterior face: */
delaunay_face_t *EXT = NULL; /* Created by the first call to {deldebug_make_edge}. */
    
/* List of edges created so far: */
static int32_t nec = 0; /* Number of edges created so far. */
static quad_arc_vec_t ec; /* Edges created so far are {ex.e[0..nec-1]}. */

quad_arc_t deldebug_make_edge(void)
  { 
    quad_arc_t e = quad_make_edge();
    if (DEBUGGING_DELAUNAY)
      { /* Set the exterior face: */
        if (EXT == NULL) { EXT = notnull(malloc(sizeof(delaunay_face_t)), "NO MEM"); }
        deldebug_set_all_left_faces(e, EXT);
      }
    if (DEBUGGING_QUAD)
      { /* Save the edge for further checking: */
        if (nec == 0) { ec = quad_arc_vec_new(100); }
        quad_arc_vec_expand(&ec, nec);
        ec.e[nec] = e;
        nec++;
        deldebug_check_quad_all_edges();
      }
    return e;
  }

void deldebug_set_all_left_faces(quad_arc_t e, delaunay_face_t *F)
  {
    quad_arc_t a = e;
    do { SET_LEFT(a, F); a = quad_lnext(a); } while(a != e);
  }

void deldebug_destroy_edge(quad_arc_t e)
  {
    if (DEBUGGING_DELAUNAY)
      { /* If only one of the faces of {e} is exterior, we must set the other too: */
        assert(EXT != NULL);
        if ((LEFT(e) == EXT) != (RITE(e) == EXT))
          { if (LEFT(e) != EXT) { deldebug_set_all_left_faces(e, EXT); }
            if (RITE(e) != EXT) { deldebug_set_all_left_faces(quad_sym(e), EXT); }
          }
      }
    if (DEBUGGING_QUAD)
      { deldebug_check_quad_all_edges();
        /* Disconnect the edge but do not destroy it: */
        if (quad_onext(e) != e) { quad_splice(e, quad_onext(e)); }
        if (quad_dnext(e) != e) { quad_splice(e, quad_dnext(e)); }
        assert(quad_onext(e) == e);
        assert(quad_dnext(e) == e);
        deldebug_set_all_left_faces(e, EXT);
        /* Remove {e} from the list of created edges: */
        int32_t m = 0;
        for (uint32_t i = 0;  i < nec; i++) 
          { if (quad_edge(e) == quad_edge(ec.e[i])) 
              { ec.e[i] = quad_arc_NULL; m++; }
          }
        assert(m == 1);
        deldebug_check_quad_all_edges();
      }
    else
      { /* Disconnect and destroy: */
        quad_destroy_edge(e);
      }
     return;
  }

void deldebug_check_left_triangle(quad_arc_t b, int32_t depth)
  {
    if (DEBUGGING_DELAUNAY)
      { demand(LEFT(b) != EXT, "non-interior face");
        deldebug_check_interior_face(b);
        if (depth > 0) 
          { quad_arc_t e = quad_onext(b);
            while(e != quad_sym(b))
              { if (LEFT(e) != EXT) 
                  { deldebug_check_delaunay_edge_property(e);
                    deldebug_check_left_triangle(e, depth-1);
                  }
                e = quad_rprev(e);
              }
          }
      }
  }
  
void deldebug_check_delaunay_edge_property(quad_arc_t e)
  {
    if (DEBUGGING_DELAUNAY)
      { demand((LEFT(e) != EXT) && (RITE(e) != EXT), "not an internal edge");
        delaunay_site_t *A = ORG(e);
        delaunay_site_t *C = DST(e);
        delaunay_site_t *B = DST(quad_onext(e));
        delaunay_site_t *D = DST(quad_oprev(e));
        demand(delaunay_leftof(B, e), "inverted triangle");
        demand(delaunay_rightof(D, e), "inverted triangle");
        demand(hi2_in_circle(&(A->pt),&(B->pt),&(C->pt),&(D->pt)), "non-delaunay edge");
      }
  }

void deldebug_check_interior_face(quad_arc_t e)
  {
    int32_t n = 0;
    quad_arc_t a = e;
    do
      { assert(LEFT(a) != EXT);
        /* Check {LEFT(a)} for convexity: */
        quad_arc_t b = quad_lnext(a);
        demand((b != a) && (b != quad_sym(a)), "face is not a polygon");
        b = quad_lnext(b);
        demand((b != a) && (b != quad_sym(a)), "face is not a polygon");
        while (b != a)
          { demand(delaunay_leftof(ORG(b), a), "interior face is not convex");
            b = quad_lnext(b);
          }
        n++;
        a = quad_lnext(a);
      } 
    while(a != e);
    demand(n >= 3, "too few sides in face");
    
    if (n >= 4)
      { /* Check that all sites are cocircular */
        a = e;
        do
          { /* Get the next 4 sites: */
            quad_arc_t b = quad_lnext(a); 
            quad_arc_t c = quad_lnext(b); 
            quad_arc_t d = quad_lnext(c); 
            
            assert((b != a) && (b != quad_sym(a)));
            assert((c != a) && (c != quad_sym(a)));
            assert((c != b) && (c != quad_sym(b)));
            assert((d != a) && (d != quad_sym(a)));
            assert((d != b) && (d != quad_sym(b)));
            assert((d != c) && (d != quad_sym(c)));
            
            hi2_point_t *pa = &(ORG(a)->pt);
            hi2_point_t *pb = &(ORG(b)->pt);
            hi2_point_t *pc = &(ORG(c)->pt);
            hi2_point_t *pd = &(ORG(d)->pt);
            demand(hi2_in_circle(pa, pb, pc, pd) == 0, "face sites are not cocircular");
            a = quad_lnext(a);
          }
        while(a != e);
      }
  }

/* Debugging: */

void deldebug_print_site (char *msg, delaunay_site_t *a)
  { fprintf(stderr, "%s site[%d] =", msg, a->index);
    i3_gen_print(stderr, &(a->pt.c), ("%+" delaunay_COORD_FMT), "[ ", " ", " ]");
  }

void deldebug_print_edge (char *msg, quad_arc_t e)
  { fprintf(stderr, "%s ", msg);
    quad_write_arc(stderr, e, 1);
    deldebug_print_site("  ORG =", ORG(e)); 
    deldebug_print_site("  DST =", DST(e)); 
  }

void deldebug_check_quad_edge(quad_arc_t e)
  {
    for (uint32_t k = 0;  k < 4; k++)
      { assert(! quad_arc_is_null(quad_onext(e))); 
        e = quad_rot(e);
      }
  }

void deldebug_check_quad_all_edges(void)
  {
    for (uint32_t i = 0;  i < nec; i++)
      { quad_arc_t e = ec.e[i];
        if (! quad_arc_is_null(e)) { deldebug_check_quad_edge(e); }
      }
  }

