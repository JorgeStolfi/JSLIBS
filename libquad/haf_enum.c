/* See {haf_enum.h}. */
/* Last edited on 2023-10-05 19:37:33 by stolfi */
 
#define haf_enum_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>

#include <haf.h>
#include <haf_enum.h>

/* INTERNAL PROTOTYPES */

void haf_enum_cycles(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, uint64_t cid[], bool_t face);
  /* If {face} is true, does {haf_get_faces(ne,a,cid)}, else does 
    {haf_get_verts(ne,a,cid)}. */

/* IMPLEMENTATIONS */

void haf_enum_edges
  ( haf_arc_count_t nr,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_arc_vec_t *E,
    haf_arc_vec_t *C
  )
  { 
    /* Provide an internal edge table if {E} is {NULL} */
    haf_arc_vec_t E_local = haf_arc_vec_new(0);
    if (E == NULL) { E = &E_local; }
      
    haf_edge_count_t ne = 0; /* The edges seen so far are {E.e[0..ne-1]}. */
    haf_edge_count_t nc = 0; /* Number of connected components found so far. */ 
   
    auto bool_t edge_is_seen(haf_arc_t a);
      /* True iff the arc {a} or its opposite have been saved in {E.e[0..ne-1]}. */
    
    int32_t pe = 0;  /* The arcs {E.e[0..pe-1]} and their opposites have been processed. */
    haf_arc_count_t kr = 0; /* The arcs {root[0..kr-1]} have been processed. */
    while ((kr < nr) || (pe < ne))
      { /* At this point {E.e[0..ne-1]} are the even-numbered arcs of
          all edges (opposite arc pairs) we have seen (renumbered), and
          {E.e[ke].aid = 2*(ke + eid0)} for all {ke} in {0..ne-1}. Those edges
          include all edges of the arcs {root[0..kr-1]},
          {E.e[0..pe-1].lnext}, and {E.e[0..pe-1].sym.lnext}. The arcs
          {root[kr..nr-1]}, {E.e[pe..ne-1].lnext} and
          {E.e[pe..ne-1].sym.lnext} may still be on yet-unnumbered
          edges. */
          
        /* Get the next arc {a} that is still possibly un-numbered: */
        haf_arc_t a;
        if (pe < ne)
          { /* Get the next edge from the queue {E.e[pe..ne-1]}: */
            haf_arc_t c = E->e[pe];
            assert(c != NULL); /* We never put {NULL} there. */
            haf_arc_t a = haf_lnext(c); demand(a != NULL, "{.lnext} pointer is {NULL}");
            if (! edge_is_seen(a)) { break; /* Note the other {.lnext} may be unseen too. */ }
            a = haf_lnext(haf_sym(c)); demand(a != NULL, "{.sym.lnext} pointer is {NULL}");
            if (! edge_is_seen(a)) { pe++; break; }
            a = NULL;
          }
        else 
          { /* Get the next root edge: */
            assert(kr < nr);
            a = root[kr];
            demand(a != NULL, "root arc is {NULL}");
            if (! edge_is_seen(a))
              { if (C != NULL)
                  { /* Another connected component: */
                    haf_arc_vec_expand(C, (vec_index_t)nc);
                    C->e[nc] = a;
                    nc++;
                  }
              }
            else
              { a = NULL; }
            kr++;
          }
        /* Now {a} is either {NULL} or an arc of a still-unseen edge: */
        if (a != NULL)
          { /* New edge; renumber and store in {E}: */
            demand(ne <= haf_edge_count_MAX, "too many edges");
            demand(eid0 <= haf_edge_id_MAX - ne, "edge identifier got too big");
            haf_set_edge_id(a, eid0 + ne);
            haf_arc_vec_expand(E, (vec_index_t)ne);
            E->e[ne] = a;
            ne++;
          }
      }
    if (E == &E_local)
      { /* Free storage: */ haf_arc_vec_trim(E, 0); }
    else
      { /* Trim table for caller: */ haf_arc_vec_trim(E, (vec_index_t)ne); }
    if (C != NULL) { haf_arc_vec_trim(C, (vec_index_t)nc); }
    return;
      
    bool_t edge_is_seen(haf_arc_t a)
      { assert(a != NULL); /* Already checked, but just in case... */
        haf_edge_id_t eid = haf_edge_id(a);
        if (eid < eid0) { return FALSE; }
        haf_edge_count_t ke = eid - eid0;
        if (ke > ne) { return FALSE; }
        return (haf_edge(E->e[ke]) == haf_edge(a));
      }
  }
  
void haf_enum_faces(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, haf_face_id_t fid[])
  { haf_enum_cycles(ne, a, eid0, fid, FALSE); }
  
void haf_enum_verts(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, haf_vert_id_t vid[])
  { haf_enum_cycles(ne, a, eid0, vid, TRUE); }
  
void haf_enum_cycles(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, uint64_t cid[], bool_t face)
  { 
    demand(ne <= haf_edge_count_MAX, "too many edges");
    demand(eid0 <= haf_edge_id_MAX - ne + 1, "edge identifiers too big");
    
    haf_arc_count_t na = 2*ne;   /* Number of arcs in structure and size of {cid}. */
    haf_arc_id_t aid0 = 2*eid0;  /* Lowest arc ID. */
    for (haf_arc_count_t ka = 0; ka < na; ka++) { cid[ka] = UINT64_MAX; }
    
    int32_t nc = 0;  /* Count of cycles identified and labeled so far. */
    int32_t pe = 0;  /* Count of edges whose cycles have been identified and labeled. */
    
    while (pe < ne)
      { /* At this point the cycles of both arcs on the edges
          {a[0..pe-1]} have been processed. Identifiers {0..nc-1} have
          been assigned to those cycles, and stored in {cid[u.aid]} for
          every arc {u} in those cycles. The arcs {a[pe..ne-1]} and/or their
          opposites may belong to cycles that have not been labeled and
          traced yet. If the cycle of arc {u} has not been processed
          yet, then {cid[u.aid - aid0]} is {UINT64_MAX}. */
          
        /* Get the next candidate arc {a} whose cycle may be unseen: */
        haf_arc_t b = a[2*pe];
        demand(b != NULL, "given arc is {NULL}");
        /* Check {b} and its opposite: */
        for (int32_t kd = 0; kd < 2; kd++)
          { haf_edge_id_t beid = haf_edge_id(b);
            demand(beid < haf_edge_id_MAX, "invalid edge id");
            demand((beid >= eid0) && (beid < eid0 + ne), "edge id out of range");
            haf_arc_count_t bi = haf_arc_id(b) - aid0; /* Index of {b} in {cid}. */
            assert(bi < na);
            if (cid[bi] == UINT64_MAX)
              { /* The cycle of {b} is new; assign number {nc} and set {cid[]}: */
                haf_arc_t u = b;
                haf_edge_id_t ueid = beid;
                haf_arc_count_t ui = bi;
                while (TRUE) 
                  { assert(ui < na);
                    cid[ui] = nc;
                    u = (face ? haf_lnext(u) : haf_lnext(haf_sym(u)));
                    if (u == b) { break; }
                    ueid = haf_edge_id(u);
                    demand(ueid < haf_edge_id_MAX, "invalid edge id");
                    demand((ueid >= eid0) && (ueid < eid0 + ne), "edge id out of range");
                    uint64_t ui = haf_arc_id(u) - aid0;
                    demand(cid[ui] == UINT64_MAX, "malformed {.lnext} links");
                  }
                nc++;
              }
            b = haf_sym(b);
          }
        assert(b == a[2*pe]); /* {.sym} is always an involution. */
      }
  }
