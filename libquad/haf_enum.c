/* See {haf_enum.h}. */
/* Last edited on 2024-12-05 10:38:51 by stolfi */
 
#define haf_enum_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

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
    
#define debug FALSE

/* IMPLEMENTATIONS */

void haf_enum_edges
  ( haf_arc_count_t nr,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_arc_vec_t *E_P,
    haf_arc_vec_t *C_P
  )
  { 
    /* Provide internal edge table if {E_P} is {NULL}: */
    haf_arc_vec_t E_local = haf_arc_vec_new(0);
    haf_arc_vec_t *E = (E_P != NULL ? E_P : &E_local);
    haf_arc_vec_t *C = C_P; 
    
    haf_edge_count_t ne = 0; /* The edges seen so far are {E.e[0..ne-1]}. */
    haf_edge_count_t nc = 0; /* Number of connected components found so far. */ 
   
    auto haf_arc_t get_next_cand_edge(void);
      /* Returns a next candidate for the next unseen edge.
        If the queue is not empty, namely {pe < ne}, and {a = E.e[pe].lnext}
        is not seen, returns {a}; else returns {b = E[pe].sym.next} and increments {pe}.
        If the queue is empty but the root list is not exhausted, namely {kr < nr},
        return {root[kr]} and increments {kr}. . */
        
    auto bool_t edge_is_seen(haf_arc_t a);
      /* True iff the arc {a} or its opposite have been saved in {E.e[0..ne-1]}. */
    
    int32_t pe = 0;  /* The arcs {E.e[0..pe-1]} and their opposites have been processed. */
    haf_arc_count_t kr = 0; /* The arcs {root[0..kr-1]} have been processed. */
    while ((kr < nr) || (pe < ne))
      { if (debug) { fprintf(stderr, "  kr = %lu pe = %d\n", kr, pe); } 
        /* At this point {E.e[0..ne-1]} are the even-numbered arcs of
          all edges (opposite arc pairs) we have seen (renumbered), and
          {E.e[ke].aid = 2*(ke + eid0)} for all {ke} in {0..ne-1}. Those edges
          include all edges of the arcs {root[0..kr-1]},
          {E.e[0..pe-1].lnext}, and {E.e[0..pe-1].sym.lnext}. The arcs
          {root[kr..nr-1]}, {E.e[pe..ne-1].lnext} and
          {E.e[pe..ne-1].sym.lnext} may still be on yet-unnumbered
          edges. */
          
        /* Get the next arc {a} that is still possibly un-numbered: */
        haf_arc_t a = get_next_cand_edge();
        demand(a != NULL, "{get_next_cand_edge} returned {NULL}");
        if (! edge_is_seen(a)) 
          { /* New edge; renumber and store in {E}: */
            demand(ne <= haf_edge_count_MAX, "too many edges");
            demand(eid0 <= haf_edge_id_MAX - ne, "edge identifier got too big");
            haf_set_edge_id(a, eid0 + ne);
            haf_arc_vec_expand(E, (vec_index_t)ne);
            if (debug) { fprintf(stderr, "    unseen, saving in {E[%lu]}\n", ne); }
            E->e[ne] = haf_base_arc(a);
            ne++;
          }
        else
          { if (debug) { fprintf(stderr, "    already seen, skipped\n"); }
          }
      }
      
    /* Trim or recycle the edge table: */
    if (E == &E_local) 
      { /* Edge table is local, reclaim: */
        assert(E_P == NULL);
        free(E->e);
      } 
    else 
      { /* Edge table is client's, trim: */
        assert(E == E_P);
        assert(E_local.ne == 0);
        haf_arc_vec_trim(E, (vec_index_t)ne);
      }
      
    /* Trim component table, if any: */
    assert(C == C_P);
    if (C != NULL) { haf_arc_vec_trim(C, (vec_index_t)nc); }
    
    return;
      
    bool_t edge_is_seen(haf_arc_t a)
      { assert(a != NULL); /* Already checked, but just in case... */
        haf_edge_id_t eid = haf_edge_id(a);
        if (eid < eid0) { return FALSE; }
        haf_edge_count_t ke = eid - eid0;
        if (ke >= ne) { return FALSE; }
        return (haf_edge(E->e[ke]) == haf_edge(a));
      }
      
    haf_arc_t get_next_cand_edge(void)
      { if (pe < ne)
          { /* Get the next edge from the queue {E.e[pe..ne-1]}: */
            haf_arc_t c = E->e[pe];
            assert(c != NULL); /* We never put {NULL} there. */
            haf_arc_t a = haf_lnext(c); demand(a != NULL, "{.lnext} pointer is {NULL}");
            if (debug) { fprintf(stderr, "  trying E[%u].lnext = %lu:%u\n", pe, haf_edge_id(a), haf_dir_bit(a)); } 
            if (! edge_is_seen(a)) { return a; /* Note the other {.lnext} may be unseen too. */ }
            if (debug) { fprintf(stderr, "    already seen, skipped\n"); }
            haf_arc_t b = haf_lnext(haf_sym(c)); demand(b != NULL, "{.sym.lnext} pointer is {NULL}");
            if (debug) { fprintf(stderr, "  trying E[%d].sym.lnext = %lu:%u\n", pe, haf_edge_id(b), haf_dir_bit(b)); } 
            pe++; 
            return b;
          }
        else 
          { /* Get the next root edge: */
            assert(kr < nr);
            haf_arc_t r = root[kr];
            demand(r != NULL, "root arc is {NULL}");
            if (debug) { fprintf(stderr, "  trying root[%lu] = %lu:%u\n", kr, haf_edge_id(r), haf_dir_bit(r)); } 
            if (! edge_is_seen(r))
              { if (C_P != NULL)
                  { /* Another connected component: */
                    if (debug) { fprintf(stderr, "    new component, saved in {C[%lu]}\n", nc); }
                    haf_arc_vec_expand(C_P, (vec_index_t)nc);
                    /* Pick the base edge: */
                    C_P->e[nc] = haf_base_arc(r);
                    nc++;
                  }
              }
            kr++;
            return r;
          }
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
        for (uint32_t kd = 0;  kd < 2; kd++)
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
