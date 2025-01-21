/* See {haf_enum.h}. */
/* Last edited on 2025-01-10 00:16:10 by stolfi */
 
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

void haf_enum_cycles(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, uint64_t cid[], bool_t face);
  /* If {face} is true, does {haf_get_faces(NE,a,cid)}, else does 
    {haf_get_verts(NE,a,cid)}. */
    
#define debug FALSE

/* IMPLEMENTATIONS */

void haf_enum_edges
  ( haf_arc_count_t NR,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_edge_vec_t *E_P,
    haf_edge_id_vec_t *oid_P,
    haf_arc_vec_t *C_P
  )
  { 
    /* Provide internal edge tables if args are {NULL}: */
    haf_edge_vec_t E_local = haf_edge_vec_new(0);
    haf_edge_vec_t *E = (E_P != NULL ? E_P : &E_local);
    
    haf_arc_vec_t C_local = haf_arc_vec_new(0);
    haf_arc_vec_t *C = (C_P != NULL ? C_P : &C_local);
    
    haf_edge_id_vec_t oid_local = haf_edge_id_vec_new(0);
    haf_edge_id_vec_t *oid = (oid_P != NULL ? oid_P : &oid_local);
    
    haf_edge_count_t NE = 0; /* The edges seen so far are {E.e[0..NE-1]}. */
    haf_edge_count_t NC = 0; /* Number of connected components found so far. */ 
   
    auto haf_arc_t get_next_cand_edge(void);
      /* Returns a next candidate for the next unseen edge.
        If the queue is not empty, namely {pe < NE}, and {a = E.e[pe].lnext}
        is not seen, returns {a}; else returns {b = E[pe].sym.next} and increments {pe}.
        If the queue is empty but the root list is not exhausted, namely {kr < NR},
        return {root[kr]} and increments {kr}. . */
        
    auto bool_t edge_is_seen(haf_arc_t a);
      /* True iff the arc {a} or its opposite have been saved in {E.e[0..NE-1]}. */
    
    uint32_t pe = 0;  /* The arcs {E.e[0..pe-1]} and their opposites have been processed. */
    haf_arc_count_t kr = 0; /* The arcs {root[0..kr-1]} have been processed. */
    while ((kr < NR) || (pe < NE))
      { if (debug) { fprintf(stderr, "  kr = %lu pe = %d\n", kr, pe); } 
        /* At this point {E.e[0..NE-1]} are the even-numbered arcs of
          all edges (opposite arc pairs) we have seen (renumbered), and
          {E.e[ke].aid = 2*(ke + eid0)} for all {ke} in {0..NE-1}. Those edges
          include all edges of the arcs {root[0..kr-1]},
          {E.e[0..pe-1].lnext}, and {E.e[0..pe-1].sym.lnext}. The arcs
          {root[kr..NR-1]}, {E.e[pe..NE-1].lnext} and
          {E.e[pe..NE-1].sym.lnext} may still be on yet-unnumbered
          edges. */
          
        /* Get the next arc {a} that is still possibly un-numbered: */
        haf_arc_t a = get_next_cand_edge();
        demand(a != NULL, "{get_next_cand_edge} returned {NULL}");
        if (! edge_is_seen(a)) 
          { /* New edge; renumber and store in {E}: */
            demand(NE <= haf_edge_count_MAX, "too many edges");
            demand(eid0 <= haf_edge_id_MAX - NE, "edge identifier got too big");
            haf_edge_t ed = haf_edge(a);
            /* Save old edge id: */
            haf_edge_id_vec_expand(oid, (vec_index_t)NE);
            oid->e[NE] = haf_edge_id(ed);
            /* Renumber edge: */
            haf_set_edge_id(ed, eid0 + NE);
            /* Append to{E}: */
            haf_edge_vec_expand(E, (vec_index_t)NE);
            if (debug) { fprintf(stderr, "    unseen, saving in {E[%lu]}\n", NE); }
            E->e[NE] = ed;
            NE++;
          }
        else
          { if (debug) { fprintf(stderr, "    already seen, skipped\n"); } }
      }
      
    /* Trim or recycle the edge table: */
    if (E == &E_local) 
      { assert(E_P == NULL); free(E->e); } 
    else 
      { assert(E == E_P);  assert(E_local.e == NULL);
        haf_edge_vec_trim(E, (vec_size_t)NE);
      }
       
    /* Trim or recycle the old edge id table: */
    if (oid == &oid_local) 
      { assert(oid_P == NULL); free(oid->e); } 
    else 
      { assert(oid == oid_P); assert(oid_local.e == NULL);
        haf_edge_id_vec_trim(oid, (vec_size_t)NE);
      }
     
    /* Trim or recycle the component table: */
    if (C == &C_local) 
      { assert(C_P == NULL); free(C->e); } 
    else 
      { assert(C == C_P); assert(C_local.e == NULL);
        haf_arc_vec_trim(C, (vec_size_t)NC);
      }
    
    return;
      
    bool_t edge_is_seen(haf_arc_t a)
      { assert(a != NULL); /* Already checked, but just in case... */
        haf_edge_id_t eid = haf_edge_id(haf_edge(a));
        if (eid < eid0) { return FALSE; }
        haf_edge_count_t ke = eid - eid0;
        if (ke >= NE) { return FALSE; }
        return E->e[ke] == haf_edge(a);
      }
      
    haf_arc_t get_next_cand_edge(void)
      { if (pe < NE)
          { /* Get the next edge from the queue {E.e[pe..NE-1]}: */
            haf_edge_t edp = E->e[pe];
            assert(edp != NULL); /* We never put {NULL} there. */
            haf_arc_t a = haf_lnext(haf_base_arc(edp)); demand(a != NULL, "{.lnext} pointer is {NULL}");
            if (debug) { fprintf(stderr, "  trying E[%u].lnext = %lu:%u\n", pe, haf_edge_id(haf_edge(a)), haf_dir_bit(a)); } 
            if (! edge_is_seen(a)) { return a; /* Note the other {.lnext} may be unseen too. */ }
            if (debug) { fprintf(stderr, "    already seen, skipped\n"); }
            haf_arc_t b = haf_lnext(haf_sym(a)); demand(b != NULL, "{.sym.lnext} pointer is {NULL}");
            if (debug) { fprintf(stderr, "  trying E[%d].sym.lnext = %lu:%u\n", pe, haf_edge_id(haf_edge(b)), haf_dir_bit(b)); } 
            pe++; 
            return b;
          }
        else 
          { /* Get the next root edge: */
            assert(kr < NR);
            haf_arc_t r = root[kr];
            demand(r != NULL, "root arc is {NULL}");
            if (debug) { fprintf(stderr, "  trying root[%lu] = %lu:%u\n", kr, haf_edge_id(haf_edge(r)), haf_dir_bit(r)); } 
            if (! edge_is_seen(r))
              { if (C_P != NULL)
                  { /* Another connected component: */
                    if (debug) { fprintf(stderr, "    new component, saved in {C[%lu]}\n", NC); }
                    haf_arc_vec_expand(C_P, (vec_index_t)NC);
                    /* Pick the base edge: */
                    C_P->e[NC] = haf_base_arc(haf_edge(r));
                    NC++;
                  }
              }
            kr++;
            return r;
          }
      }
  }
  
void haf_enum_faces(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, haf_face_id_t fid[])
  { haf_enum_cycles(NE, a, eid0, fid, FALSE); }
  
void haf_enum_verts(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, haf_vert_id_t vid[])
  { haf_enum_cycles(NE, a, eid0, vid, TRUE); }
  
void haf_enum_cycles(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, uint64_t cid[], bool_t face)
  { 
    demand(NE <= haf_edge_count_MAX, "too many edges");
    demand(eid0 <= haf_edge_id_MAX - NE + 1, "edge identifiers too big");
    
    haf_arc_count_t NA = 2*NE;   /* Number of arcs in structure and size of {cid}. */
    haf_arc_id_t aid0 = 2*eid0;  /* Lowest arc ID. */
    for (haf_arc_count_t ka = 0; ka < NA; ka++) { cid[ka] = UINT64_MAX; }
    
    uint32_t NC = 0;  /* Count of cycles identified and labeled so far. */
    uint32_t pe = 0;  /* Count of edges whose cycles have been identified and labeled. */
    
    while (pe < NE)
      { /* At this point the cycles of both arcs on the edges
          {a[0..pe-1]} have been processed. Identifiers {0..NC-1} have
          been assigned to those cycles, and stored in {cid[u.aid]} for
          every arc {u} in those cycles. The arcs {a[pe..NE-1]} and/or their
          opposites may belong to cycles that have not been labeled and
          traced yet. If the cycle of arc {u} has not been processed
          yet, then {cid[u.aid - aid0]} is {UINT64_MAX}. */
          
        /* Get the next candidate arc {a} whose cycle may be unseen: */
        haf_arc_t b = a[2*pe];
        demand(b != NULL, "given arc is {NULL}");
        /* Check {b} and its opposite: */
        for (uint32_t kd = 0;  kd < 2; kd++)
          { haf_edge_id_t beid = haf_edge_id(haf_edge(b));
            demand(beid < haf_edge_id_MAX, "invalid edge id");
            demand((beid >= eid0) && (beid < eid0 + NE), "edge id out of range");
            haf_arc_count_t bi = haf_arc_id(b) - aid0; /* Index of {b} in {cid}. */
            assert(bi < NA);
            if (cid[bi] == UINT64_MAX)
              { /* The cycle of {b} is new; assign number {NC} and set {cid[]}: */
                haf_arc_t u = b;
                haf_edge_id_t ueid = beid;
                haf_arc_count_t ui = bi;
                while (TRUE) 
                  { assert(ui < NA);
                    cid[ui] = NC;
                    u = (face ? haf_lnext(u) : haf_lnext(haf_sym(u)));
                    if (u == b) { break; }
                    ueid = haf_edge_id(haf_edge(u));
                    demand(ueid < haf_edge_id_MAX, "invalid edge id");
                    demand((ueid >= eid0) && (ueid < eid0 + NE), "edge id out of range");
                    uint64_t ui = haf_arc_id(u) - aid0;
                    demand(cid[ui] == UINT64_MAX, "malformed {.lnext} links");
                  }
                NC++;
              }
            b = haf_sym(b);
          }
        assert(b == a[2*pe]); /* {.sym} is always an involution. */
      }
  }
