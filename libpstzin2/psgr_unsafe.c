/* See {psgr_unsafe.h} */
/* Last edited on 2025-04-27 02:58:06 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <filefmt.h>
#include <jsstring.h>
#include <affirm.h>

#include <psgr_types.h>
#include <psgr.h>

#include <psgr_unsafe.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED
  
void psgr_unsafe_arc_set_onext(psgr_t *gr, psgr_arc_t a, psgr_arc_t b)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    ed->enext[db] = b;
  }
   
void psgr_unsafe_arc_set_oprev(psgr_t *gr, psgr_arc_t a, psgr_arc_t b)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    ed->eprev[db] = b;
  }

void psgr_unsafe_arc_set_org(psgr_t *gr, psgr_arc_t a, psgr_vertex_t org)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    ed->org[db] = org;
  }
  
void psgr_unsafe_vertex_set_out_arc(psgr_t *gr, psgr_vertex_t v, psgr_arc_t a)
  { if (a.ix != ANONE.ix) { demand(psgr_arc_org(gr, a).ix == v.ix, "origin or arc {a} is not {v}"); }
    gr->vdata[v.ix].aout = a;
  }
  
void psgr_arc_unsafe_split
  ( psgr_t *gr,
    psgr_arc_t a,
    psgr_edge_t *e_P,
    psgr_dir_bit_t *db_P,
    psgr_edge_data_t **ed_P
  )
  { demand(a.ix != ANONE.ix, "arc is null");
    psgr_edge_t e = psgr_arc_edge(a);
    psgr_dir_bit_t db = psgr_arc_dir_bit(a);
    psgr_edge_data_t *ed = NULL;
    if (gr != NULL)
      { demand(e.ix < gr->NE, "invalid edge ID");
        ed = &(gr->edata[e.ix]);
      }
    if (e_P != NULL) { (*e_P) = e; }
    if (db_P != NULL) { (*db_P) = db; }
    if (ed_P != NULL) { (*ed_P) = ed; }
  }  

psgr_arc_t psgr_unsafe_add_unlinked_edge(psgr_t *gr, psgr_vertex_t org, psgr_vertex_t dst, psgr_path_t P)
  { demand(org.ix < gr->NV, "invalid vertex {org}");
    demand(dst.ix < gr->NV, "invalid vertex {dst}");
    if (gr->NE == gr->NE_max)
      { /* Reallocate the edge vector with more space: */ 
        gr->NE_max = gr->NE_max*2 + 5;
        gr->edata = retalloc(gr->edata, gr->NE_max, psgr_edge_data_t);
      }
    assert(gr->NE < gr->NE_max);

    psgr_edge_t e = (psgr_edge_t){ gr->NE };
    psgr_edge_data_t *ed = &(gr->edata[e.ix]);

    gr->NE++;

    ed->org[0] = org;
    ed->org[1] = dst;
    
    ed->emark = CLEAR;
    ed->delta = NAN;
    ed->weight = NAN;
    ed->path = P;

    psgr_arc_t a = psgr_orient_edge(e, 0);
    return a;
  }
