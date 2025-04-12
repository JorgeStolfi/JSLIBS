/* See {pst_gr.h} */
/* Last edited on 2025-03-15 00:14:26 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <filefmt.h>
#include <jsstring.h>
#include <affirm.h>

#include <pst_gr.h>

void pst_gr_split_arc
  ( pst_gr_t *gr,
    pst_gr_arc_t ai,
    pst_gr_edge_t *ei_P,
    pst_gr_dir_bit_t *db_P,
    pst_gr_edge_data_t **ed_P
  );
  /* Given an arc index {ai}, obtains the index {ei} of its underlying
    undirected edge, its direction bit {db}, and (if {gr} is not {NULL})
    the address of the corresponding edge data record {gr.edata[ei]}.
    These values are returned in {*ei_P}, {*db_P}, and {*ed_P}, if they
    are not {NULL}. Also checks that {ai} is not {NONE} and (if {gr} is
    not {NULL}) that {ei} is a valid edge index for {gr}. */

pst_gr_edge_data_t pst_gr_edge_data_create(void);
  /* Returns an edge data record with default information. */

void pst_gr_find_onext_oprev
  ( pst_gr_t *gr,
    pst_gr_vertex_t vi,
    r2_t *unew,
    pst_gr_arc_t *oprev_P,
    pst_gr_arc_t *onext_P
  );
  /* Among the ring of arcs out of vertex {vi}, finds the two arcs {oprev}, {onext}
    such that the vector {unew} lies between the starting directions {uprev,unext}
    of {oprev,onext} as they leave {vi}. Returns the results in {*oprev_P}
    and {*onext_P}.  Fails if {vi} has no outgoing edges. */

void pst_gr_set_org(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_vertex_t org);
  /* Sets the vertex with index {org} as the origin of the arc {ai}
    (which must not be {NULL}. */
 
void pst_gr_set_onext(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_arc_t bi);
  /* Sets the {enext} field of the edge of {ai} so that {pst_gr_arc_onext(ai)} becomes {bi}. 
    Other {enext} and {oprev} fields must be set too to keep the structure consistent. */
    
void pst_gr_set_oprev(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_arc_t bi);
  /* Sets the {enext} field of the edge of {ai} so that {pst_gr_arc_oprev(ai)} becomes {bi}. 
    Other {enext} and {oprev} fields must be set too to keep the structure consistent. */

#define NONE pst_gr_NONE

#define UNMARKED pst_gr_UNMARKED

/* IMPLEMENTATIONS */

pst_gr_t *pst_gr_new(uint32_t NV_max, uint32_t NE_max, uint32_t NX, uint32_t NY)
  { demand((NX >= 0) && (NX <= pst_gr_MAX_IMG_SIZE), "bad {NX}");
    demand((NY >= 0) && (NY <= pst_gr_MAX_IMG_SIZE), "bad {NY}");

    pst_gr_t *gr = talloc(1, pst_gr_t);
    gr->NE_max = NE_max;
    gr->NV_max = NV_max;
    gr->NE = 0;
    gr->NV = 0;
    
    gr->vdata = talloc(NV_max, pst_gr_vertex_data_t);
    gr->edata = talloc(NE_max, pst_gr_edge_data_t);
    
    gr->NX = NX;
    gr->NY = NY;
    if ((NX != 0) && (NY != 0))
      { uint32_t NXY = NX*NY;
        gr->vix = talloc(NXY, pst_gr_vertex_t);
        for (uint32_t xy = 0; xy < NXY; xy++) { gr->vix[xy] = NONE; }
      }
    else
      { demand((NX == 0) && (NY == 0), "inconsistent {NX,NY}");
        gr->vix = NULL; 
      }
    return gr;
  }

pst_gr_edge_data_t pst_gr_edge_data_create(void)
  { pst_gr_edge_data_t ed;
    ed.org[0] = ed.org[1] = UINT32_MAX;  /*Index to the vextex in the list*/
    ed.delta = 0; /*Derivative*/
    ed.weight = 0; /* weight*/
    ed.emark = 0; /*marking number*/
    ed.path = pst_gr_path_NULL;
    return ed;
  }

pst_gr_edge_t pst_gr_arc_edge(pst_gr_arc_t ai)
  { if (ai == NONE) { return NONE; }
    return (pst_gr_edge_t)(((uint32_t)ai) / 2);
  }

pst_gr_dir_bit_t pst_gr_arc_dir_bit(pst_gr_arc_t ai)
  { demand(ai != NONE, "arc is {NONE}");
    return (pst_gr_dir_bit_t)(((uint32_t)ai) % 2);
  }
  
pst_gr_arc_t pst_gr_orient_edge(pst_gr_edge_t ei, pst_gr_dir_bit_t db)
  { demand(ei != NONE, "edge is {NONE}");
    demand((db == 0) || (db == 1), "invalid arc direction bit");
    return (pst_gr_arc_t)(2*((uint32_t)ei) + db);
  }

pst_gr_edge_t pst_gr_arc_sym(pst_gr_arc_t ai)
  { demand(ai != NONE, "arc is {NONE}");
    pst_gr_edge_t ei = pst_gr_arc_edge(ai);
    pst_gr_dir_bit_t db = pst_gr_arc_dir_bit(ai);
    return pst_gr_orient_edge(ei, 1-db);
  }
  
void pst_gr_split_arc
  ( pst_gr_t *gr,
    pst_gr_arc_t ai,
    pst_gr_edge_t *ei_P,
    pst_gr_dir_bit_t *db_P,
    pst_gr_edge_data_t **ed_P
  )
  { demand(ai != NONE, "arc is {NONE}");
    pst_gr_edge_t ei = pst_gr_arc_edge(ai);
    pst_gr_dir_bit_t db = pst_gr_arc_dir_bit(ai);
    pst_gr_edge_data_t *ed = NULL;
    if (gr != NULL)
      { demand(ei < gr->NE, "invalid edge ID");
        ed = &(gr->edata[ei]);
      }
    if (ei_P != NULL) { (*ei_P) = ei; }
    if (db_P != NULL) { (*db_P) = db; }
    if (ed_P != NULL) { (*ed_P) = ed; }
  }  

pst_gr_vertex_t pst_gr_arc_org(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return ed->org[db];
  }

pst_gr_vertex_t pst_gr_arc_dst(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return ed->org[1-db];
  }

pst_gr_edge_t pst_gr_arc_onext(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return ed->enext[db];
  }
  
pst_gr_edge_t pst_gr_arc_oprev(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return ed->eprev[db];
  }

pst_gr_arc_t pst_gr_arc_lnext(pst_gr_t *gr, pst_gr_arc_t ai)
  { return pst_gr_arc_oprev(gr, pst_gr_arc_sym(ai));
  }

double pst_gr_arc_weight(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, NULL, &ed);
    return ed->weight;
  }

double pst_gr_arc_delta(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return (db == 0 ? +1.0 : -1.0)*ed->delta;
  }

pst_gr_path_t pst_gr_arc_path(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    return (db == 0 ? ed->path : pst_gr_path_reverse(ed->path));
  }

r2_t pst_gr_arc_start_dir(pst_gr_t *gr, pst_gr_arc_t ai)
  { pst_gr_vertex_t org = pst_gr_arc_org(gr, ai);
    r2_t *po = &(gr->vdata[org].coords);
    pst_gr_vertex_t dst = pst_gr_arc_dst(gr, ai);
    r2_t *pd = &(gr->vdata[dst].coords);
    pst_gr_path_t P = pst_gr_arc_path(gr, ai);
    return pst_gr_path_start_dir(po, P, pd);
  }

pst_gr_vertex_t pst_gr_add_vertex
  ( pst_gr_t *gr,
    int32_t x,
    int32_t y,
    r2_t coords
  )
  { if (gr->NV == gr->NV_max)
      { gr->NV_max = gr->NV_max*2 + 5;
        gr->vdata = retalloc(gr->vdata, gr->NV_max, pst_gr_vertex_data_t);
      }
    assert(gr->NV < gr->NV_max);
    pst_gr_vertex_t vi = gr->NV;
    pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);

    demand((x == -1) == (y == -1), "inconsistent {x,y}");
    vd->x = x;
    vd->y = y;
    if ((x != -1) && (y != -1))
      { demand((gr->NX != 0) && (gr->NY != 0), "graph has no associated map");
        demand((x >= 0) && (x < gr->NX), "invalid {x} coordinate");
        demand((y >= 0) && (y < gr->NY), "invalid {y} coordinate");
        uint32_t xy = gr->NX*(uint32_t)y + (uint32_t)x;
        gr->vix[xy] = vi;
      }
        
    vd->vmark = UNMARKED;
    vd->aout = NONE;
    vd->coords = coords;
    gr->NV++;
    return vi;
  }

pst_gr_arc_t pst_gr_add_edge
  ( pst_gr_t *gr,
    pst_gr_vertex_t org,  
    pst_gr_vertex_t dst,
    double d,
    double w,
    pst_gr_path_t P,
    bool_t setLinks
  )
  { demand(org < gr->NV, "invalid vertex {org}");
    demand(dst < gr->NV, "invalid vertex {dst}");
    if (gr->NE == gr->NE_max)
      { gr->NE_max = gr->NE_max*2 + 5;
        gr->edata = retalloc(gr->edata, gr->NE_max, pst_gr_edge_data_t);
      }
    assert(gr->NE < gr->NE_max);
     
    auto void set_enext_eprev(pst_gr_arc_t ai);
      /* Sets {ed.onext[db]} and {ei.oprev[db]} where {ed} is
        {gr.edgata[ei]}, {ei} is the unoriented edge index of arc {ai},
        and {db} is its direction bit. Assumes that arc {ai} has not been
        inserted yet in the ring of arcs out of the vertex 
        {pst_gr_arc_org(gr,ai)}. */

    pst_gr_edge_t ei = gr->NE;
    pst_gr_edge_data_t *edi = &(gr->edata[ei]);
    
    pst_gr_arc_t ai = 2*ei + 0;
    pst_gr_arc_t bi = pst_gr_arc_sym(ai);
    gr->NE++;

    edi->org[0] = org;
    edi->org[1] = dst;
    
    edi->emark = UNMARKED;
    edi->delta = d;
    edi->weight = w;
    edi->path = P;
    
    if (setLinks)
      { set_enext_eprev(ai);
        set_enext_eprev(bi);

        if (gr->vdata[org].aout == NONE) { gr->vdata[org].aout = ai; }
        if (gr->vdata[dst].aout == NONE) { gr->vdata[dst].aout = bi; }
      }

    return ai;
    
    void set_enext_eprev(pst_gr_arc_t ai)
      { pst_gr_edge_t ei = pst_gr_arc_edge(ai);
        pst_gr_edge_data_t *ed = &(gr->edata[ei]);
        pst_gr_dir_bit_t db = pst_gr_arc_dir_bit(ai);
        pst_gr_vertex_t ovi = pst_gr_arc_org(gr, ai);
        pst_gr_vertex_t dvi = pst_gr_arc_dst(gr, ai);
        
        pst_gr_arc_t aout = gr->vdata[ovi].aout;
        if (aout == NONE)
          { /* Arc {ai} will be first out of {ovi}: */
            ed->enext[db] = ed->eprev[db] = ai;
          }
        else
          { /* Origin vertex {ovi} has degree 1 or more, must find proper place: */
            r2_t *po = &(gr->vdata[ovi].coords);
            r2_t *pd = &(gr->vdata[dvi].coords);
            r2_t uai = pst_gr_path_start_dir(po, P, pd);
            /* Find position of {uai} in the ring of arcs around {ovi} */
            pst_gr_arc_t oprev, onext;
            pst_gr_find_onext_oprev(gr, ovi, &uai, &oprev, &onext);
            pst_gr_set_onext(gr, oprev, ai);
            pst_gr_set_onext(gr, ai, onext);
            pst_gr_set_oprev(gr, onext, ai);
            pst_gr_set_oprev(gr, ai, oprev);
          }
      }
  }
            
void pst_gr_find_onext_oprev
  ( pst_gr_t *gr,
    pst_gr_vertex_t vi,
    r2_t *unew,
    pst_gr_arc_t *oprev_P,
    pst_gr_arc_t *onext_P
  )
  { bool_t debug = FALSE;
    pst_gr_arc_t onext, oprev;
    pst_gr_arc_t aout = gr->vdata[vi].aout;
    if (pst_gr_arc_onext(gr, aout) == aout)
      { /* Origin vertex {ovi} has degree 1, insert anywhere: */
        oprev = aout; onext = aout; 
      }
    else
      { onext = aout; /* Next arc out of {vi} to check. */
        oprev = pst_gr_arc_oprev(gr, aout); /* Previous arc. */
        sign_t sense = 0; /* Set to {+1} when {unew} is between {uprev} and {unext} in CCW order. */
        do
          { r2_t uprev = pst_gr_arc_start_dir(gr, oprev);
            r2_t unext = pst_gr_arc_start_dir(gr, onext);
            sense = r2_cyclic_order(&uprev, unew, &unext);
            if (debug)
              { r2_gen_print(stderr, &uprev, "%8.4f", "uprev = ( ", " ", " ) ");
                pst_gr_arc_print(stderr, gr, oprev); fputc('\n', stderr);
                r2_gen_print(stderr, unew,   "%8.4f", "unew =  ( ", " ", " )\n");
                r2_gen_print(stderr, &unext, "%8.4f", "unext = ( ", " ", " )  ");
                pst_gr_arc_print(stderr, gr, onext); fputc('\n', stderr);
                fprintf(stderr, "sense = %+d\n", sense);
                fprintf(stderr, "\n");
              }
            demand(sense != 0, "collinear arcs");
            if (sense == +1) { break; }
            oprev = onext; 
            onext = pst_gr_arc_onext(gr, onext);
          }
        while (onext != aout);      
        assert(oprev != NONE);
        affirm(sense == +1, "{r2_dir_sense} failed");
      }
    (*oprev_P) = oprev;
    (*onext_P) = onext;
   }

uint32_t pst_gr_outdegree(pst_gr_t *gr, pst_gr_vertex_t vi)
  { demand(vi < gr->NV, "invalid vertex index {vi}");
    pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
    if (vd->aout == NONE) return 0;
    uint32_t count = 1;
    pst_gr_arc_t ai = vd->aout;
    for (ai = pst_gr_arc_onext(gr, vd->aout); ai!= vd->aout; ai = pst_gr_arc_onext(gr, ai)) { count++; }
    return count;
  }

pst_gr_arc_t pst_gr_get_connecting_arc(pst_gr_t *gr, pst_gr_vertex_t vi0, pst_gr_vertex_t vi1)
  { demand (vi0 < gr->NV, "invalid vertex index {vi0}");
    demand (vi1 < gr->NV, "invalid vertex index {vi1}");
    pst_gr_vertex_data_t *vd0 = &(gr->vdata[vi0]);
    pst_gr_vertex_data_t *vd1 = &(gr->vdata[vi1]);
    if ((vd0->aout == NONE) || (vd1->aout == NONE))
      { return NONE; }
    pst_gr_arc_t ai = vd0->aout;
    do
      { pst_gr_vertex_t dst = pst_gr_arc_dst(gr, ai);
        if (dst == vi1) { return ai; }
        ai = pst_gr_arc_onext(gr, ai);
      } while(ai != vd0->aout);
    return NONE;
  }

void pst_gr_check_consistency(pst_gr_t *gr)
  {
    demand(gr->NV <= gr->NV_max, "NV > NV_MAX"); 
    demand(gr->NE <= gr->NE_max, "NE > NE_MAX"); 
    
    demand((gr->NX >= 0) && (gr->NX <= pst_gr_MAX_IMG_SIZE), "bad NX");
    demand((gr->NY >= 0) && (gr->NY <= pst_gr_MAX_IMG_SIZE), "bad NY");
    
    auto void check_vertex(pst_gr_vertex_t vi); 
      /* Runs some consistency checks on the vertex with index {vi} (which must not be {NONE})
        of the graph {gr}.  Bombs out with message if any test fails. */

    auto void check_edge(pst_gr_edge_t ei);
      /* Runs some consistency checks on the edgex with index {ei} (which must not be {NONE})
        of the graph {gr}.  Bombs out with message if any test fails. */

    auto void check_vix(void);
      /* Checks whether the {vix} table is consistent with the indices {x,y}
        in the vertex records.  Assumes that these coordinates have been 
        already checked to be either {-1,-1} or in the ranges {0..gr.NX-1},
        {0..gr.NY-1}. */
    
    for (pst_gr_vertex_t vi = 0; vi < gr->NV; vi++)
      { check_vertex(vi); }
    
    for (pst_gr_edge_t ei = 0; ei < gr->NE; ei++)
      { check_edge(ei); }
      
    check_vix();
    
    return;

    void check_vertex(pst_gr_vertex_t vi)  
      { demand(vi != NONE, "vertex index is {NONE}");
        pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
        demand((vd->x == -1) == (vd->y == -1), "inconsistent grid indices");
        if ((vd->x != -1) || (vd->y != -1))
          { demand((gr->NX != 0) && (gr->NY != 0), "graph has no associated map");
            demand((vd->x >= 0) && (vd->x < gr->NX), "invalid {x} coordinate");
            demand((vd->y >= 0) && (vd->y < gr->NY), "invalid {y} coordinate");
            uint32_t xy = gr->NX*(uint32_t)(vd->y) + (uint32_t)(vd->x);
            demand(gr->vix[xy] == vi, "inconsistent pixel-to-vertex table {vix}");
          }
        pst_gr_edge_t ai = vd->aout;
        if (ai != NONE)
          { demand(pst_gr_arc_org(gr, ai) == vi, "inconsistent outgoing arc"); }
      }

    void check_edge(pst_gr_edge_t ei) 
      { demand(ei != NONE, "edge index is {NONE}");
        pst_gr_edge_data_t *ed = &(gr->edata[ei]);
        demand(isfinite(ed->weight) && (ed->weight >= 0), "invalid edge weight");
        demand(isfinite(ed->delta), "invalid edge delta");
        for (uint32_t db = 0; db <= 1; db++)
          { pst_gr_arc_t ai = pst_gr_orient_edge(ei, (pst_gr_dir_bit_t)db);
            pst_gr_vertex_t vi = ed->org[db];
            demand(vi < gr->NV, "edge {e} has invalid origin or destination");
            demand(ed->enext[db] != NONE, "invalid null {enext}");
            demand(ed->eprev[db] != NONE, "invalid null {enext}");
            demand(pst_gr_arc_onext(gr, pst_gr_arc_oprev(gr, ai)) == ai, "inconsistent {enext(eprev(ai))}");
            demand(pst_gr_arc_oprev(gr, pst_gr_arc_onext(gr, ai)) == ai, "inconsistent {eprev(enext(ai))}");
          }
        /* !!! should check {ed->path} !!! */
      }
      
    void check_vix(void)
      { if ((gr->NX != 0) && (gr->NY != 0))
         { demand(gr->vix != NULL, "missing pixel-to-vertex table {vix}");
           uint32_t NXY = gr->NX*gr->NY;
           for (uint32_t xy = 0; xy < NXY; xy++)
             { pst_gr_vertex_t vi = gr->vix[xy];
               if (vi != NONE)
                 { pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
                   int32_t x = vd->x, y = vd->y;
                   demand((x != -1) && (y != -1), "spurious {vix} entry");
                   assert((x >= 0) && (x < gr->NX)); /* Should have been checked. */
                   assert((y >= 0) && (y < gr->NY)); /* Should have been checked. */
                   demand(xy == (gr->NX*(uint32_t)y + (uint32_t)x), "inconsistent {vix} entry");
                 }
             }
         }
      }
  }

bool_t pst_gr_equal(pst_gr_t *gr, pst_gr_t *hr)
  {
    bool_t ok = TRUE;
    
    auto void compare_int32(int32_t gval, int32_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_uint32(uint32_t gval, uint32_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_r2(r2_t *gval, r2_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_edge_data(pst_gr_edge_data_t *gval, pst_gr_edge_data_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_vertex_data(pst_gr_vertex_data_t *gval, pst_gr_vertex_data_t *hval, char *name1, int32_t ix, char *name2);
    
    compare_uint32(gr->NV, hr->NV, "NV", -1, NULL);
    compare_uint32(gr->NE, hr->NE, "NE", -1, NULL);
    
    compare_uint32(gr->NX, hr->NX, "NX", -1, NULL);
    compare_uint32(gr->NY, hr->NY, "NY", -1, NULL);
    
    for (int32_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_vertex_data_t *gvd = &(gr->vdata[vi]);
        pst_gr_vertex_data_t *hvd = &(hr->vdata[vi]);
        compare_vertex_data(gvd, hvd, "vertex", vi, NULL);
      }
      
    for (int32_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_edge_data_t *ged = &(gr->edata[ei]);
        pst_gr_edge_data_t *hed = &(hr->edata[ei]);
        compare_edge_data(ged, hed, "edge", ei, NULL);
      }
      
    uint32_t NXY = gr->NX*gr->NY;
    for (int32_t xy = 0; xy < NXY; xy++)
      { compare_uint32(gr->vix[xy], hr->vix[xy], "vix", xy, NULL); }
    
    return ok;
    
    void blast(char *name1, int32_t ix, char *name2)
      { fprintf(stderr, "** %s", name1);
        if (ix >= 0) { fprintf(stderr, " at index %d", ix); }
        if (name2 != NULL) { fprintf(stderr, " %s", name2); }
        fprintf(stderr, " differ:");
      }
    
    void compare_int32(int32_t gval, int32_t hval, char *name1, int32_t ix, char *name2)
      { if (gval != hval)
          { blast(name1, ix, name2);
            fprintf(stderr, " %d %d\n", gval, hval);
            ok = FALSE;
          }
      }
    
    void compare_uint32(uint32_t gval, uint32_t hval, char *name1, int32_t ix, char *name2)
      { if (gval != hval)
          { blast(name1, ix, name2);
            fprintf(stderr, " %d %d\n", gval, hval);
            ok = FALSE;
          }
      }
    
    void compare_double(double gval, double hval, char *name1, int32_t ix, char *name2)
      { double tol = 2.0e-6; /* Considering the format of {pst_gr_write}. */
        if (fabs(gval - hval) > tol)
          { blast(name1, ix, name2);
            fprintf(stderr, " %24.16e %24.16e\n", gval, hval);
            ok = FALSE;
          }
      }
      
    void compare_r2(r2_t *gval, r2_t *hval, char *name1, int32_t ix, char *name2)
      { double d = r2_dist(gval, hval);
        double tol = 2.0e-6; /* Considering the format of {pst_gr_write}. */
        if (d > tol)
          { blast(name1, ix, name2);
            r2_print(stderr, gval);
            fprintf(stderr, " ");
            r2_print(stderr, hval);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
      
    void compare_edge_data(pst_gr_edge_data_t *gval, pst_gr_edge_data_t *hval, char *name1, int32_t ix, char *name2)
      { for (uint32_t db = 0; db <= 1; db++)
          { compare_uint32(gval->org[db], hval->org[db], name1, ix, "{org[db]}");
            compare_uint32(gval->enext[db], hval->enext[db], name1, ix, "{enext[db]}");
            compare_uint32(gval->eprev[db], hval->eprev[db], name1, ix, "{eprev[db]}");
          }
        compare_double(gval->delta, hval->delta, name1, ix, "{delta}");
        compare_double(gval->weight, hval->weight, name1, ix, "{weight}");
        /* !!! Should compare the paths. !!! */
      }
      
    void compare_vertex_data(pst_gr_vertex_data_t *gval, pst_gr_vertex_data_t *hval, char *name1, int32_t ix, char *name2)
      { 
        compare_int32(gval->x, hval->x, "vertex", ix, "{x}");
        compare_int32(gval->y, hval->y, "vertex", ix, "{y}");
        /* Don't care about the {mark}. */
        compare_r2(&(gval->coords), &(hval->coords), "vertex", ix, "{coords}");
        compare_uint32(gval->aout, hval->aout, "edge of vertex", ix, NULL);
      }
      
  }

pst_gr_arc_t pst_gr_find_enclosing_face(pst_gr_t *gr, int32_t x, int32_t y)
  { r2_t p = (r2_t){{ x + 0.5, y + 0.5 }};
    /* Tries to locate the face that contains {p}. */
    /* !!! Should do this right !!! */
    uint32_t vi = pst_gr_find_nearest_vertex(gr, &p);
    pst_gr_vertex_data_t *v = &(gr->vdata[vi]);
    assert(vi < gr->NV);
    pst_gr_arc_t a = v->aout;
    return a;
  }

pst_gr_vertex_t pst_gr_find_nearest_vertex(pst_gr_t *gr, r2_t *p)
  { double min_dist_sqr = +INF;
    pst_gr_vertex_t min_vi = UINT32_MAX;
    for (pst_gr_vertex_t vi = 0; vi < gr->NV; vi++)
      { pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
        double d = r2_dist_sqr(p, &(vd->coords));
        if (d < min_dist_sqr) { min_dist_sqr = d; min_vi = vi; }
      }
    demand(min_vi != UINT32_MAX, "graph has no valid vertices");
    return min_vi;
  }

double pst_gr_compute_left_face_curl(pst_gr_t *gr, pst_gr_arc_t ai)
  { double curl = 0;
    pst_gr_arc_t bi = ai;
    do
      { curl += pst_gr_arc_delta(gr, bi);
        bi = pst_gr_arc_lnext(gr, bi);
      } while (bi != ai);
    return curl;
  }

r2_t pst_gr_left_face_barycenter(pst_gr_t *gr, pst_gr_arc_t ai)
  {
    r2_t sum_p = (r2_t) {{0,0}};
    uint32_t count = 0;
    pst_gr_arc_t bi = ai;
    do 
      { pst_gr_vertex_t borg = pst_gr_arc_org(gr, bi);
        assert(borg < gr->NV);
        pst_gr_vertex_data_t *vd = &(gr->vdata[borg]);
        r2_add(&sum_p, &(vd->coords), &sum_p);
        count ++;
        bi = pst_gr_arc_lnext(gr, bi);
      } while(bi != ai);
      
    r2_t bar; r2_scale(1.0/(double)count, &sum_p, &bar);
    return bar;
  }

pst_gr_vertex_t pst_gr_get_vertex_from_map_indices(pst_gr_t *gr, int32_t x, int32_t y)
  { demand((gr->NX > 0) && (gr->NY > 0), "graph has no associated map");
    demand((x >= 0) && (x < gr->NX), "invalid map index {x}");
    demand((y >= 0) && (y < gr->NY), "invalid map index {y}");
    uint32_t xy = gr->NX*(uint32_t)y + (uint32_t)x;
    return gr->vix[xy];
  }

void pst_gr_get_map_indices_from_vertex(pst_gr_t *gr, pst_gr_vertex_t vi, int32_t *x_P, int32_t *y_P)
  { demand((gr->NX > 0) && (gr->NY > 0), "graph has no associated map");
    demand((vi != NONE) && (vi >= 0) && (vi < gr->NV), "invalid vertex index {vi}");
    pst_gr_vertex_data_t *vd = &(gr->vdata[vi]);
    (*x_P) = vd->x;
    (*y_P) = vd->y;
  }

void pst_gr_free(pst_gr_t *gr)
  {
    /* Free the edge-data records and the edge list: */
    for (int32_t ei = 0; ei < gr->NE; ei++)
      { pst_gr_edge_data_t *ed = &(gr->edata[ei]);
        pst_gr_path_free(ed->path);
      }
    free(gr->edata);
    free(gr->vdata);
    free(gr);
  }

void pst_gr_arc_print(FILE *wr, pst_gr_t *gr, pst_gr_arc_t ai)
  { if (ai == NONE)
      { fprintf(wr, "-1"); }
    else
      { pst_gr_dir_bit_t db; pst_gr_edge_t ei;
        pst_gr_split_arc(gr, ai, &ei, &db, NULL);
        fprintf(wr, "%d:%d", ei, db);
      }
  }
  
void pst_gr_set_onext(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_arc_t bi)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    ed->enext[db] = bi;
  }
   
void pst_gr_set_oprev(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_arc_t bi)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    ed->eprev[db] = bi;
  }

void pst_gr_set_org(pst_gr_t *gr, pst_gr_arc_t ai, pst_gr_vertex_t org)
  { pst_gr_dir_bit_t db; pst_gr_edge_data_t *ed;
    pst_gr_split_arc(gr, ai, NULL, &db, &ed);
    ed->org[db] = org;
  }
