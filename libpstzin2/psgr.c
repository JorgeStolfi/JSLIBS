/* See {psgr.h} */
/* Last edited on 2025-04-27 12:11:31 by stolfi */
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
#include <r2.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr_unsafe.h>

#include <psgr.h>

/* INTERNAL PROTOTYPES */

psgr_edge_data_t psgr_edge_data_create(void);
  /* Returns an edge data record with default information. */

void psgr_find_onext_oprev
  ( psgr_t *gr,
    psgr_vertex_t v,
    r2_t *unew,
    psgr_arc_t *oprev_P,
    psgr_arc_t *onext_P
  );
  /* Among the ring of arcs out of vertex {v}, finds the two arcs {oprev}, {onext}
    such that the vector {unew} lies between the starting directions {uprev,unext}
    of {oprev,onext} as they leave {v}. Returns the results in {*oprev_P}
    and {*onext_P}.  Fails if {v} has no outgoing edges.
    
    This option requires that the {.enext} and {.eprev} links of all
    edges incident to {v} proprly describe the onext ring of {v}. */

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

/* IMPLEMENTATIONS */

psgr_t *psgr_new(uint32_t NV_max, uint32_t NE_max, uint32_t NX, uint32_t NY)
  { demand((NV_max <= psgr_MAX_VERTICES), "bad {NV_max}");
    demand((NE_max <= psgr_MAX_EDGES), "bad {NE_max}");
    
    demand((NX > 0) && (NX <= psgr_MAX_IMG_SIZE), "bad {NX}");
    demand((NY > 0) && (NY <= psgr_MAX_IMG_SIZE), "bad {NY}");

    psgr_t *gr = talloc(1, psgr_t);
    gr->NE_max = NE_max;
    gr->NV_max = NV_max;
    gr->NE = 0;
    gr->NV = 0;
    
    gr->vdata = talloc(NV_max, psgr_vertex_data_t);
    gr->edata = talloc(NE_max, psgr_edge_data_t);
    
    gr->NX = NX;
    gr->NY = NY;
    
    uint32_t NXY = NX*NY;
    gr->vix = talloc(NXY, psgr_vertex_t);
    for (uint32_t xy = 0; xy < NXY; xy++) { gr->vix[xy] = VNONE; }
    return gr;
  }

void psgr_free(psgr_t *gr)
  {
    /* Free the edge-data records and the edge list: */
    for (uint32_t e_ix  = 0; e_ix < gr->NE; e_ix++)
      { psgr_edge_data_t *ed = &(gr->edata[e_ix]);
        psgr_path_free(ed->path);
      }
    free(gr->edata);
    free(gr->vdata);
    free(gr);
  }
uint32_t psgr_num_vertices(psgr_t *gr)
  { return gr->NV; }
    
uint32_t psgr_num_edges(psgr_t *gr)
  { return gr->NE; }

i2_t psgr_image_size(psgr_t *gr)
  {
    i2_t sz = (i2_t){{ (int32_t)gr->NX, (int32_t)gr->NY }};
    return sz;
  }
  
psgr_edge_data_t psgr_edge_data_create(void)
  { psgr_edge_data_t ed;
    ed.org[0] = ed.org[1] = VNONE;  /*Index to the vextex in the list*/
    ed.delta = 0; /*Derivative*/
    ed.weight = 0; /* weight*/
    ed.emark = 0; /*marking number*/
    ed.path = psgr_path_NULL;
    return ed;
  }

psgr_edge_t psgr_arc_edge(psgr_arc_t a)
  { if (a.ix == ANONE.ix) { return ENONE; }
    return (psgr_edge_t){ ((uint32_t)a.ix) / 2 };
  }

psgr_dir_bit_t psgr_arc_dir_bit(psgr_arc_t a)
  { demand(a.ix != ANONE.ix, "arc is null");
    return (psgr_dir_bit_t)(((uint32_t)a.ix) % 2);
  }
  
psgr_arc_t psgr_orient_edge(psgr_edge_t e, psgr_dir_bit_t db)
  { demand(e.ix != ENONE.ix, "edge is null");
    demand((db == 0) || (db == 1), "invalid arc direction bit");
    return (psgr_arc_t){ 2*((uint32_t)e.ix) + db };
  }

psgr_arc_t psgr_arc_sym(psgr_arc_t a)
  { demand(a.ix != ANONE.ix, "arc is null");
    psgr_edge_t e = psgr_arc_edge(a);
    psgr_dir_bit_t db = psgr_arc_dir_bit(a);
    return psgr_orient_edge(e, 1-db);
  }
  
void psgr_arc_split
  ( psgr_arc_t a,
    psgr_edge_t *e_P,
    psgr_dir_bit_t *db_P
  )
  { demand(a.ix != ANONE.ix, "arc is null");
    psgr_edge_t e = psgr_arc_edge(a);
    psgr_dir_bit_t db = psgr_arc_dir_bit(a);
    if (e_P != NULL) { (*e_P) = e; }
    if (db_P != NULL) { (*db_P) = db; }
  }  

psgr_vertex_t psgr_arc_org(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return ed->org[db];
  }

psgr_vertex_t psgr_arc_dst(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return ed->org[1-db];
  }

psgr_arc_t psgr_arc_onext(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return ed->enext[db];
  }
  
psgr_arc_t psgr_arc_oprev(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return ed->eprev[db];
  }

psgr_arc_t psgr_arc_dnext(psgr_t *gr, psgr_arc_t a)
  { return psgr_arc_sym(psgr_arc_onext(gr,psgr_arc_sym(a))); }

psgr_arc_t psgr_arc_dprev(psgr_t *gr, psgr_arc_t a)
  { return psgr_arc_sym(psgr_arc_oprev(gr,psgr_arc_sym(a))); }

psgr_arc_t psgr_arc_lnext(psgr_t *gr, psgr_arc_t a)
  { return psgr_arc_oprev(gr, psgr_arc_sym(a)); }

psgr_arc_t psgr_arc_lprev(psgr_t *gr, psgr_arc_t a)
  { return psgr_arc_sym(psgr_arc_onext(gr, a)); }

double psgr_arc_delta(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return (db == 0 ? +1.0 : -1.0)*ed->delta;
  }

void psgr_arc_set_delta(psgr_t *gr, psgr_arc_t a, double d)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    ed->delta = (db == 0 ? +1.0 : -1.0)*d;
  }

double psgr_arc_weight(psgr_t *gr, psgr_arc_t a)
  { psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, NULL, &ed);
    return ed->weight;
  }

void psgr_arc_set_weight(psgr_t *gr, psgr_arc_t a, double w)
  { psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, NULL, &ed);
    ed->weight = w;
  }

i2_t psgr_vertex_pixel(psgr_t *gr, psgr_vertex_t v)
  { assert((gr->NX > 0) && (gr->NY > 0));
    demand((v.ix != VNONE.ix) && (v.ix  >= 0) && (v.ix  < gr->NV), "invalid vertex index {v}");
    return gr->vdata[v.ix].pixel; 
  }

psgr_path_t psgr_arc_path(psgr_t *gr, psgr_arc_t a)
  { psgr_dir_bit_t db; psgr_edge_data_t *ed;
    psgr_arc_unsafe_split(gr, a, NULL, &db, &ed);
    return (db == 0 ? ed->path : psgr_path_reverse(ed->path));
  }

psgr_mark_t psgr_vertex_mark(psgr_t *gr, psgr_vertex_t v)
  { return gr->vdata[v.ix].vmark; }
  
psgr_mark_t psgr_edge_mark(psgr_t *gr, psgr_edge_t e)
  { return gr->edata[e.ix].emark; }

void psgr_vertex_set_mark(psgr_t *gr, psgr_vertex_t v, psgr_mark_t mrk)
  { gr->vdata[v.ix].vmark = mrk; }
  
void psgr_edge_set_mark(psgr_t *gr, psgr_edge_t e, psgr_mark_t mrk)
  { gr->edata[e.ix].emark = mrk; }

r2_t psgr_arc_starting_dir(psgr_t *gr, psgr_arc_t a)
  { psgr_vertex_t org = psgr_arc_org(gr, a);
    i2_t *po = &(gr->vdata[org.ix].pixel);
    psgr_vertex_t dst = psgr_arc_dst(gr, a);
    i2_t *pd = &(gr->vdata[dst.ix].pixel);
    psgr_path_t P = psgr_arc_path(gr, a);
    return psgr_path_starting_dir(po, P, pd);
  }

psgr_vertex_t psgr_add_vertex(psgr_t *gr, i2_t *ip)
  { /* Resize the vertex table if needed: */
    if (gr->NV == gr->NV_max)
      { gr->NV_max = gr->NV_max*2 + 5;
        gr->vdata = retalloc(gr->vdata, gr->NV_max, psgr_vertex_data_t);
      }
    assert(gr->NV < gr->NV_max);
    psgr_vertex_t v = (psgr_vertex_t){ gr->NV };
    psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
    int32_t x = ip->c[0], y = ip->c[1];
    demand((gr->NX != 0) && (gr->NY != 0), "graph has no associated map");
    demand((x >= 0) && (x < gr->NX), "invalid {x} pixel index");
    demand((y >= 0) && (y < gr->NY), "invalid {y} pixel index");
    vd->pixel = *ip;
    vd->vmark = CLEAR;
    vd->aout = ANONE;
    gr->NV++;

    uint32_t xy = gr->NX*(uint32_t)y + (uint32_t)x;
    gr->vix[xy] = v;
        
    return v;
  }

psgr_arc_t psgr_add_edge
  ( psgr_t *gr,
    psgr_vertex_t org,  
    psgr_vertex_t dst,
    double d,
    double w,
    psgr_path_t P
  )
  { 
    psgr_arc_t a = psgr_unsafe_add_unlinked_edge(gr, org, dst, P);
    psgr_arc_t b = psgr_arc_sym(a);
    
    auto void set_onext_oprev(psgr_arc_t a);
      /* Sets {ed.onext[db]} and {e.oprev[db]} where {ed} is
        {gr.edgata[e]}, {e} is the unoriented edge index of arc {a},
        and {db} is its direction bit. Assumes that arc {a} has not been
        inserted yet in the ring of arcs out of the vertex 
        {psgr_arc_org(gr,a)}. */

    psgr_arc_set_delta(gr, a, d);
    psgr_arc_set_weight(gr, a, w);
    
    set_onext_oprev(a);
    set_onext_oprev(b);
    
    /* Set the {.aout} of {org} and {dst} if needed: */
    if (psgr_vertex_out_arc(gr, org).ix == ANONE.ix)
      { psgr_unsafe_vertex_set_out_arc(gr, org, a); }
      
    if (psgr_vertex_out_arc(gr, dst).ix == ANONE.ix)
      { psgr_unsafe_vertex_set_out_arc(gr, dst, b); }

    return a;
    
    void set_onext_oprev(psgr_arc_t a)
      { psgr_vertex_t ovi = psgr_arc_org(gr, a);
        psgr_vertex_t dvi = psgr_arc_dst(gr, a);
        
        psgr_arc_t aout = gr->vdata[ovi.ix].aout;
        if (aout.ix == ANONE.ix)
          { /* Arc {a} will be first out of {ovi}: */
            psgr_unsafe_arc_set_onext(gr, a, a);
            psgr_unsafe_arc_set_oprev(gr, a, a);
          }
        else
          { /* Origin vertex {ovi} has degree 1 or more, must find proper place: */
            i2_t *ipo = &(gr->vdata[ovi.ix].pixel);
            i2_t *ipd = &(gr->vdata[dvi.ix].pixel);
            r2_t uai = psgr_path_starting_dir(ipo, P, ipd);
            /* Find position of {uai} in the ring of arcs around {ovi} */
            psgr_arc_t oprev, onext;
            psgr_find_onext_oprev(gr, ovi, &uai, &oprev, &onext);
            psgr_unsafe_arc_set_onext(gr, oprev, a);
            psgr_unsafe_arc_set_onext(gr, a, onext);
            psgr_unsafe_arc_set_oprev(gr, onext, a);
            psgr_unsafe_arc_set_oprev(gr, a, oprev);
          }
      }
  }
            
void psgr_find_onext_oprev
  ( psgr_t *gr,
    psgr_vertex_t v,
    r2_t *unew,
    psgr_arc_t *oprev_P,
    psgr_arc_t *onext_P
  )
  { bool_t debug = FALSE;
    psgr_arc_t onext, oprev;
    psgr_arc_t aout = gr->vdata[v.ix].aout;
    if (psgr_arc_onext(gr, aout).ix == aout.ix)
      { /* Origin vertex {ovi} has degree 1, insert anywhere: */
        oprev = aout; onext = aout; 
      }
    else
      { onext = aout; /* Next arc out of {v} to check. */
        oprev = psgr_arc_oprev(gr, aout); /* Previous arc. */
        sign_t sense = 0; /* Set to {+1} when {unew} is between {uprev} and {unext} in CCW order. */
        do
          { r2_t uprev = psgr_arc_starting_dir(gr, oprev);
            r2_t unext = psgr_arc_starting_dir(gr, onext);
            sense = r2_cyclic_order(&uprev, unew, &unext);
            if (debug)
              { r2_gen_print(stderr, &uprev, "%8.4f", "uprev = ( ", " ", " ) ");
                psgr_arc_print(stderr, gr, oprev); fputc('\n', stderr);
                r2_gen_print(stderr, unew,   "%8.4f", "unew =  ( ", " ", " )\n");
                r2_gen_print(stderr, &unext, "%8.4f", "unext = ( ", " ", " )  ");
                psgr_arc_print(stderr, gr, onext); fputc('\n', stderr);
                fprintf(stderr, "sense = %+d\n", sense);
                fprintf(stderr, "\n");
              }
            demand(sense != 0, "collinear arcs");
            if (sense == +1) { break; }
            oprev = onext; 
            onext = psgr_arc_onext(gr, onext);
          }
        while (onext.ix != aout.ix);      
        assert(oprev.ix != ANONE.ix);
        affirm(sense == +1, "{r2_dir_sense} failed");
      }
    (*oprev_P) = oprev;
    (*onext_P) = onext;
   }

psgr_arc_t psgr_vertex_out_arc(psgr_t *gr, psgr_vertex_t v)
  { return gr->vdata[v.ix].aout; }

uint32_t psgr_vertex_outdegree(psgr_t *gr, psgr_vertex_t v)
  { demand(v.ix != VNONE.ix, "invalid null vertex");
    demand(v.ix < gr->NV, "invalid vertex index");
    psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
    if (vd->aout.ix == ANONE.ix) return 0;
    uint32_t count = 1;
    psgr_arc_t a = vd->aout;
    for (a = psgr_arc_onext(gr, vd->aout); a.ix != vd->aout.ix; a = psgr_arc_onext(gr, a)) { count++; }
    return count;
  }

psgr_arc_t psgr_get_connecting_arc(psgr_t *gr, psgr_vertex_t vi0, psgr_vertex_t vi1)
  { demand(vi0.ix != VNONE.ix, "invalid null vertex {vi0}");
    demand (vi0.ix < gr->NV, "invalid vertex index {vi0}");
    demand(vi1.ix != VNONE.ix, "invalid null vertex {vi1}");
    demand (vi1.ix < gr->NV, "invalid vertex index {vi1}");
    psgr_vertex_data_t *vd0 = &(gr->vdata[vi0.ix]);
    psgr_vertex_data_t *vd1 = &(gr->vdata[vi1.ix]);
    if ((vd0->aout.ix == ANONE.ix) || (vd1->aout.ix == ANONE.ix))
      { return ANONE; }
    psgr_arc_t a = vd0->aout;
    do
      { psgr_vertex_t dst = psgr_arc_dst(gr, a);
        if (dst.ix == vi1.ix) { return a; }
        a = psgr_arc_onext(gr, a);
      } while(a.ix != vd0->aout.ix);
    return ANONE;
  }

void psgr_check_consistency(psgr_t *gr)
  {
    demand(gr->NV <= gr->NV_max, "NV > NV_MAX"); 
    demand(gr->NE <= gr->NE_max, "NE > NE_MAX"); 
    
    demand((gr->NX >= 0) && (gr->NX <= psgr_MAX_IMG_SIZE), "bad NX");
    demand((gr->NY >= 0) && (gr->NY <= psgr_MAX_IMG_SIZE), "bad NY");
    
    auto void check_vertex(psgr_vertex_t v); 
      /* Runs some consistency checks on the vertex with index {v} (which must not be {VNONE})
        of the graph {gr}.  Bombs out with message if any test fails. */

    auto void check_edge(psgr_edge_t e);
      /* Runs some consistency checks on the edgex with index {e} (which must not be {ENONE})
        of the graph {gr}.  Bombs out with message if any test fails. */

    auto void check_vix(void);
      /* Checks whether the {vix} table is consistent with the pixel indices {x,y}
        in the vertex records.  Assumes that these indices have been 
        already checked to be either {-1,-1} or in the ranges {0..gr.NX-1},
        {0..gr.NY-1}. */
    
    for (uint32_t v = 0; v < gr->NV; v++)
      { check_vertex((psgr_vertex_t){ v }); }
    
    for (uint32_t e = 0; e < gr->NE; e++)
      { check_edge((psgr_edge_t){ e }); }
      
    check_vix();
    
    return;

    void check_vertex(psgr_vertex_t v)  
      { demand(v.ix != VNONE.ix, "vertex index is {NONE}");
        psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
        demand((gr->NX != 0) && (gr->NY != 0), "graph has no associated map");
        demand((vd->pixel.c[0] >= 0) && (vd->pixel.c[0] < gr->NX), "invalid {x} pixel index");
        demand((vd->pixel.c[1] >= 0) && (vd->pixel.c[1] < gr->NY), "invalid {y} pixel index");
        psgr_mark_t mkv = vd->vmark;
        demand((mkv == CLEAR) || (mkv == KEEP) || (mkv == KILL) || (mkv == DELETED), "invalid vertex mark");
        uint32_t xy = gr->NX*(uint32_t)(vd->pixel.c[1]) + (uint32_t)(vd->pixel.c[0]);
        demand(gr->vix[xy].ix == v.ix, "inconsistent pixel-to-vertex table {vix}");
        psgr_arc_t a = vd->aout;
        if (a.ix != ANONE.ix)
          { demand(psgr_arc_org(gr, a).ix == v.ix, "inconsistent outgoing arc"); }
      }

    void check_edge(psgr_edge_t e) 
      { demand(e.ix != ENONE.ix, "edge is null");
        psgr_edge_data_t *ed = &(gr->edata[e.ix]);
        demand(isfinite(ed->weight) && (ed->weight >= 0), "invalid edge weight");
        demand(isfinite(ed->delta), "invalid edge delta");
        psgr_mark_t mke = ed->emark;
        demand((mke == CLEAR) || (mke == DELETED), "invalid edge mark");
        for (uint32_t db = 0; db <= 1; db++)
          { psgr_arc_t a = psgr_orient_edge(e, (psgr_dir_bit_t)db);
            psgr_vertex_t v = ed->org[db];
            psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
            demand(v.ix < gr->NV, "edge {e} has invalid origin or destination");
            demand(ed->enext[db].ix != ANONE.ix, "invalid null {enext}");
            demand(ed->eprev[db].ix != ANONE.ix, "invalid null {eprev}");
            demand(psgr_arc_onext(gr, psgr_arc_oprev(gr, a)).ix == a.ix, "inconsistent {enext(eprev(a))}");
            demand(psgr_arc_oprev(gr, psgr_arc_onext(gr, a)).ix == a.ix, "inconsistent {eprev(enext(a))}");
            psgr_mark_t mkv = vd->vmark;
            if (mkv == DELETED)
              { demand(mke == DELETED, "non-deleted edge incident to deleted vertex"); }
          }
        /* !!! should check {ed->path} !!! */
      }
      
    void check_vix(void)
      { if ((gr->NX != 0) && (gr->NY != 0))
         { demand(gr->vix != NULL, "missing pixel-to-vertex table {vix}");
           uint32_t NXY = gr->NX*gr->NY;
           for (uint32_t xy = 0; xy < NXY; xy++)
             { psgr_vertex_t v = gr->vix[xy];
               if (v.ix != VNONE.ix)
                 { psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
                   int32_t x = vd->pixel.c[0], y = vd->pixel.c[1];
                   demand((x != -1) && (y != -1), "spurious {vix} entry");
                   assert((x >= 0) && (x < gr->NX)); /* Should have been checked. */
                   assert((y >= 0) && (y < gr->NY)); /* Should have been checked. */
                   demand(xy == (gr->NX*(uint32_t)y + (uint32_t)x), "inconsistent {vix} entry");
                 }
             }
         }
      }
  }

bool_t psgr_equal(psgr_t *gr, psgr_t *hr)
  {
    bool_t ok = TRUE;
    
    auto void compare_uint32(uint32_t gval, uint32_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_i2(i2_t *gval, i2_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_edge_data(psgr_edge_data_t *gval, psgr_edge_data_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_vertex_data(psgr_vertex_data_t *gval, psgr_vertex_data_t *hval, char *name1, int32_t ix, char *name2);
    
    compare_uint32(gr->NV, hr->NV, "NV", -1, NULL);
    compare_uint32(gr->NE, hr->NE, "NE", -1, NULL);
    
    compare_uint32(gr->NX, hr->NX, "NX", -1, NULL);
    compare_uint32(gr->NY, hr->NY, "NY", -1, NULL);
    
    for (int32_t v = 0; v < gr->NV; v++)
      { psgr_vertex_data_t *gvd = &(gr->vdata[v]);
        psgr_vertex_data_t *hvd = &(hr->vdata[v]);
        compare_vertex_data(gvd, hvd, "vertex", v, NULL);
      }
      
    for (int32_t e = 0; e < gr->NE; e++)
      { psgr_edge_data_t *ged = &(gr->edata[e]);
        psgr_edge_data_t *hed = &(hr->edata[e]);
        compare_edge_data(ged, hed, "edge", e, NULL);
      }
      
    uint32_t NXY = gr->NX*gr->NY;
    for (int32_t xy = 0; xy < NXY; xy++)
      { compare_uint32(gr->vix[xy].ix, hr->vix[xy].ix, "vix", xy, NULL); }
    
    return ok;
    
    void blast(char *name1, int32_t ix, char *name2)
      { fprintf(stderr, "** %s", name1);
        if (ix >= 0) { fprintf(stderr, " at index %d", ix); }
        if (name2 != NULL) { fprintf(stderr, " %s", name2); }
        fprintf(stderr, " differ:");
      }
    
    void compare_uint32(uint32_t gval, uint32_t hval, char *name1, int32_t ix, char *name2)
      { if (gval != hval)
          { blast(name1, ix, name2);
            fprintf(stderr, " %d %d\n", gval, hval);
            ok = FALSE;
          }
      }
    
    void compare_double(double gval, double hval, char *name1, int32_t ix, char *name2)
      { double tol = 2.0e-6; /* Considering the format of {psgr_write}. */
        if (fabs(gval - hval) > tol)
          { blast(name1, ix, name2);
            fprintf(stderr, " %24.16e %24.16e\n", gval, hval);
            ok = FALSE;
          }
      }
       
    void compare_i2(i2_t *gval, i2_t *hval, char *name1, int32_t ix, char *name2)
      { if ((gval->c[0] != hval->c[0]) || (gval->c[1] != hval->c[1]))
          { blast(name1, ix, name2);
            i2_print(stderr, gval);
            fprintf(stderr, " ");
            i2_print(stderr, hval);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
      
    void compare_edge_data(psgr_edge_data_t *gval, psgr_edge_data_t *hval, char *name1, int32_t ix, char *name2)
      { for (uint32_t db = 0; db <= 1; db++)
          { compare_uint32(gval->org[db].ix, hval->org[db].ix, name1, ix, "{org[db]}");
            compare_uint32(gval->enext[db].ix, hval->enext[db].ix, name1, ix, "{enext[db]}");
            compare_uint32(gval->eprev[db].ix, hval->eprev[db].ix, name1, ix, "{eprev[db]}");
          }
        compare_double(gval->delta, hval->delta, name1, ix, "{delta}");
        compare_double(gval->weight, hval->weight, name1, ix, "{weight}");
        /* !!! Should compare the paths. !!! */
      }
      
    void compare_vertex_data(psgr_vertex_data_t *gval, psgr_vertex_data_t *hval, char *name1, int32_t ix, char *name2)
      { 
        /* Don't care about the {mark}. */
        compare_i2(&(gval->pixel), &(hval->pixel), "vertex", ix, "{pixel}");
        compare_uint32(gval->aout.ix, hval->aout.ix, "arc out of vertex", ix, NULL);
      }
      
  }

psgr_arc_t psgr_find_enclosing_face(psgr_t *gr, r2_t *p)
  { /* Tries to locate the face that contains {p}. */
    /* !!! Should do this right !!! */
    psgr_vertex_t v = psgr_find_nearest_vertex(gr, p);
    psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
    assert(v.ix < gr->NV);
    psgr_arc_t a = vd->aout;
    return a;
  }

psgr_vertex_t psgr_find_nearest_vertex(psgr_t *gr, r2_t *p)
  { double min_dist_sqr = +INF;
    uint32_t min_v_ix = UINT32_MAX;
    for (uint32_t v_ix = 0; v_ix < gr->NV; v_ix++)
      { psgr_vertex_data_t *vd = &(gr->vdata[v_ix]);
        r2_t q = (r2_t){{ vd->pixel.c[0], vd->pixel.c[1] }};
        double d = r2_dist_sqr(p, &q);
        if (d < min_dist_sqr) { min_dist_sqr = d; min_v_ix = v_ix; }
      }
    demand(min_v_ix != UINT32_MAX, "graph has no valid vertices");
    return (psgr_vertex_t){ min_v_ix };
  }

double psgr_arc_left_face_curl(psgr_t *gr, psgr_arc_t a)
  { double curl = 0;
    psgr_arc_t b = a;
    do
      { curl += psgr_arc_delta(gr, b);
        b = psgr_arc_lnext(gr, b);
      } while (b.ix != a.ix);
    return curl;
  }

r2_t psgr_arc_left_face_barycenter(psgr_t *gr, psgr_arc_t a)
  {
    r2_t sum_p = (r2_t) {{0,0}};
    uint32_t count = 0;
    psgr_arc_t b = a;
    do 
      { psgr_vertex_t borg = psgr_arc_org(gr, b);
        assert(borg.ix < gr->NV);
        psgr_vertex_data_t *vd = &(gr->vdata[borg.ix]);
        r2_t p = (r2_t){{ vd->pixel.c[0], vd->pixel.c[1] }};
        r2_add(&sum_p, &p, &sum_p);
        count ++;
        b = psgr_arc_lnext(gr, b);
      } while(b.ix != a.ix);
      
    r2_t bar; r2_scale(1.0/(double)count, &sum_p, &bar);
    return bar;
  }

psgr_vertex_t psgr_vertex_from_pixel(psgr_t *gr, i2_t *ip)
  { demand((gr->NX > 0) && (gr->NY > 0), "graph has no associated map");
    int32_t x = ip->c[0], y = ip->c[1];
    demand((x >= 0) && (x < gr->NX), "invalid map index {x}");
    demand((y >= 0) && (y < gr->NY), "invalid map index {y}");
    uint32_t xy = gr->NX*(uint32_t)y + (uint32_t)x;
    return gr->vix[xy];
  }

void psgr_arc_print(FILE *wr, psgr_t *gr, psgr_arc_t a)
  { if (a.ix == ANONE.ix)
      { fprintf(wr, "-1"); }
    else
      { psgr_dir_bit_t db; psgr_edge_t e;
        psgr_arc_split(a, &e, &db);
        fprintf(wr, "%d:%d", e.ix, db);
      }
  }

void psgr_edge_remove(psgr_t *gr, psgr_edge_t e, int32_t indent)
  { demand(e.ix != ENONE.ix, "invalid null edge {e}");
    
    psgr_arc_t a0 = psgr_orient_edge(e, 0);
    psgr_arc_t a1 = psgr_orient_edge(e, 1);
    
    psgr_edge_data_t *ed = &(gr->edata[e.ix]);
    
    auto void disconnect_org(psgr_arc_t a);
      /* If the ring of arcs around the origin {org} of {a} is not
        just {a}, removes the arc {a} from that ring, by fixing the
        {onext} and {oprev} links or adjacent arcs.  Also sets 
        the {onext} and {oprev} links of {a} to {a} itself.  In any case,
        if {org.aout} is {a}, replaces it by some other arc 
        of that ring, or by {ANONE} if that was the only arc. */

    disconnect_org(a0);
    disconnect_org(a1);

    assert(psgr_arc_dprev(gr, a0).ix == a0.ix);
    assert(psgr_arc_dnext(gr, a0).ix == a0.ix);
    assert(psgr_arc_lnext(gr, a0).ix == a1.ix);
    assert(psgr_arc_lprev(gr, a0).ix == a1.ix);

    /* Clear out the edge data record {ed}: */
    ed->emark = DELETED;
    return;
    
    void disconnect_org(psgr_arc_t at)
      { assert(at.ix != ANONE.ix);
        
        psgr_edge_t et;
        psgr_dir_bit_t st;
        psgr_edge_data_t *edt;
        psgr_arc_unsafe_split(gr, at, &et, &st, &edt);
        assert(et.ix == e.ix);
        assert(edt == ed);
        psgr_arc_t onext_at = psgr_arc_onext(gr, at);
        psgr_arc_t oprev_at = psgr_arc_oprev(gr, at);
        psgr_vertex_t org_at = ed->org[st];
        gr->vdata[org_at.ix].aout = (onext_at.ix == at.ix ? ANONE : onext_at);
        
        if (at.ix != onext_at.ix) 
          { psgr_unsafe_arc_set_onext(gr, oprev_at, onext_at);
            psgr_unsafe_arc_set_oprev(gr, onext_at, oprev_at);
            psgr_unsafe_arc_set_onext(gr, at, at);
            psgr_unsafe_arc_set_oprev(gr, at, at);
          }
      }
  }

void psgr_vertex_remove(psgr_t *gr, psgr_vertex_t v, int32_t indent)
  { demand((v.ix != VNONE.ix) && (v.ix < gr->NV), "invalid vertex {v}");
    psgr_vertex_data_t *vd = &(gr->vdata[v.ix]);
    demand(vd->aout.ix == ANONE.ix, "vertex is still connected");
    demand(vd->vmark != DELETED, "vertex was already deleted");
    demand(vd->vmark != KEEP, "vertex should not be deleted");
    vd->vmark = DELETED;
  }
