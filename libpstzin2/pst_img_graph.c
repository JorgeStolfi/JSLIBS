/* See {pst_img_graph.h} */
/* Last edited on 2025-01-13 08:25:56 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <float_image.h>
#include <jsfile.h>
#include <filefmt.h>
#include <affirm.h>
#include <rn.h>

#include <pst_img_graph.h>

pst_img_graph_edge_data_t pst_img_graph_edge_data_create(void);
  /* Returns an edge data record with default information. */

void pst_img_graph_check_consistency_edge(pst_img_graph_t *g, haf_edge_t e);
  /* Runs some consistency checks on the edge {e} (which must be non-null)
    of the graph {g}.  Bombs out with message if any test fails. */
   
void pst_img_graph_edge_free(pst_img_graph_t *g, haf_edge_t e);
  /* Reclaims the heap data associated with edge {e}, including
    its label and path, and the {had_edge_rec_t} pointed to by {e}.
    Does NOT free the edge data record {g.edata[haf_edge_id(e)]}
    since that record is not heal-allocated. */ 
 
#define NONE pst_img_graph_mark_NONE


/* IMPLEMENTATIONS */

pst_img_graph_t *pst_img_graph_new(uint32_t NV_max, uint32_t NE_max)
  { pst_img_graph_t *g = talloc(1, pst_img_graph_t);
    g->NE_max = NE_max;
    g->NV_max = NV_max;
    g->NE = 0;
    g->NV = 0;
    
    g->vdata = talloc(NV_max, pst_img_graph_vertex_data_t);
    g->edata = talloc(NE_max, pst_img_graph_edge_data_t);
    g->hedge = talloc(NE_max, haf_edge_t);
    
    /* Initialization */
    for (int32_t i = 0; i < NE_max; i++) 
      { g->hedge[i] = NULL; g->edata[i].emark = NONE; }
      
    for (int32_t i = 0; i < NV_max; i++)
      { pst_img_graph_vertex_data_t *vd = &(g->vdata[i]);
        vd->aout = NULL;
        vd->x = -1;
        vd->y = -1;
        vd->vmark = pst_img_graph_mark_DELETED;
      }

    return g;
  }

pst_img_graph_edge_data_t pst_img_graph_edge_data_create(void)
  { pst_img_graph_edge_data_t ed;
    ed.org[0] = ed.org[1] = UINT32_MAX;  /*Index to the vextex in the list*/
    ed.delta = 0; /*Derivative*/
    ed.weight = 0; /* weight*/
    ed.emark = 0; /*marking number*/
    ed.label = NULL;
    ed.path = pst_path_create_empty();
    return ed;
  }

uint32_t pst_img_graph_add_vertex
  ( pst_img_graph_t *g,
    int32_t x,
    int32_t y,
    haf_arc_t aout,
    r2_t coords
  )
  { assert(g->NV < g->NV_max);
    uint32_t ix = g->NV;
    pst_img_graph_vertex_data_t *vd = &(g->vdata[ix]);

    vd->x = x;
    vd->y = y;
    vd->vmark = pst_img_graph_mark_NONE;
    vd->aout = aout;
    vd->coords = coords;
    g->NV++;
    return ix;
  }

int32_t pst_img_graph_get_arc_id(haf_arc_t a)
  { if (a == NULL) { return  -1; }
    return (int32_t)haf_arc_id(a);
  }

int32_t pst_img_graph_get_edge_id(haf_arc_t a)
  { if (a == NULL) { return -1; }
    return (int32_t)haf_edge_id(haf_edge(a));
  }

uint32_t pst_img_graph_get_arc_origin(pst_img_graph_t *g, haf_arc_t a)
  { demand(a != NULL, "arc is null");
    uint32_t eid = (uint32_t)haf_edge_id(haf_edge(a));
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    return ed->org[lbit];
  }

void pst_img_graph_set_arc_origin(pst_img_graph_t *g, haf_arc_t a, uint32_t org)
  { demand(a != NULL, "arc is null");
    demand(org < g->NV, "invalid vertex index");
    uint32_t eid = (uint32_t)haf_edge_id(haf_edge(a));
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    ed->org[lbit] = org;
  }

double pst_img_graph_get_edge_weight(pst_img_graph_t *g, haf_arc_t a)
  { uint32_t eid = (uint32_t)haf_edge_id(haf_edge(a));
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    return ed->weight;
  }

void pst_img_graph_set_edge_weight(pst_img_graph_t *g, haf_arc_t a, double w)
  { haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e); 
    demand(eid < g->NE, "invalid edge ID");
    assert((w == 0) || (w > 1.0e-100));
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    ed->weight = w;
  }

double pst_img_graph_get_arc_delta(pst_img_graph_t *g, haf_arc_t a)
  { haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e);
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    return (lbit == 0 ? +1.0 : -1.0)*ed->delta;
  }

void pst_img_graph_set_arc_delta(pst_img_graph_t *g, haf_arc_t a, double delta)
  { haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e);
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    ed->delta = (lbit == 0 ? +1.0 : -1.0)*delta;
  }

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t *g, haf_arc_t a)
  { haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e);
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    return (lbit == 0 ? ed->path : pst_path_reverse(ed->path));
  }

void pst_img_graph_set_edge_path(pst_img_graph_t *g, haf_arc_t a, pst_path_t p)
  { haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e);
    uint32_t lbit = haf_dir_bit(a);
    demand(eid < g->NE, "invalid edge ID");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    ed->path = (lbit == 0 ? p : pst_path_reverse(p));
  }

pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t *g, uint32_t vi0, uint32_t vi1, uint32_t vi2)
  { demand(vi0 < g->NV, "invalid vertex {vi0}");
    demand(vi1 < g->NV, "invalid vertex {vi1}");
    demand(vi2 < g->NV, "invalid vertex {vi2}");
    
    r2_t pi = g->vdata[vi0].coords;
    r2_t pe = g->vdata[vi1].coords;
    r2_t pf = g->vdata[vi2].coords;
    r2_t bar =
      (r2_t)
        {{  (pi.c[0] + pe.c[0] + pf.c[0] )/3.0, 
            (pi.c[1] + pe.c[1] + pf.c[1] )/3.0
        }};
    return pst_path_create_single(bar);
  }

haf_arc_t pst_img_graph_edge_add
  ( pst_img_graph_t *g,
    uint32_t org,  
    uint32_t dst,
    double d,
    double w,
    char *label,
    pst_path_t path
  )
  { demand(org < g->NV, "invalid vertex {org}");
    demand(dst < g->NV, "invalid vertex {dst}");
    if (g->NE == g->NE_max)
      { g->NE_max = (uint32_t)ceil(g->NE_max*1.5);
        g->hedge = retalloc(g->hedge, g->NE_max, haf_edge_t);
        g->edata = retalloc(g->edata, g->NE_max, pst_img_graph_edge_data_t);
      }

    uint32_t eid = g->NE;
    haf_arc_t a = haf_make_stick((haf_edge_id_t)eid);
    haf_edge_t e = haf_edge(a);
    assert(haf_dir_bit(a) == 0);
    assert(haf_edge_id(e) == eid);
    g->hedge[eid] = e;
    g->edata[eid] = pst_img_graph_edge_data_create();
    g->NE++;
    
    pst_img_graph_set_arc_origin(g, a, org);
    pst_img_graph_set_arc_origin(g, haf_sym(a), dst);
    g->vdata[org].aout = a;
    g->vdata[dst].aout = haf_sym(a);

    g->edata[eid].emark = NONE;
    g->edata[eid].weight = w;
    g->edata[eid].delta = d;
    g->edata[eid].label = label;
    g->edata[eid].path = path;

    return a;
  }

uint32_t pst_img_graph_vertex_out_degree(pst_img_graph_t *g, uint32_t vi)
  { demand(vi < g->NV, "invalid vertex index {vi}");
    pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
    if (vd->aout == NULL) return 0;
    uint32_t count = 1;
    haf_arc_t a = vd->aout;
    for (a = haf_onext(vd->aout); a!= vd->aout; a = haf_onext(a)) { count++; }
    return count;
  }

haf_arc_t pst_img_graph_get_connecting_arc(pst_img_graph_t *g, uint32_t vi0, uint32_t vi1)
  { demand (vi0 < g->NV, "invalid vertex index {vi0}");
    demand (vi1 < g->NV, "invalid vertex index {vi1}");
    pst_img_graph_vertex_data_t *vd0 = &(g->vdata[vi0]);
    pst_img_graph_vertex_data_t *vd1 = &(g->vdata[vi1]);
    if ((vd0->aout == NULL) || (vd1->aout == NULL))
      { return NULL; }
    haf_arc_t a = vd0->aout;
    do
      { uint32_t dst = pst_img_graph_get_arc_origin(g, haf_sym(a));
        if (dst == vi1) { return a; }
        a = haf_onext(a);
      } while(a != vd0->aout);
    return NULL;
  }

haf_arc_t pst_img_graph_find_rightmost_arc(pst_img_graph_t *g, haf_arc_t a)
  { demand(a != NULL, "invalid null arc {a}");
    haf_edge_t e = haf_edge(a);
    uint32_t eid = (uint32_t)haf_edge_id(e);
    demand(eid < g->NE, "invalid arc {a}");
    uint32_t org = pst_img_graph_get_arc_origin(g, a);
    assert(org < g->NV);
    haf_arc_t a_max = NULL;
    double cx_max = -INF;
    haf_arc_t b = a;
    do 
      { uint32_t dst = pst_img_graph_get_arc_origin(g,haf_sym(b));
        pst_img_graph_vertex_data_t *vd_dst = &(g->vdata[dst]);
        double cx = vd_dst->coords.c[0];
        if (cx > cx_max) { a_max = b; cx_max = cx; }
        b = haf_onext(b);
      } while(b != a);
    return a_max;
  }

void pst_img_graph_check_consistency(pst_img_graph_t *g)
  {
    demand(g->NV <= g->NV_max, "NV > NV_MAX"); 
    demand(g->NE <= g->NE_max, "NE > NE_MAX"); 
    for (uint32_t vi = 0; vi < g->NV; vi++)
      { pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
        demand((vd->x == -1) == (vd->y == -1), "inconsistent grid indices");
        haf_edge_t eout = haf_edge(vd->aout);
        pst_img_graph_check_consistency_edge(g, eout);
      }
    
    for (uint32_t eid = 0; eid < g->NE; eid++)
      { haf_edge_t e = g->hedge[eid];
        demand(e != NULL, "invalid null edge");
        demand(eid == (uint32_t)haf_edge_id(e), "haf edge {e} with inconsistent id");
        pst_img_graph_check_consistency_edge(g, e);
      }
  }

void pst_img_graph_check_consistency_edge(pst_img_graph_t *g, haf_edge_t e)  
  { demand(e != NULL, "invalid null edge {e}");
    uint32_t eid = (uint32_t)haf_edge_id(e);
    demand(eid < g->NE, "edge {e} has invalid id");
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    demand(isfinite(ed->weight) && (ed->weight >= 0), "invalid edge weight");
    demand(isfinite(ed->delta), "invalid edge delta");
    for (int32_t dbit = 0; dbit <= 1; dbit++)
      { uint32_t kv = ed->org[dbit];
        demand(kv < g->NV, "edge {e} has invalid origin or destination");
      }
    /* !!! should check {ed->path} !!! */
  }

bool_t pst_img_graph_equal(pst_img_graph_t *g, pst_img_graph_t *h)
  {
    bool_t ok = TRUE;
    
    auto void compare_int32(int32_t gval, int32_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_uint32(uint32_t gval, uint32_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_r2(r2_t *gval, r2_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_string(char *gval, char *hval, char *name1, int32_t ix, char *name2);
    auto void compare_r2(r2_t *gval, r2_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_arc(haf_arc_t gval, haf_arc_t hval, char *name1, int32_t ix, char *name2);
    auto void compare_edge_data(pst_img_graph_edge_data_t *gval, pst_img_graph_edge_data_t *hval, char *name1, int32_t ix, char *name2);
    auto void compare_vertex_data(pst_img_graph_vertex_data_t *gval, pst_img_graph_vertex_data_t *hval, char *name1, int32_t ix, char *name2);
    
    compare_uint32(g->NV, h->NV, "NV", -1, NULL);
    compare_uint32(g->NE, h->NE, "NE", -1, NULL);
    
    for (int32_t iv = 0; iv < g->NV; iv++)
      { pst_img_graph_vertex_data_t *gv = &(g->vdata[iv]);
        pst_img_graph_vertex_data_t *hv = &(h->vdata[iv]);
        compare_vertex_data(gv, hv, "vertex", iv, NULL);
      }
      
    for (int32_t ie = 0; ie < g->NE; ie++)
      { pst_img_graph_edge_data_t *ged = &(g->edata[ie]);
        pst_img_graph_edge_data_t *hed = &(h->edata[ie]);
        compare_edge_data(ged, hed, "edge", ie, NULL);
      }
    
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
      { double tol = 2.0e-6; /* Considering the format of {pst_img_graph_write}. */
        if (fabs(gval - hval) > tol)
          { blast(name1, ix, name2);
            fprintf(stderr, " %24.16e %24.16e\n", gval, hval);
            ok = FALSE;
          }
      }
    
    void compare_string(char *gval, char *hval, char *name1, int32_t ix, char *name2)
      { if (strcmp(gval, hval) != 0)
          { blast(name1, ix, name2);
            fprintf(stderr, " «%s» «%s»\n", gval, hval);
            ok = FALSE;
          }
      }
      
    void compare_r2(r2_t *gval, r2_t *hval, char *name1, int32_t ix, char *name2)
      { double d = r2_dist(gval, hval);
        double tol = 2.0e-6; /* Considering the format of {pst_img_graph_write}. */
        if (d > tol)
          { blast(name1, ix, name2);
            r2_print(stderr, gval);
            fprintf(stderr, " ");
            r2_print(stderr, hval);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
      
    void compare_arc(haf_arc_t gval, haf_arc_t hval, char *name1, int32_t ix, char *name2)
      { 
        haf_arc_id_t gaid = haf_arc_id(gval);
        haf_arc_id_t haid = haf_arc_id(hval);
        compare_uint32((uint32_t)gaid, (uint32_t)haid, name1, ix, "arc id");
      }
      
    void compare_edge_data(pst_img_graph_edge_data_t *gval, pst_img_graph_edge_data_t *hval, char *name1, int32_t ix, char *name2)
      { compare_uint32(gval->org[0], hval->org[0], name1, ix, "{org[0]}");
        compare_uint32(gval->org[1], hval->org[1], name1, ix, "{org[1]}");
        compare_double(gval->delta, hval->delta, name1, ix, "{delta}");
        compare_double(gval->weight, hval->weight, name1, ix, "{weight}");
        compare_string(gval->label, hval->label, name1, ix, "{label}");
        /* !!! Should check the paths. !!! */
      }
      
    void compare_vertex_data(pst_img_graph_vertex_data_t *gval, pst_img_graph_vertex_data_t *hval, char *name1, int32_t ix, char *name2)
      { 
        compare_int32(gval->x, hval->x, "vertex", ix, "{x}");
        compare_int32(gval->y, hval->y, "vertex", ix, "{y}");
        /* Don't care about the {mark}. */
        compare_r2(&(gval->coords), &(hval->coords), "vertex", ix, "{coords}");
        compare_arc(gval->aout, hval->aout, "edge of vertex", ix, NULL);
      }
        
  }

uint32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t *g, r2_t *p)
  { double min_dist_sqr = +INF;
    uint32_t min_vi = UINT32_MAX;
    for (uint32_t vi = 0; vi < g->NV; vi++)
      { pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
        double d = r2_dist_sqr(p, &(vd->coords));
        if (d < min_dist_sqr) { min_dist_sqr = d; min_vi = vi; }
      }
    demand(min_vi != UINT32_MAX, "graph has no valid vertices");
    return min_vi;
  }

double pst_img_graph_compute_left_face_curl(pst_img_graph_t *g, haf_arc_t a)
  { double curl = 0;
    haf_arc_t b = a;
    do
      { curl += pst_img_graph_get_arc_delta(g, b);
        b = haf_lnext(b);
      } while (b != a);
    return curl;
  }

r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t *g, haf_arc_t a)
  {
    r2_t sum_p = (r2_t) {{0,0}};
    uint32_t count = 0;
    haf_arc_t b = a;
    do 
      { uint32_t borg = pst_img_graph_get_arc_origin(g, b);
        assert(borg < g->NV);
        pst_img_graph_vertex_data_t *vd = &(g->vdata[borg]);
        r2_add(&sum_p, &(vd->coords), &sum_p);
        count ++;
        b = haf_lnext(b);
      } while(b != a);
      
    r2_t bar; r2_scale(1.0/(double)count, &sum_p, &bar);
    return bar;
  }

void pst_img_graph_put_curl_into_image(pst_img_graph_t *g, float_image_t *U)
  { int32_t NC_U, NX_U, NY_U;
    float_image_get_size(U, &NC_U, &NX_U, &NY_U);
    demand(NC_U == 1, "invalid channels in curl image");
    for (int32_t y = 0; y < NY_U; y++)
      { for (int32_t x = 0; x < NX_U; x++)
          { r2_t p = (r2_t){{ x + 0.5, y + 0.5 }};
            /* Tries to locate the face that contains {p}. */
            /* !!! Should do this right !!! */
            uint32_t vi = pst_img_graph_find_nearest_vertex(g, &p);
            pst_img_graph_vertex_data_t *v = &(g->vdata[vi]);
            assert(vi < g->NV);
            haf_arc_t a = v->aout;
            assert(a != NULL);
            double curl = pst_img_graph_compute_left_face_curl(g, a);
            float_image_set_sample(U, 0,x,y, (float)curl);
          }
      }
  }

void pst_img_graph_free(pst_img_graph_t *g)
  {
    /* Free the edge-data records and the edge list: */
    for (int32_t ke = 0; ke < g->NE; ke++)
      { haf_edge_t e = g->hedge[ke];
        if (e != NULL) pst_img_graph_edge_free(g, e);
      }
    free(g->hedge);
    free(g->edata);
    free(g->vdata);

    /* Free the graph: */
    free(g);
  }
     
void pst_img_graph_edge_free(pst_img_graph_t *g, haf_edge_t e)  
  { if (e == NULL) { return; }
    uint32_t eid = (uint32_t)haf_edge_id(e);
    assert(e == g->hedge[eid]);

    /* Trash contents of edge data record: */
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);
    ed->org[0] = UINT32_MAX;
    ed->org[1] = UINT32_MAX;
    ed->weight = 0;
    ed->delta = 0;
    free(ed->label); ed->label = NULL;
    pst_path_t *p = &(ed->path);
    if (p->v != NULL) { free(p->v); p->v = NULL; }

    haf_edge_free(e);
  }

