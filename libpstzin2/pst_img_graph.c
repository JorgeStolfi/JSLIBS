/* See {pst_img_graph.h} */
/* Last edited on 2025-01-05 12:22:10 by stolfi */
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
#include <affirm.h>
#include <rn.h>

#include <pst_img_graph.h>

void pst_img_graph_edge_data_free(pst_img_graph_edge_data_t* ed);

pst_img_graph_edge_data_t* pst_img_graph_edge_data_create(void);

void pst_img_graph_read_vertex(FILE* rd, pst_img_graph_t * g);

void pst_img_graph_read_edge(FILE* rd, pst_img_graph_t *g, int32_t list_onext[], int32_t list_long_bit[]);

void pst_img_graph_print_vertex(FILE* wr, pst_img_graph_vertex_data_t* v);

void pst_img_graph_print_edge(FILE* wr,pst_img_graph_t* g, haf_arc_t e);

void pst_img_graph_write_vertex(FILE* wr, pst_img_graph_vertex_data_t* v);

void pst_img_graph_write_edge(FILE* wr, pst_img_graph_t *g, haf_arc_t e);

#define MARK_VERTEX_NONE (-1)
  /* A null value for the {mark} field of a {pst_img_graph_edge_data_t} record. */

/* IMPLEMENTATIONS */

pst_img_graph_t *pst_img_graph_new(int32_t NV_max, int32_t NE_max)
  { pst_img_graph_t *g = talloc(1, pst_img_graph_t);
    g->NE_max = NE_max;
    g->NV_max = NV_max;
    g->NE = 0;
    g->NV = 0;
    g->NV_valid = 0;
    g->NE_valid = 0;
    
    g->vertex = talloc(NV_max, pst_img_graph_vertex_data_t);
    g->edge = talloc(NE_max, pst_img_graph_edge_t);
    
    /* Initialization */
    for (int32_t i = 0; i < NE_max; i++) 
      { g->edge[i].edge = NULL; g->edge[i].data = NULL; }
      
    for (int32_t i = 0; i < NV_max; i++)
      { g->vertex[i].edge = NULL;
        g->vertex[i].id = -1;
        g->vertex[i].mark = MARK_VERTEX_NONE;
      }

    return g;
  }

void pst_img_graph_edge_data_free(pst_img_graph_edge_data_t *ed)
  { free(ed); }

pst_img_graph_edge_data_t *pst_img_graph_edge_data_create(void)
  { pst_img_graph_edge_data_t *ed = talloc(1, pst_img_graph_edge_data_t);
    ed->id = 0;  /*Id - unique*/
    ed->org[0] = ed->org[1] = -1;  /*Index to the vextex in the list*/
    ed->delta[0] = ed->delta[1] = 0; /*Derivative*/
    ed->weight = 0; /* weight*/
    ed->mark = 0; /*marking number*/
    ed->label = NULL;
    ed->path[0] = ed->path[1] = pst_path_create_empty();
    return ed;
  }

int32_t  pst_img_graph_add_vertex
  ( pst_img_graph_t *g,
    int32_t id,
    haf_arc_t edge,
    r2_t coords
  )
  {  assert(g->NV < g->NV_max);
    int32_t ix = g->NV;
    pst_img_graph_vertex_data_t *v = &(g->vertex[ix]);

    v->id = id;
    v->mark = MARK_VERTEX_NONE;
    v->edge = edge;
    v->coords = coords;
    g->NV++;
    g->NV_valid++;
    return ix;
  }

int32_t pst_img_graph_get_dir_edge_num(haf_arc_t e)
  { if (e == NULL) { return  -1; }
    return (int32_t)haf_arc_id(e);
  }

int32_t pst_img_graph_get_edge_num(haf_arc_t e)
  { if (e == NULL) { return -1; }
    return (int32_t)haf_edge_id(e);
  }

int32_t pst_img_graph_get_edge_origin(pst_img_graph_t *g, haf_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->NE_max)));
    return g->edge[ind].data->org[lbit];
  }

void pst_img_graph_set_edge_origin(pst_img_graph_t *g, haf_arc_t e,int32_t org)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind < (2*g->NE_max)));
    g->edge[ind].data->org[lbit] = org;
  }

double pst_img_graph_get_edge_weight(pst_img_graph_t *g, haf_arc_t e)
  { int32_t ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->NE_max));
    return g->edge[ind].data->weight;
  }

void pst_img_graph_set_edge_weight(pst_img_graph_t *g, haf_arc_t e,double w)
  { int32_t ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->NE_max));
    assert((w == 0) || (w > 1.0e-100));
    g->edge[ind].data->weight = w;
  }

double pst_img_graph_get_edge_delta(pst_img_graph_t *g, haf_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->NE_max)));
    return g->edge[ind].data->delta[lbit];
  }

void pst_img_graph_set_edge_delta(pst_img_graph_t *g,haf_arc_t e,double delta)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->NE_max)));
    g->edge[ind].data->delta[lbit] = delta;
    g->edge[ind].data->delta[!lbit] = -delta;
  }

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t *g, haf_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->NE_max)));
    return g->edge[ind].data->path[lbit];
  }

void pst_img_graph_set_edge_path(pst_img_graph_t *g,haf_arc_t e,pst_path_t p)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert( (ind_dir >= 0) && (ind_dir < (2*g->NE_max)));
    g->edge[ind].data->path[lbit] = p;
    g->edge[ind].data->path[!lbit] = pst_path_reverse(p);
  }

pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t *g,int32_t vi, int32_t ve, int32_t vf)
  { r2_t pi = g->vertex[vi].coords;
    r2_t pe= g->vertex[ve].coords;
    r2_t pf = g->vertex[vf].coords;
    r2_t NE =
      (r2_t)
        {{  (pi.c[0] + pe.c[0] + pf.c[0] )/3.0, 
            (pi.c[1] + pe.c[1] + pf.c[1] )/3.0
        }};
    return pst_path_create_single(NE);
  }

haf_arc_t pst_img_graph_add_edge
  ( pst_img_graph_t *g,
    int32_t org,  
    int32_t dst,
    double d,
    double w,
    char *label,
    pst_path_t path
  )
  {
    if (g->NE == g->NE_max)
      { g->NE_max = (int32_t)ceil(g->NE_max*1.5);
        g->edge = retalloc(g->edge, g->NE_max, pst_img_graph_edge_t);
      }

    haf_arc_t e = haf_make_stick((haf_edge_id_t)g->NE);
    int32_t ind = pst_img_graph_get_edge_num(e);
    g->edge[ind].edge = e;
    g->edge[ind].data = pst_img_graph_edge_data_create();
    pst_img_graph_set_edge_origin(g,e,org);
    pst_img_graph_set_edge_origin(g,haf_sym(e),dst);
    g->vertex[org].edge = e;
    g->vertex[dst].edge = haf_sym(e);

    pst_img_graph_set_edge_weight(g,e,w);
    pst_img_graph_set_edge_delta(g,e,d);

    g->edge[ind].data->label = label;
    pst_img_graph_set_edge_path(g,e,path);
    g->NE++;
    g->NE_valid++;
    return e;
  }

void pst_img_graph_edge_remove(pst_img_graph_t *g,haf_arc_t e)
  { int32_t org_e = pst_img_graph_get_edge_origin(g,e);
    int32_t dst_e = pst_img_graph_get_edge_origin(g,haf_sym(e));

    haf_arc_t a = haf_oprev(e);
    haf_arc_t b = haf_oprev(haf_sym(e));

    pst_img_graph_set_edge_origin(g,e,-1);
    pst_img_graph_set_edge_origin(g,haf_sym(e),-1);
    pst_img_graph_set_edge_weight(g,e,0);
    pst_img_graph_set_edge_delta(g,e,0);

    pst_path_t p = pst_img_graph_get_edge_path(g,e);
    if (p.v != NULL) { free(p.v); }

    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t ind_e0 = pst_img_graph_get_dir_edge_num(e);
    int32_t ind_e1 = pst_img_graph_get_dir_edge_num(haf_sym(e));
    if (e != a) haf_splice(e,a);
    if (haf_sym(e) != b) { haf_splice(haf_sym(e),b); }

    bool_t test_dprev = ( pst_img_graph_get_dir_edge_num(haf_dprev(e)) == ind_e0 );
    bool_t test_dnext = ( pst_img_graph_get_dir_edge_num(haf_dnext(e)) == ind_e0 );
    bool_t test_lnext = ( pst_img_graph_get_dir_edge_num(haf_lnext(e)) == ind_e1 );
    bool_t test_lprev = ( pst_img_graph_get_dir_edge_num(haf_lprev(e)) == ind_e1 );

    if (!(test_dprev && test_dnext && test_lnext && test_lprev))
      { fprintf
          ( stderr,
            "REMOVAL FAILED AT [%d] E0 %d E1 %d\ndp %d dn %d  lp %d ln %d\nop %d on %d  rp %d rn %d\n",
            ind, ind_e0,ind_e1,
            pst_img_graph_get_dir_edge_num(haf_dprev(e)),
            pst_img_graph_get_dir_edge_num(haf_dnext(e)),
            pst_img_graph_get_dir_edge_num(haf_lprev(e)),
            pst_img_graph_get_dir_edge_num(haf_lnext(e)),
            pst_img_graph_get_dir_edge_num(haf_oprev(e)),
            pst_img_graph_get_dir_edge_num(haf_onext(e)),
            pst_img_graph_get_dir_edge_num(haf_rprev(e)),
            pst_img_graph_get_dir_edge_num(haf_rnext(e))
          );
      }
    assert(test_dprev && test_dnext && test_lnext && test_lprev);
    g->vertex[org_e].edge = (a == e ? NULL: a);
    g->vertex[dst_e].edge = (b == haf_sym(e) ? NULL: b);
    haf_free_edge(e);
    g->edge[ind].edge = NULL;
    pst_img_graph_edge_data_free(g->edge[ind].data);
    g->edge[ind].data = NULL;
    g->NE_valid--;
  }


void pst_img_graph_print_vertex(FILE *wr, pst_img_graph_vertex_data_t *v)
  { int32_t ind_edge_dir = pst_img_graph_get_dir_edge_num(v->edge);
    int32_t ind_edge = pst_img_graph_get_edge_num(v->edge);
    fprintf(wr,"ID: %d (%f,%f)  EDGE: %d (%d) \n",v->id,v->coords.c[0],v->coords.c[1],ind_edge_dir, ind_edge);
  }

void pst_img_graph_print_edge(FILE *wr,pst_img_graph_t *g, haf_arc_t e)
  { if (e == NULL) 
      { fprintf(wr," EMPTY "); }
    else
      {  int32_t ind = pst_img_graph_get_edge_num(e);
        double   w   = pst_img_graph_get_edge_weight(g,e);
        double   d   = pst_img_graph_get_edge_delta(g,e);
        double   rd  = pst_img_graph_get_edge_delta(g,haf_sym(e));
        int32_t org = pst_img_graph_get_edge_origin(g,e);
        int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
        pst_path_t p = pst_img_graph_get_edge_path(g,e);
        char *label = (g->edge[ind].data->label == NULL ? "": g->edge[ind].data->label);
        fprintf(wr,"ID: %d ORG: %d DEST: %d Delta: %lf RDelta %lf W: %lf LABEL %s PATH: %d ",ind,org,dst,d,rd,w,label,p.n);
        /*Neighbours*/
        fprintf(wr,"ONEXTS ");
        for (haf_arc_t oe = haf_onext(e); oe != e; oe = haf_onext(oe))
          { fprintf(wr, "%d:%d ",pst_img_graph_get_edge_num(oe),pst_img_graph_get_dir_edge_num(oe)&1); }
      }
    fprintf(wr,"\n");
  }

void pst_img_graph_print(FILE *wr,pst_img_graph_t *g)
  {
    for (int32_t i = 0; i < g->NV; i++)
      { fprintf(wr,"[%d] ",i);
        pst_img_graph_print_vertex(wr,&(g->vertex[i]));
      }
    fprintf(wr,"\n");
    for (int32_t i = 0; i < g->NE; i++)
      { fprintf(wr,"[%d] ",i);
        pst_img_graph_print_edge(wr,g,g->edge[i].edge);
      }
    fprintf(wr,"\n");
  }

int32_t pst_img_graph_vertex_count_neighbours(pst_img_graph_t *g, int32_t vi)
  { pst_img_graph_vertex_data_t *v = &(g->vertex[vi]);
    if (v->edge == NULL) return 0;
    int32_t count = 1;
    haf_arc_t e = v->edge;
    for (e = haf_onext(v->edge); e!= v->edge; e = haf_onext(e)) { count++; }
    return count;
  }

haf_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t *g,int32_t vi0, int32_t vi1)
  { pst_img_graph_vertex_data_t *v0 = &(g->vertex[vi0]);
    if ((v0->edge == NULL) || (v0->edge == NULL))
      { return NULL; }
    haf_arc_t e = v0->edge;
    do
      { int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
        if (dst == vi1) { return e; }
        e = haf_onext(e);
      } while(e != v0->edge);
    return NULL;
  }


haf_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t *g,haf_arc_t e0)
  {
    haf_arc_t e = e0;
    assert(e0 != NULL);
    int32_t org = pst_img_graph_get_edge_origin(g,e0);
    pst_img_graph_vertex_data_t *v_org = &(g->vertex[org]);
    haf_arc_t le = NULL;
    double cmax = -INF;
    do 
      { int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
        pst_img_graph_vertex_data_t *v_dst = &(g->vertex[dst]);
        r2_t d;
        r2_sub(&(v_dst->coords),&(v_org->coords),&d);
        (void) r2_dir(&d,&d);
        double c = d.c[0];
        if (c > cmax) { le = e; cmax = c; }
        e = haf_onext(e);
      } while(e != e0);
    return le;
  }
 
void debug_vertex_remove_edge(pst_img_graph_t *g, haf_arc_t e)
  {
    int32_t ind = pst_img_graph_get_edge_num(e);
    double   w   = pst_img_graph_get_edge_weight(g,e);
    double   d   = pst_img_graph_get_edge_delta(g,e);
    int32_t org = pst_img_graph_get_edge_origin(g,e);
    int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
    fprintf(stderr,"ID: %5d ORG: %5d DST: %5d D: %12.6lf W: %12.6lf",ind,org,dst,d,w);
  }

void pst_img_graph_check_consistency(pst_img_graph_t *g)
  {
    bool_t test = TRUE;
    for (int32_t i = 0; i < g->NV; i++)
      { pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (v->edge != NULL)
          { if (v->id == -1)
              { fprintf(stderr,"Dead-alive vertex [%d](%d) with edge not NULL\n",i,v->id);
                test = FALSE;
              }
            else
              { int32_t ind_edge = pst_img_graph_get_edge_num(v->edge);
                if (g->edge[ind_edge].edge == NULL)
                  { fprintf(stderr,"False link at vertex [%d](%d) pointing to [%d]\n",i,v->id,ind_edge);
                    test = FALSE;
                  }
                else
                  { int32_t orig_e = pst_img_graph_get_edge_origin(g,v->edge);
                    int32_t orig_ed = pst_img_graph_get_edge_origin(g,g->edge[ind_edge].edge);
                    int32_t orig_ed_s = pst_img_graph_get_edge_origin(g,haf_sym(g->edge[ind_edge].edge));
                    if ((orig_e != orig_ed) && (orig_e != orig_ed_s))
                      { fprintf(stderr,"Inconsitency of org-dest at vertex [%d](%d) pointing to [%d]\n",i,v->id,ind_edge);
                        test = FALSE;
                      }
                  }
              }
          }
    }
    if ( test)  fprintf(stderr,"Result - OK\n");
    else fprintf(stderr,"Result - FAIL\n");
  }
bool_t pst_img_graph_compare(pst_img_graph_t* g, pst_img_graph_t* h)
  {
    bool_t ok = TRUE;
    
    auto void check_uint32(int32_t ga, int32_t ha, char *name1, int32_t ix, char *name2);
    auto void check_r2(r2_t *ga, r2_t *ha, char *name1, int32_t ix, char *name2);
    auto void check_string(char *ga, char *ha, char *name1, int32_t ix, char *name2);
    auto void check_r2(r2_t *ga, r2_t *ha, char *name1, int32_t ix, char *name2);
    auto void check_edge_data(pst_img_graph_edge_data_t *ga, pst_img_graph_edge_data_t *ha, char *name1, int32_t ix, char *name2);
    
    
    check_uint32(g->NV, h->NV, "NV", -1, NULL);
    check_uint32(g->NE, h->NE, "NE", -1, NULL);
    
    for (int32_t iv = 0; iv < g->NV; iv++)
      { pst_img_graph_vertex_data_t *gv = &(g->vertex[iv]);
        pst_img_graph_vertex_data_t *hv = &(h->vertex[iv]);
        check_uint32(gv->id, hv->id, "vertex", iv, "{id}");
        /* Don't care about the {mark}. */
        check_r2(&(gv->coords), &(hv->coords), "vertex", iv, "{coords}");
        int32_t geid = pst_img_graph_get_dir_edge_num(gv->edge);
        int32_t heid = pst_img_graph_get_dir_edge_num(hv->edge);
        pst_img_graph_edge_data_t *ged = g->edge[geid].data;
        pst_img_graph_edge_data_t *hed = h->edge[heid].data;
        check_edge_data(ged, hed, "edge of vertex", iv, NULL);
      }
      
    for (int32_t ie = 0; ie < g->NE; ie++)
      { 
        pst_img_graph_edge_data_t *ged = g->edge[ie].data;
        pst_img_graph_edge_data_t *hed = h->edge[ie].data;
        check_edge_data(ged, hed, "edge", ie, NULL);
      }
    
    return ok;
    
    auto void blast(char *name1, int32_t ix, char *name2)
      { fprintf(stderr, "** %s", name1);
        if (ix >= 0) { fprintf(stderr, " at index %d", ix); }
        if (name2 != NULL) { fprintf(stderr, " %s", name2); }
        fprintf(stderr, " differ:");
      }
    
    auto void check_uint32(int32_t ga, int32_t ha, char *name1, int32_t ix, char *name2)
      { if (ga != ha)
          { blast(name1, ix, name2);
            fprintf(stderr, " %d %d\n", ga, ha);
            ok = FALSE;
          }
      }
    
    auto void check_double(double ga, double ha, char *name1, int32_t ix, char *name2)
      { double tol = 2.0e-6; /* Considering the format of {pst_img_graph_write}. */
        if (fabs(ga - ha) > tol)
          { blast(name1, ix, name2);
            fprintf(stderr, " %24.16e %24.16e\n", ga, ha);
            ok = FALSE;
          }
      }
    
    auto void check_string(char *ga, char *ha, char *name1, int32_t ix, char *name2)
      { if (strcmp(ga, ha) != 0)
          { blast(name1, ix, name2);
            fprintf(stderr, " «%s» «%s»\n", ga, ha);
            ok = FALSE;
          }
      }
      
    auto void check_r2(r2_t *ga, r2_t *ha, char *name1, int32_t ix, char *name2)
      { double d = r2_dist(ga, ha);
        double tol = 2.0e-6; /* Considering the format of {pst_img_graph_write}. */
        if (d > tol)
          { blast(name1, ix, name2);
            r2_print(stderr, ga);
            fprintf(stderr, " ");
            r2_print(stderr, ha);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
      
    auto void check_edge_data(pst_img_graph_edge_data_t *ga, pst_img_graph_edge_data_t *ha, char *name1, int32_t ix, char *name2)
      { check_uint32(ga->id, ha->id, name1, ix, "{id}");
        
        check_uint32(ga->org[0], ha->org[0], name1, ix, "{org[0]}");
        check_uint32(ga->dst[0], ha->dst[0], name1, ix, "{dst[0]}");
        check_double(ga->delta[0], ha->delta[0], name1, ix, "{delta[0]}");

        check_uint32(ga->org[1], ha->org[1], name1, ix, "{org[1]}");
        check_uint32(ga->dst[1], ha->dst[1], name1, ix, "{dst[1]}");
        check_double(ga->delta[1], ha->delta[1], name1, ix, "{delta[1]}");
        
        check_double(ga->weight, ha->weight, name1, ix, "{weight}");
        check_string(ga->label, ha->label, name1, ix, "{label}");
        /* !!! Should check the paths. !!! */
      }
        
  }
  
void pst_img_graph_reorganize(pst_img_graph_t *g)
  {
    demand(FALSE, "pst_img_graph_reorganize not implemented");
  }

pst_img_graph_t *pst_img_graph_copy(pst_img_graph_t *g)
  {
    /*First, count only valids*/

    auto haf_arc_t linha( haf_arc_t ed); 

    int32_t *eq_vector_vt = talloc(g->NV, int32_t);
    int32_t valid_NV = 0;
    for (int32_t i = 0; i < g->NV;i++)
      { if (g->vertex[i].id != -1)
          { eq_vector_vt[i] = valid_NV;
            valid_NV++;
          }
        else
          { eq_vector_vt[i] = -1; }
      }
    int32_t *eq_vector_ed = talloc(g->NE, int32_t);
    haf_arc_t *ref_tab = talloc(2*g->NE, haf_arc_t);
    int32_t valid_NE  = 0;
    for (int32_t i = 0; i < g->NE;i++)
      { if (g->edge[i].edge != NULL)
          { eq_vector_ed[i] = valid_NE;
            valid_NE++;
          }
        else
          { eq_vector_ed[i] = -1; }
      }
    pst_img_graph_t  *ng = pst_img_graph_new(valid_NV, valid_NE);
    /*Easy part, copy the raw data*/
    for (int32_t i = 0; i < g->NV;i++)
      { pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (eq_vector_vt[i] != -1)
          {  pst_img_graph_add_vertex(ng,v->id,NULL, v->coords ); }
      }
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->edge[i].edge;
        if (eq_vector_ed[i] != -1)
          { int32_t ind_e = pst_img_graph_get_edge_num(e);
            int32_t org =  pst_img_graph_get_edge_origin(g,e);
            int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
            pst_path_t p = pst_img_graph_get_edge_path(g,e);
            double delta = pst_img_graph_get_edge_delta(g,e);
            double weight = pst_img_graph_get_edge_weight(g,e);

            char *label = g->edge[ind_e].data->label;
            char *nl = (label == NULL ? NULL : talloc(strlen(label)+1, char));
            if (label != NULL ) strcpy(nl,label);

            pst_path_t np = p;
            if (p.v != NULL )
              { np.v = talloc(p.n, r2_t);
                memcpy(np.v,p.v,sizeof(r2_t)*(p.n));
              }
            /*We dont do much with it now...*/
            int32_t new_org = eq_vector_vt[org];
            int32_t new_dst = eq_vector_vt[dst];
            haf_arc_t new_e =  pst_img_graph_add_edge(ng,new_org,new_dst,delta,weight,nl,np);
            int32_t ind_e_dir = pst_img_graph_get_dir_edge_num(e);
            int32_t ind_e_sym = pst_img_graph_get_dir_edge_num(haf_sym(e));
            ref_tab[ind_e_dir] = new_e;
            ref_tab[ind_e_sym] = haf_sym(new_e);
          }
      }
      
    haf_arc_t linha( haf_arc_t ed)
      { int32_t ind =  pst_img_graph_get_dir_edge_num(ed);
        return ref_tab[ind];
      }

    /*Hard part, assemble the graph data*/
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->edge[i].edge;
        int32_t ind_el = eq_vector_ed[i];
        if (ind_el != -1)
          { haf_splice(haf_oprev(linha(e)),linha(haf_oprev(e)));
            haf_splice(haf_oprev(haf_sym(linha(e))),linha(haf_oprev(haf_sym(e))));
          }
      }

    free(eq_vector_vt);
    free(eq_vector_ed);
    free(ref_tab);
    return ng;
  }

void pst_img_graph_write_vertex(FILE *wr, pst_img_graph_vertex_data_t *v)
  {
    int32_t ind_edge = (v->edge == NULL ? -1 : pst_img_graph_get_dir_edge_num(v->edge));
    fprintf(wr,"%09d %d %9.6f %9.6f %09d\n",v->id,v->mark,v->coords.c[0],v->coords.c[1],ind_edge);
  }

void pst_img_graph_write_edge(FILE *wr, pst_img_graph_t *g, haf_arc_t e)
  {
    if (e == NULL)
      { fprintf(wr,"-1\n");
        return;
      }
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t org = pst_img_graph_get_edge_origin(g,e);
    int32_t dst = pst_img_graph_get_edge_origin(g,haf_sym(e));
    double delta = pst_img_graph_get_edge_delta(g,e);
    double weight = pst_img_graph_get_edge_weight(g,e);
    char *label = g->edge[ind].data->label;
    pst_path_t p = pst_img_graph_get_edge_path(g,e);

    fprintf(wr,"%09d %09d %09d %09d %9.6e %9.6e ",ind,ind_dir,org,dst,delta,weight);

    haf_arc_t e_next = haf_onext(e);
    int32_t ind_e_next = pst_img_graph_get_edge_num(e_next);
    int32_t ind_dir_e_next = pst_img_graph_get_dir_edge_num(e_next);
    int32_t long_bit_enext = ind_dir_e_next&1;

    haf_arc_t e_next_sym = haf_onext(haf_sym(e));

    int32_t ind_e_next_sym = pst_img_graph_get_edge_num(e_next_sym);
    int32_t ind_dir_e_next_sym = pst_img_graph_get_dir_edge_num(e_next_sym);
    int32_t long_bit_enext_sym = ind_dir_e_next_sym&1;

    fprintf(wr,"%d %d %d %d ",ind_e_next,long_bit_enext,ind_e_next_sym,long_bit_enext_sym);

    if ((label == NULL) || (strlen(label) == 0))
      { fprintf(wr,"0 "); }
    else
      { fprintf(wr, "%lu ", strlen(label));
        fprintf(wr, " %s", label);
      }
    pst_img_graph_write_path(wr,p);
    fprintf(wr,"\n");
  }

void pst_img_graph_write(FILE *wr, pst_img_graph_t *g)
  {
    /*First, write static data*/
    fprintf(wr,"%d %d %d %d \n",g->NV_valid, g->NE_valid,g->NV, g->NE);
    for (int32_t i = 0; i < g->NV; i++)
      {  pst_img_graph_write_vertex(wr,&(g->vertex[i]));  }
    for (int32_t i = 0; i < g->NE; i++)
      { pst_img_graph_write_edge(wr,g,g->edge[i].edge); }
    fprintf(wr,"\n");
  }

void pst_img_graph_read_vertex(FILE *wr, pst_img_graph_t  *g)
  {
    int32_t id;
    int32_t mark;
    r2_t coords;
    int32_t ind_edge;
    int32_t nsc = fscanf(wr,"%d %d %lf %lf %d",&id,&mark,&(coords.c[0]),&(coords.c[1]),&ind_edge);
    demand(nsc == 5, "Cannot read vertex");
    /*Who is responsible to fill the  real edge's value is the pst_img_read_edges*/
    int32_t ix = pst_img_graph_add_vertex(g,id,NULL,coords);
    g->vertex[ix].mark = mark;
  }

void pst_img_graph_read_edge
  ( FILE *wr,
    pst_img_graph_t *g,
    int32_t list_onext[],
    int32_t list_long_bit[]
  )
  { int32_t ind;
    int32_t ind_dir;
    int32_t org ;
    int32_t dst;
    double delta;
    double weight;

    demand(fscanf(wr,"%d",&ind) == 1, "Cannot read edge index");
    if (ind != -1)
      {
        int32_t ns1 = fscanf(wr,"%d %d %d %lf %lf",&ind_dir,&org,&dst,&delta,&weight);
        demand(ns1 == 5, "Cannot read edge data");
        int32_t ns2 = fscanf(wr,"%d %d %d %d",&(list_onext[0]),&(list_long_bit[0]),&(list_onext[1]),&(list_long_bit[1]));
        demand(ns2 == 4, "Cannot read edge connectivity");

        int32_t label_len;
        int32_t ns3 = fscanf(wr,"%d",&(label_len));
        demand(ns3 == 1, "Cannot read edge label lenght");
        demand(label_len >= 0, "invalid label length");
        char *label = (label_len == 0 ? NULL : talloc(label_len +1, char));
        if (label_len > 0)
          { fgetc(wr); /*skip the first blank space*/
            for (int32_t i = 0; i < label_len; i++)
              { fscanf(wr, "%c",&(label[i])); }
            label[label_len] = '\0';
          }
        pst_path_t p = pst_img_graph_read_path(wr);
        pst_img_graph_add_edge(g,org,dst,delta,weight,label,p);
      }
    else
      { g->edge[g->NE].edge = NULL;
        g->NE++;
      }
  }

pst_img_graph_t *pst_img_graph_read(FILE *wr)
  {
    int32_t NV_valid, NE_valid, NV, NE;
    demand( fscanf(wr,"%d %d %d %d", &NV_valid, &NE_valid, &NV, &NE) == 4, "parse of vertex and edge counts failed");
    demand((NV_valid >= 0) && (NE_valid >= 0) && (NV >= 0) && (NE >=0), "invalid vertex or edge counts");

    pst_img_graph_t *g = pst_img_graph_new(NV, NE);

    for (int32_t i = 0; i < NV ; i++) { pst_img_graph_read_vertex(wr,g); }
    int32_t *list_next = talloc(2*NE, int32_t);
    int32_t *list_long_bit = talloc(2*NE, int32_t);
    for (int32_t i = 0; i < NE ; i++)
      { pst_img_graph_read_edge(wr,g,&(list_next[2*i]),&(list_long_bit[2*i])); }

    auto haf_arc_t recover_edge(int32_t index, int32_t long_bit);

    haf_arc_t  recover_edge(int32_t index, int32_t long_bit)
      { assert((index >= 0 ) && (index < g->NE));
        haf_arc_t e = g->edge[index].edge;
        assert(e != NULL);
        if (long_bit == 1) { e = haf_sym(e); }
        return e;
      }

    /*Hard part*/
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->edge[i].edge;
        if (e != NULL)
          { haf_arc_t e_next = recover_edge(list_next[2*i],list_long_bit[2*i]);
            haf_arc_t e_next_sym = recover_edge(list_next[(2*i)+1],list_long_bit[(2*i)+1]);

            haf_splice( haf_oprev(e_next),e);
            haf_splice( haf_oprev(e_next_sym),haf_sym(e));
          }
      }
    free(list_next);
    free(list_long_bit);
    return g;
  }

int32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t *g, r2_t p)
  { double min_dist_sqr = +INF;
    int32_t min_i = -1;
    for (int32_t i = 0; i < g->NV; i++)
      {
        pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (v->id == -1) continue;
        double d= r2_dist_sqr(&p,&(v->coords));
        if (d < min_dist_sqr)
          { min_dist_sqr = d; min_i = i; }
      }

    return min_i;
  }

double pst_img_graph_compute_left_face_curl(pst_img_graph_t *g, haf_arc_t e0)
  { haf_arc_t e = e0;
    double curl = 0;
    do
      { curl+= pst_img_graph_get_edge_delta(g,e);
        e = haf_lnext(e);
      } while(e != e0);
    return curl;
  }

r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t *g, haf_arc_t e0)
  {
    haf_arc_t e = e0;
    r2_t p = (r2_t) {{0,0}};
    int32_t count= 0;
    do 
      { int32_t ind_org = pst_img_graph_get_edge_origin(g,e);
        assert(ind_org != -1);
        pst_img_graph_vertex_data_t *v = &(g->vertex[ind_org]);
        assert(v->id != -1);
        r2_add(&p,&(v->coords),&p);
        count++;
        e = haf_lnext(e);
      } while(e != e0);
    r2_scale(1.0/(double)count,&p,&p);
    return p;
  }

void pst_img_graph_put_curl_into_image(pst_img_graph_t *g, float_image_t *OZ)
  { int32_t NC, NX, NY;
    float_image_get_size(OZ, &NC, &NX, &NY);
    demand(NC == 1, "invalid channels in curl image");
    for (int32_t i = 0; i < g->NE; i++)
      { haf_arc_t e = g->edge[i].edge;
        if (e != NULL)
          { for (int32_t j = 0; j < 2; j++)
              { double curl = pst_img_graph_compute_left_face_curl(g,e);
                r2_t p = pst_img_graph_left_face_baricenter(g,e);
                int32_t x, y;
                pst_img_graph_get_vertex_image_indices(&p, NX, NY, &x, &y);
                float_image_set_sample(OZ, 0,x,y, (float)curl);
                e = haf_sym(e);
              }
          }
      }
  }

void pst_img_graph_free(pst_img_graph_t *g)
  {
    /* Free the edge-data records and the edge list: */
    for (int32_t i = 0; i < g->NE; i++)
      { haf_arc_t e = g->edge[i].edge;
        if (e != NULL) pst_img_graph_edge_remove(g,e);
        pst_img_graph_edge_data_t *edi = g->edge[i].data;
        if (edi != NULL)
          { pst_path_free(edi->path[0]);
            pst_path_free(edi->path[1]);
            pst_img_graph_edge_data_free(g->edge[i].data); 
          }
      }
    free(g->edge);
    
    /* Free the vertex list: */
    free(g->vertex);

    /* Free the graph: */
    free(g);
  }
