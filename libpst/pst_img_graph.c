/* See {pst_img_graph.h} */
/* Last edited on 2024-12-23 05:47:23 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <oct.h>
#include <float_image.h>
#include <jsfile.h>
#include <affirm.h>
#include <rn.h>

#include <pst_img_graph.h>

void pst_edge_data_free(pst_edge_data_t *ed)
  { free(ed); }

pst_edge_data_t *pst_edge_data_create(void)
  { pst_edge_data_t *ed = (pst_edge_data_t*)malloc(sizeof(pst_edge_data_t));
    ed->id = 0;  /*Id - unique*/
    ed->org[0] = ed->org[1] = -1;  /*Index to the vextex in the list*/
    ed->delta[0] = ed->delta[1] = 0; /*Derivative*/
    ed->weight = 0; /* weight*/
    ed->mark = 0; /*marking number*/
    ed->label = NULL;
    ed->path[0] = ed->path[1] = pst_path_create_empty();
    return ed;
  }

int32_t pst_img_graph_get_vertex_index_from_image_indices(int32_t ix, int32_t iy, int32_t NX, int32_t NY)
  { if ((iy < 0) || (iy > NY)) { return -1; }
    if ((ix < 0) || (ix > NX)) { return -1; }
    return iy*NX + ix;
  }

void pst_img_graph_get_vertex_image_indices(r2_t *p,int32_t NX, int32_t NY, int32_t *ix, int32_t *iy)
  { (*ix) = (int32_t)floor(p->c[0] + 0.5);
    (*iy) = (int32_t)floor(p->c[1] + 0.5);
    if ((*ix) >= NX ) { (*ix) = NX -1; }
    if ((*iy) >= NY ) { (*iy) = NY -1; }
    if ((*ix) < 0 ) { (*ix) = 0; }
    if ((*iy) < 0 ) { (*iy) = 0; }
  }

int32_t  pst_img_graph_vertex_add(pst_img_graph_t *g, int32_t id,oct_arc_t edge,r2_t coords)
  {  assert(g->n < g->max_n);
    int32_t ix = g->n;
    pst_vertex_data_t *v = &(g->vertex[ix]);

    v->id = id;
    v->mark = MARK_VERTEX_NONE;
    v->edge = edge;
    v->coords = coords;
    g->n++;
    g->n_valid++;
    return ix;
  }

int32_t pst_img_graph_get_dir_edge_num(oct_arc_t e)
  { if (e == oct_arc_NULL) { return  -1; }
    return (2*oct_edge_id(oct_edge(e))) + oct_lon_bit(e);
  }

int32_t pst_img_graph_get_edge_num(oct_arc_t e)
  { if (e == oct_arc_NULL) { return -1; }
    return oct_edge_id(oct_edge(e));
  }

int32_t pst_img_graph_get_edge_origin(pst_img_graph_t *g, oct_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->org[lbit];
  }

void pst_img_graph_set_edge_origin(pst_img_graph_t *g, oct_arc_t e,int32_t org)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind < (2*g->max_m)));
    g->edge[ind].data->org[lbit] = org;
  }

double pst_img_graph_get_edge_weight(pst_img_graph_t *g, oct_arc_t e)
  { int32_t ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->max_m));
    return g->edge[ind].data->weight;
  }

void pst_img_graph_set_edge_weight(pst_img_graph_t *g, oct_arc_t e,double w)
  { int32_t ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->max_m));
    assert((w == 0) || (w > 1.0e-100));
    g->edge[ind].data->weight = w;
  }

double pst_img_graph_get_edge_delta(pst_img_graph_t *g, oct_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->delta[lbit];
  }

void pst_img_graph_set_edge_delta(pst_img_graph_t *g,oct_arc_t e,double delta)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    g->edge[ind].data->delta[lbit] = delta;
    g->edge[ind].data->delta[!lbit] = -delta;
  }

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t *g, oct_arc_t e)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->path[lbit];
  }

void pst_img_graph_set_edge_path(pst_img_graph_t *g,oct_arc_t e,pst_path_t p)
  { int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t lbit = ind_dir&1;
    assert( (ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    g->edge[ind].data->path[lbit] = p;
    g->edge[ind].data->path[!lbit] = pst_path_reverse(p);
  }

pst_path_t pst_path_create_empty(void)
  { pst_path_t p = (pst_path_t){ .n=0, .v=NULL, .reverse=FALSE }; 
    return p;
  }

pst_path_t pst_path_create_single(r2_t coords)
  { r2_t *v = (r2_t*)malloc(sizeof(r2_t));
    v[0] = coords;
    pst_path_t p = (pst_path_t){ .n=1, .v=v, .reverse=FALSE }; 
    return p;
  }

pst_path_t pst_path_reverse(pst_path_t p)
  { p.reverse = !p.reverse;
    return p;
  }

r2_t pst_path_get_vertex(pst_path_t p, int32_t i)
  { return ( p.reverse ? p.v[p.n - i -1 ] : p.v[i] );  }

pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1)
  { int32_t n0 = p0.n;
    int32_t n1 = p1.n;
    int32_t n = n0 + n1 + 1;

    r2_t *v = (r2_t*)malloc(sizeof(r2_t)*n);
    int32_t i;
    for (uint32_t i = 0; i < n; i++)
      { r2_t c;
        if (i < n0)
          { c = pst_path_get_vertex(p0,i); }
        else if (i == n0)
          { c = coords; }
        else
          { int32_t j = i - n0 - 1;
            c = pst_path_get_vertex(p1,j);
          }
        v[i] = c;
      }
    return (pst_path_t) { .n=n, .v=v, .reverse=FALSE };
  }

pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t *g,int32_t vi, int32_t ve, int32_t vf)
  { r2_t pi = g->vertex[vi].coords;
    r2_t pe= g->vertex[ve].coords;
    r2_t pf = g->vertex[vf].coords;
    r2_t m =
      (r2_t)
        {{  (pi.c[0] + pe.c[0] + pf.c[0] )/3.0, 
            (pi.c[1] + pe.c[1] + pf.c[1] )/3.0
        }};
    return pst_path_create_single(m);
  }

oct_arc_t pst_img_graph_edge_insert
  ( pst_img_graph_t *g,
    int32_t org,int32_t dst,
    double d,
    double w,
    char *label,
    pst_path_t path
  )
  {
    if (g->m == g->max_m)
      { g->max_m = (g->max_m*1.5);
        g->edge = (pst_edge_t*)notnull(realloc(g->edge,sizeof(pst_edge_t)*(g->max_m)),"bug");
      }

    oct_arc_t e = oct_make_edge();
    oct_set_edge_num(oct_edge(e), g->m);
    int32_t ind = pst_img_graph_get_edge_num(e);
    g->edge[ind].edge = e;
    g->edge[ind].data = pst_edge_data_create();
    pst_img_graph_set_edge_origin(g,e,org);
    pst_img_graph_set_edge_origin(g,oct_sym(e),dst);
    g->vertex[org].edge = e;
    g->vertex[dst].edge = oct_sym(e);

    pst_img_graph_set_edge_weight(g,e,w);
    pst_img_graph_set_edge_delta(g,e,d);

    g->edge[ind].data->label = label;
    pst_img_graph_set_edge_path(g,e,path);
    g->m++;
    g->m_valid++;
    return e;
  }

void pst_img_graph_edge_remove(pst_img_graph_t *g,oct_arc_t e)
  { int32_t org_e = pst_img_graph_get_edge_origin(g,e);
    int32_t dst_e = pst_img_graph_get_edge_origin(g,oct_sym(e));

    oct_arc_t a = oct_oprev(e);
    oct_arc_t b = oct_oprev(oct_sym(e));

    pst_img_graph_set_edge_origin(g,e,-1);
    pst_img_graph_set_edge_origin(g,oct_sym(e),-1);
    pst_img_graph_set_edge_weight(g,e,0);
    pst_img_graph_set_edge_delta(g,e,0);

    pst_path_t p = pst_img_graph_get_edge_path(g,e);
    if (p.v != NULL) { free(p.v); }

    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t ind_e0 = pst_img_graph_get_dir_edge_num(e);
    int32_t ind_e1 = pst_img_graph_get_dir_edge_num(oct_sym(e));
    if (e != a) oct_splice(e,a);
    if (oct_sym(e) != b) { oct_splice(oct_sym(e),b); }

    bool_t test_dprev = ( pst_img_graph_get_dir_edge_num(oct_dprev(e)) == ind_e0 );
    bool_t test_dnext = ( pst_img_graph_get_dir_edge_num(oct_dnext(e)) == ind_e0 );
    bool_t test_lnext = ( pst_img_graph_get_dir_edge_num(oct_lnext(e)) == ind_e1 );
    bool_t test_lprev = ( pst_img_graph_get_dir_edge_num(oct_lprev(e)) == ind_e1 );

    if (!(test_dprev && test_dnext && test_lnext && test_lprev))
      { fprintf
          ( stderr,
            "REMOVAL FAILED AT [%ld] E0 %ld E1 %ld\ndp %ld dn %ld  lp %ld ln %ld\nop %ld on %ld  rp %ld rn %ld\n",
            ind, ind_e0,ind_e1,
            pst_img_graph_get_dir_edge_num(oct_dprev(e)),
            pst_img_graph_get_dir_edge_num(oct_dnext(e)),
            pst_img_graph_get_dir_edge_num(oct_lprev(e)),
            pst_img_graph_get_dir_edge_num(oct_lnext(e)),
            pst_img_graph_get_dir_edge_num(oct_oprev(e)),
            pst_img_graph_get_dir_edge_num(oct_onext(e)),
            pst_img_graph_get_dir_edge_num(oct_rprev(e)),
            pst_img_graph_get_dir_edge_num(oct_rnext(e))
          );
      }
    assert(test_dprev && test_dnext && test_lnext && test_lprev);
    g->vertex[org_e].edge = (a == e ? oct_arc_NULL: a);
    g->vertex[dst_e].edge = (b == oct_sym(e) ? oct_arc_NULL: b);
    oct_destroy_edge(e);
    g->edge[ind].edge = oct_arc_NULL;
    pst_edge_data_free(g->edge[ind].data);
    g->edge[ind].data = NULL;
    g->m_valid--;
  }

pst_img_graph_t *pst_img_graph_create(int32_t max_m,int32_t max_n)
  { pst_img_graph_t *g = (pst_img_graph_t*)malloc(sizeof(pst_img_graph_t));
    g->max_m = max_m;
    g->max_n = max_n;
    g->m = 0;
    g->n = 0;
    g->n_valid = 0;
    g->m_valid = 0;
    
    g->vertex = (pst_vertex_data_t*)malloc(sizeof(pst_vertex_data_t)*(g->max_n));
    g->edge = (pst_edge_t*)malloc(sizeof(pst_edge_t)*(g->max_m));
    
    /* Initialization */
    int32_t i;
    for (uint32_t i = 0; i < max_m; i++) { g->edge[i].edge = oct_arc_NULL; g->edge[i].data = NULL; }
    for (uint32_t i = 0; i < max_n; i++)
      { g->vertex[i].edge = oct_arc_NULL;
        g->vertex[i].id = -1;
        g->vertex[i].mark = MARK_VERTEX_NONE;
      }

    return g;
  }

void pst_img_graph_print_vertex(FILE *wr, pst_vertex_data_t *v)
  { int32_t ind_edge_dir = pst_img_graph_get_dir_edge_num(v->edge);
    int32_t ind_edge = pst_img_graph_get_edge_num(v->edge);
    fprintf(wr,"ID: %ld (%f,%f)  EDGE: %ld (%ld) \n",v->id,v->coords.c[0],v->coords.c[1],ind_edge_dir, ind_edge);
  }

void pst_img_graph_print_edge(FILE *wr,pst_img_graph_t *g, oct_arc_t e)
  { if (e == oct_arc_NULL) 
      { fprintf(wr," EMPTY "); }
    else
      {  int32_t ind = pst_img_graph_get_edge_num(e);
        double   w   = pst_img_graph_get_edge_weight(g,e);
        double   d   = pst_img_graph_get_edge_delta(g,e);
        double   rd  = pst_img_graph_get_edge_delta(g,oct_sym(e));
        int32_t org = pst_img_graph_get_edge_origin(g,e);
        int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
        pst_path_t p = pst_img_graph_get_edge_path(g,e);
        char *label = (g->edge[ind].data->label == NULL ? "": g->edge[ind].data->label);
        fprintf(wr,"ID: %ld ORG: %ld DEST: %ld Delta: %lf RDelta %lf W: %lf LABEL %s PATH: %ld ",ind,org,dst,d,rd,w,label,p.n);
        /*Neighbours*/
        fprintf(wr,"ONEXTS ");
        oct_arc_t oe;
        for (oct_arc_t oe = oct_onext(e); oe != e; oe = oct_onext(oe))
          { fprintf(wr, "%ld:%ld ",pst_img_graph_get_edge_num(oe),pst_img_graph_get_dir_edge_num(oe)&1); }
      }
    fprintf(wr,"\n");
  }

void pst_img_graph_print(FILE *wr,pst_img_graph_t *g)
  {
    int32_t i;
    for (uint32_t i = 0; i < g->n; i++)
      { fprintf(wr,"[%ld] ",i);
        pst_img_graph_print_vertex(wr,&(g->vertex[i]));
      }
    fprintf(wr,"\n");
    for (uint32_t i = 0; i < g->m; i++)
      { fprintf(wr,"[%ld] ",i);
        pst_img_graph_print_edge(wr,g,g->edge[i].edge);
      }
    fprintf(wr,"\n");
  }

int32_t pst_img_graph_vertex_count_neighbours(pst_img_graph_t *g, int32_t vi)
  { pst_vertex_data_t *v = &(g->vertex[vi]);
    if (v->edge == oct_arc_NULL) return 0;
    int32_t count = 1;
    oct_arc_t e = v->edge;
    for (e = oct_onext(v->edge); e!= v->edge; e = oct_onext(e)) { count++; }
    return count;
  }

oct_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t *g,int32_t vi0, int32_t vi1)
  { pst_vertex_data_t *v0 = &(g->vertex[vi0]);
    if ((v0->edge == oct_arc_NULL) || (v0->edge == oct_arc_NULL))
      { return oct_arc_NULL; }
    oct_arc_t e = v0->edge;
    do
      { int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
        if (dst == vi1) { return e; }
        e = oct_onext(e);
      } while(e != v0->edge);
    return oct_arc_NULL;
  }


oct_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t *g,oct_arc_t e0)
  {
    oct_arc_t e = e0;
    assert(e0 != oct_arc_NULL);
    int32_t org = pst_img_graph_get_edge_origin(g,e0);
    pst_vertex_data_t *v_org = &(g->vertex[org]);
    oct_arc_t le = oct_arc_NULL;
    double cmax = -INF;
    do 
      { int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
        pst_vertex_data_t *v_dst = &(g->vertex[dst]);
        r2_t d;
        r2_sub(&(v_dst->coords),&(v_org->coords),&d);
        (void) r2_dir(&d,&d);
        double c = d.c[0];
        if (c > cmax) { le = e; cmax = c; }
        e = oct_onext(e);
      } while(e != e0);
    return le;
  }


void debug_vertex_remove_edge(pst_img_graph_t *g, oct_arc_t e)
  {
    int32_t ind = pst_img_graph_get_edge_num(e);
    double   w   = pst_img_graph_get_edge_weight(g,e);
    double   d   = pst_img_graph_get_edge_delta(g,e);
    int32_t org = pst_img_graph_get_edge_origin(g,e);
    int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
    fprintf(stderr,"ID: %5ld ORG: %5ld DST: %5ld D: %12.6lf W: %12.6lf",ind,org,dst,d,w);
  }

void pst_img_graph_check_consistency(pst_img_graph_t *g)
  {
    int32_t i;
    bool_t test = TRUE;
    for (uint32_t i = 0; i < g->n; i++)
      { pst_vertex_data_t *v = &(g->vertex[i]);
        if (v->edge != oct_arc_NULL)
          { if (v->id == -1)
              { fprintf(stderr,"Dead-alive vertex [%ld](%ld) with edge not NULL\n",i,v->id);
                test = FALSE;
              }
            else
              { int32_t ind_edge = pst_img_graph_get_edge_num(v->edge);
                if (g->edge[ind_edge].edge == oct_arc_NULL)
                  { fprintf(stderr,"False link at vertex [%ld](%ld) pointing to [%ld]\n",i,v->id,ind_edge);
                    test = FALSE;
                  }
                else
                  { int32_t orig_e = pst_img_graph_get_edge_origin(g,v->edge);
                    int32_t orig_ed = pst_img_graph_get_edge_origin(g,g->edge[ind_edge].edge);
                    int32_t orig_ed_s = pst_img_graph_get_edge_origin(g,oct_sym(g->edge[ind_edge].edge));
                    if ((orig_e != orig_ed) && (orig_e != orig_ed_s))
                      { fprintf(stderr,"Inconsitency of org-dest at vertex [%ld](%ld) pointing to [%ld]\n",i,v->id,ind_edge);
                        test = FALSE;
                      }
                  }
              }
          }
    }
    if ( test)  fprintf(stderr,"Result - OK\n");
    else fprintf(stderr,"Result - FAIL\n");
  }

void pst_img_graph_reorganize(pst_img_graph_t *g)
  {
    demand(FALSE, "pst_img_graph_reorganize not implemented");
  }

pst_img_graph_t *pst_img_graph_copy(pst_img_graph_t *g)
  {
    /*First, count only valids*/

    auto oct_arc_t linha( oct_arc_t ed); 

    int32_t i;
    int32_t *eq_vector_vt = (int32_t*)malloc(sizeof(int32_t)*(g->n));
    int32_t valid_n = 0;
    for (uint32_t i = 0; i < g->n;i++)
      { if (g->vertex[i].id != -1)
          { eq_vector_vt[i] = valid_n;
            valid_n++;
          }
        else
          { eq_vector_vt[i] = -1; }
      }
    int32_t *eq_vector_ed = (int32_t*)malloc(sizeof(int32_t)*(g->m));
    oct_arc_t *ref_tab = (oct_arc_t*)malloc(sizeof(oct_arc_t)*(g->m)*2);
    int32_t valid_m  = 0;
    for (uint32_t i = 0; i < g->m;i++)
      { if (g->edge[i].edge != oct_arc_NULL)
          { eq_vector_ed[i] = valid_m;
            valid_m++;
          }
        else
          { eq_vector_ed[i] = -1; }
      }
    pst_img_graph_t  *ng = pst_img_graph_create(valid_m,valid_n);
    /*Easy part, copy the raw data*/
    for (uint32_t i = 0; i < g->n;i++)
      { pst_vertex_data_t *v = &(g->vertex[i]);
        if (eq_vector_vt[i] != -1)
          {  pst_img_graph_vertex_add(ng,v->id,oct_arc_NULL, v->coords ); }
      }
    for (uint32_t i = 0; i < g->m;i++)
      { oct_arc_t e = g->edge[i].edge;
        if (eq_vector_ed[i] != -1)
          { int32_t ind_e = pst_img_graph_get_edge_num(e);
            int32_t org =  pst_img_graph_get_edge_origin(g,e);
            int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
            pst_path_t p = pst_img_graph_get_edge_path(g,e);
            double delta = pst_img_graph_get_edge_delta(g,e);
            double weight = pst_img_graph_get_edge_weight(g,e);

            char *label = g->edge[ind_e].data->label;
            char *nl = (label == NULL ? NULL : (char*)malloc(sizeof(char)*(strlen(label)+1)));
            if (label != NULL ) strcpy(nl,label);

            pst_path_t np = p;
            if (p.v != NULL )
              { np.v = (r2_t*)malloc(sizeof(r2_t)*(p.n));
                memcpy(np.v,p.v,sizeof(r2_t)*(p.n));
              }
            /*We dont do much with it now...*/
            int32_t new_org = eq_vector_vt[org];
            int32_t new_dst = eq_vector_vt[dst];
            oct_arc_t new_e =  pst_img_graph_edge_insert(ng,new_org,new_dst,delta,weight,nl,np);
            int32_t ind_e_dir = pst_img_graph_get_dir_edge_num(e);
            int32_t ind_e_sym = pst_img_graph_get_dir_edge_num(oct_sym(e));
            ref_tab[ind_e_dir] = new_e;
            ref_tab[ind_e_sym] = oct_sym(new_e);
          }
      }
      
    oct_arc_t linha( oct_arc_t ed)
      { int32_t ind =  pst_img_graph_get_dir_edge_num(ed);
        return ref_tab[ind];
      }

    /*Hard part, assemble the graph data*/
    for (uint32_t i = 0; i < g->m;i++)
      { oct_arc_t e = g->edge[i].edge;
        int32_t ind_el = eq_vector_ed[i];
        if (ind_el != -1)
          { oct_splice(oct_oprev(linha(e)),linha(oct_oprev(e)));
            oct_splice(oct_oprev(oct_sym(linha(e))),linha(oct_oprev(oct_sym(e))));
          }
      }

    free(eq_vector_vt);
    free(eq_vector_ed);
    free(ref_tab);
    return ng;
  }

void pst_img_graph_write_vertex(FILE *wr, pst_vertex_data_t *v)
  {
    int32_t ind_edge = (v->edge == oct_arc_NULL ? -1 : pst_img_graph_get_dir_edge_num(v->edge));
    fprintf(wr,"%09ld %d %9.6f %9.6f %09ld\n",v->id,v->mark,v->coords.c[0],v->coords.c[1],ind_edge);
  }

void pst_img_graph_write_path(FILE *wr, pst_path_t p)
  {
    fprintf(wr,"%ld %d ",p.n,p.reverse);
    int32_t i;
    for (uint32_t i = 0; i < p.n; i++)
      { fprintf(wr,"%9.6f %9.6f ",p.v[i].c[0],p.v[i].c[1]); }
  }

void pst_img_graph_write_edge(FILE *wr, pst_img_graph_t *g, oct_arc_t e)
  {
    if (e == oct_arc_NULL)
      { fprintf(wr,"-1\n");
        return;
      }
    int32_t ind = pst_img_graph_get_edge_num(e);
    int32_t ind_dir = pst_img_graph_get_dir_edge_num(e);
    int32_t org = pst_img_graph_get_edge_origin(g,e);
    int32_t dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
    double delta = pst_img_graph_get_edge_delta(g,e);
    double weight = pst_img_graph_get_edge_weight(g,e);
    char *label = g->edge[ind].data->label;
    pst_path_t p = pst_img_graph_get_edge_path(g,e);

    fprintf(wr,"%09ld %09ld %09ld %09ld %9.6e %9.6e ",ind,ind_dir,org,dst,delta,weight);

    oct_arc_t e_next = oct_onext(e);
    int32_t ind_e_next = pst_img_graph_get_edge_num(e_next);
    int32_t ind_dir_e_next = pst_img_graph_get_dir_edge_num(e_next);
    int32_t long_bit_enext = ind_dir_e_next&1;

    oct_arc_t e_next_sym = oct_onext(oct_sym(e));

    int32_t ind_e_next_sym = pst_img_graph_get_edge_num(e_next_sym);
    int32_t ind_dir_e_next_sym = pst_img_graph_get_dir_edge_num(e_next_sym);
    int32_t long_bit_enext_sym = ind_dir_e_next_sym&1;

    fprintf(wr,"%ld %d %ld %d ",ind_e_next,long_bit_enext,ind_e_next_sym,long_bit_enext_sym);

    if (label == NULL)
      { fprintf(wr,"0 "); }
    else
      { fprintf(wr,"%d ",strlen(label));
        fprintf(wr,"%s ",label);
      }
    pst_img_graph_write_path(wr,p);
    fprintf(wr,"\n");
  }

void pst_img_graph_write(FILE *wr, pst_img_graph_t *g)
  {
    /*First, write static data*/
    fprintf(wr,"%ld %ld %ld %ld \n",g->n_valid, g->m_valid,g->n, g->m);
    int32_t i;
    for (uint32_t i = 0; i < g->n; i++)
      {  pst_img_graph_write_vertex(wr,&(g->vertex[i]));  }
    for (uint32_t i = 0; i < g->m; i++)
      { pst_img_graph_write_edge(wr,g,g->edge[i].edge); }
    fprintf(wr,"\n");
  }

void pst_img_graph_read_vertex(FILE *wr, pst_img_graph_t  *g)
  {
    int32_t id;
    int32_t mark;
    r2_t coords;
    int32_t ind_edge;
    int32_t nsc = fscanf(wr,"%ld %d %lf %lf %ld",&id,&mark,&(coords.c[0]),&(coords.c[1]),&ind_edge)
    demand(nsc == 5, "Cannot read vertex");
    /*Who is responsible to fill the  real edge's value is the pst_img_read_edges*/
    int32_t ix = pst_img_graph_vertex_add(g,id,oct_arc_NULL,coords);
    g->vertex[ix].mark = mark;
  }

pst_path_t pst_img_graph_read_path(FILE *wr)
  {
    pst_path_t p;
    int32_t reverse;
    demand(fscanf(wr,"%ld %d",&(p.n),&(reverse)) == 2, "Cannot read path");
    p.reverse = (reverse == 1);
    p.v = NULL;
    if (p.n == 0) return p;
    p.v = (r2_t*) malloc(sizeof(r2_t)*(p.n));
    int32_t i;
    for (uint32_t i = 0; i < p.n; i++)
      { demand(fscanf(wr,"%lf %lf",&(p.v[i].c[0]),&(p.v[i].c[1])) == 2, "Cannot read path elements"); }
    return p;
  }

void pst_img_graph_read_edge(FILE *wr, pst_img_graph_t *g, int32_t list_onext[], int32_t list_long_bit[])
  { int32_t ind;
    int32_t ind_dir;
    int32_t org ;
    int32_t dst;
    double delta;
    double weight;

    demand(fscanf(wr,"%ld",&ind) == 1, "Cannot read edge index");
    if (ind != -1)
      {
        int32_t ns1 = fscanf(wr,"%ld %ld %ld %lf %lf",&ind_dir,&org,&dst,&delta,&weight);
        demand(ns == 5, "Cannot read edge data");
        int32_t ns2 = fscanf(wr,"%ld %d %ld %d",&(list_onext[0]),&(list_long_bit[0]),&(list_onext[1]),&(list_long_bit[1]));
        demand(ns2 == 4, "Cannot read edge connectivity");

        int32_t label_len;
        int32_ t ns3 = fscanf(wr,"%ld",&(label_len));
        demand(ns3 == 1, "Cannot read edge label lenght");
        char *label = (label_len == 0  ? NULL : (char*)malloc(sizeof(char)*(label_len +1)));
        int32_t i;
        fgetc(wr); /*skip the first blank space*/
        for (uint32_t i = 0; i < label_len; i++)
          { fscanf(wr,"%c",&(label[i])); }
        if (label != NULL) label[i] = '\0';
        pst_path_t p = pst_img_graph_read_path(wr);
        pst_img_graph_edge_insert(g,org,dst,delta,weight,label,p);
      }
    else
      { g->edge[g->m].edge = oct_arc_NULL;
        g->m++;
      }
  }

pst_img_graph_t *pst_img_graph_read(FILE *wr)
  {
    int32_t n_valid,m_valid, n,m;
    demand( fscanf(wr,"%ld %ld %ld %ld",&n_valid,&m_valid,&n,&m) == 4,"Cannot read graph header properly");
    assert( (n_valid >= 0) && (m_valid >= 0) && (n >= 0) && (m >=0) );

    pst_img_graph_t *g = pst_img_graph_create(m,n);

    int32_t i;
    for (uint32_t i = 0; i < n ; i++) { pst_img_graph_read_vertex(wr,g); }
    int32_t *list_next = (int32_t*)malloc(sizeof(int32_t)*2*m);
    int32_t *list_long_bit = (int32_t*)malloc(sizeof(int32_t)*2*m);
    for (uint32_t i = 0; i < m ; i++) { pst_img_graph_read_edge(wr,g,&(list_next[2*i]),&(list_long_bit[2*i])); }

    auto oct_arc_t recover_edge(int32_t index, int32_t long_bit);

    oct_arc_t  recover_edge(int32_t index, int32_t long_bit)
      { assert((index >= 0 ) && (index < g->m));
        oct_arc_t e = g->edge[index].edge;
        assert(e != oct_arc_NULL);
        if (long_bit == 1) { e = oct_sym(e); }
        return e;
      }

    /*Hard part*/
    for (uint32_t i = 0; i < g->m;i++)
      { oct_arc_t e = g->edge[i].edge;
        if (e != oct_arc_NULL)
          { oct_arc_t e_next = recover_edge(list_next[2*i],list_long_bit[2*i]);
            oct_arc_t e_next_sym = recover_edge(list_next[(2*i)+1],list_long_bit[(2*i)+1]);

            oct_splice( oct_oprev(e_next),e);
            oct_splice( oct_oprev(e_next_sym),oct_sym(e));
          }
      }
    free(list_next);
    free(list_long_bit);
    return g;
  }

int32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t *g, r2_t p)
  {
   int32_t i;
   double min_dist_sqr = +INF;
   int32_t min_i = -1;
   for (uint32_t i = 0; i < g->n; i++)
     {
       pst_vertex_data_t *v = &(g->vertex[i]);
       if (v->id == -1) continue;
       double d= r2_dist_sqr(&p,&(v->coords));
       if (d < min_dist_sqr)
         { min_dist_sqr = d; min_i = i; }
     }

    return min_i;
  }

double pst_img_graph_compute_left_face_curl(pst_img_graph_t *g, oct_arc_t e0)
  { oct_arc_t e = e0;
    double curl = 0;
    do{
      curl+= pst_img_graph_get_edge_delta(g,e);
      e = oct_lnext(e);
    }while(e != e0);
    return curl;
  }

r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t *g, oct_arc_t e0)
  {
    oct_arc_t e = e0;
    r2_t p = (r2_t) {{0,0}};
    int32_t count= 0;
    do 
      { int32_t ind_org = pst_img_graph_get_edge_origin(g,e);
        assert(ind_org != -1);
        pst_vertex_data_t *v = &(g->vertex[ind_org]);
        assert(v->id != -1);
        r2_add(&p,&(v->coords),&p);
        count++;
        e = oct_lnext(e);
      } while(e != e0);
    r2_scale(1.0/(double)count,&p,&p);
    return p;
  }

void pst_img_graph_put_curl_into_image(pst_img_graph_t *g,float_image_t *OZ)
  {
    int32_t i;
    for (uint32_t i = 0; i < g->m; i++)
      { oct_arc_t e = g->edge[i].edge;
        if (e != oct_arc_NULL)
          { int32_t j;
            for (uint32_t j = 0; j < 2; j++)
              {
                double c = pst_img_graph_compute_left_face_curl(g,e);
                r2_t p = pst_img_graph_left_face_baricenter(g,e);
                int32_t x,y;
                pst_img_graph_get_vertex_image_indices(&p,OZ->sz[1],OZ->sz[2],&x,&y);
                float_image_set_sample(OZ,0,x,y,c);
                e = oct_sym(e);
              }
          }
      }
  }

void pst_img_graph_free(pst_img_graph_t *g)
{
  int32_t i;
  for (uint32_t i = 0; i < g->m; i++)
{
    oct_arc_t e = g->edge[i].edge;
    if (e != oct_arc_NULL) pst_img_graph_edge_remove(g,e);
    if (g->edge[i].data != NULL)
{
      pst_edge_data_free(g->edge[i].data);
      
    }
  }
  free(g->edge);
  free(g->vertex);
  
  free(g);
}
