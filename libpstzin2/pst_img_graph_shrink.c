/* See {pst_img_graph_shrink.h} */
/* Last edited on 2025-01-13 08:25:59 by stolfi */

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
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <pst_img_graph.h>

#include <pst_img_graph_shrink.h>

#define MARK_VERTEX_NONE 0
#define MARK_VERTEX_REMOVED 1
#define MARK_VERTEX_PRESERVED 2

#define DELETED pst_img_graph_mark_DELETED

#define DEFAULT_WMAG 1.0
 
void pst_img_graph_edge_remove(pst_img_graph_t *g, haf_edge_t e);
  /* Disconnects the edge underlying the arcs {a} and {haf_sym(a)} from
    the graph, and reclaims the record pointed to by the corresponding
    entry of {g.hedge}, which is set to {NULL}. Its {.emark} cannot be
    {DELETED}, and is set to {DELETED} by the procedure. Reclaims its
    label and path. Does not change the IDs of other arcs. */

void pst_img_graph_vertex_remove(pst_img_graph_t *g, uint32_t vi);
  /* Removes the vertex number {vi} from the graph {g}.
    The vertex must be isolated, i.e. it must not be origin or destination
    of any edge, and its {.aout} must be {NULL}. Its {vmark} cannot be 
    {DELETED}, and is set to {DELETED} by the procedure. */ 

void pst_img_graph_mark_vertex_removal(pst_img_graph_t* g, uint32_t degree[]);

void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t* g);
  
void pst_img_graph_remove_paralel_edges(pst_img_graph_t* g);

void pst_img_graph_edge_remove(pst_img_graph_t *g, haf_edge_t e)
  { demand(e != NULL, "invalid null edge {e}");
    
    uint32_t eid = (uint32_t)haf_edge_id(e);
    assert(e == g->hedge[eid]);
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);

    /* Disconnect the edge from he half-edge mesh structure: */
    haf_arc_t a0 = haf_base_arc(e);
    haf_arc_t a1 = haf_sym(a0);

    haf_arc_t oprev_a0 = haf_oprev(a0);
    haf_arc_t oprev_a1 = haf_oprev(a1);
    
    uint32_t org_id = ed->org[0];
    uint32_t dst_id = ed->org[1];

    if (a0 != oprev_a0) { haf_splice(a0, oprev_a0); }
    if (a1 != oprev_a1) { haf_splice(a1, oprev_a1); }

    assert(haf_dprev(a0) == a0);
    assert(haf_dnext(a0) == a0);
    assert(haf_lnext(a0) == a1);
    assert(haf_lprev(a0) == a1);

    g->vdata[org_id].aout = (oprev_a0 == a0 ? NULL: oprev_a0);
    g->vdata[dst_id].aout = (oprev_a1 == a1 ? NULL: oprev_a1);

    haf_edge_free(e);
    g->hedge[eid] = NULL;

    /* Reclaim storage in the edge data record {ed}: */
    ed->org[0] = UINT32_MAX;
    ed->org[1] = UINT32_MAX;
    ed->weight = 0;
    ed->delta = 0;
    free(ed->label); ed->label = NULL;
    pst_path_t *p = &(ed->path);
    if (p->v != NULL) { free(p->v); p->v = NULL; }
    ed->emark = DELETED;
  }

void pst_img_graph_vertex_remove(pst_img_graph_t *g, uint32_t vi)
  { demand(vi < g->NV, "invalid vertex index {vi}");
    pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
    demand(vd->aout == NULL, "vertex is still connected");
    demand(vd->vmark != DELETED, "vertex was already deleted");
    vd->vmark = DELETED;
  }


void debug_vertex_remove_edge(pst_img_graph_t *g, haf_edge_t e)
  {
    uint32_t eid = (uint32_t)haf_edge_id(e);
    haf_arc_t a = haf_base_arc(e);
    double   w   = pst_img_graph_get_edge_weight(g,a);
    double   d   = pst_img_graph_get_arc_delta(g,a);
    uint32_t org = pst_img_graph_get_arc_origin(g,a);
    uint32_t dst = pst_img_graph_get_arc_origin(g,haf_sym(a));
    fprintf(stderr,"ID: %5d ORG: %5d DST: %5d D: %12.6lf W: %12.6lf",eid,org,dst,d,w);
  }

void pst_img_graph_vertex_remove_general
  ( pst_img_graph_t* g,
    uint32_t vi,
    uint32_t n_neighbours,
    double* w_i,
    double wmag,
    bool_t merge_diagonals,
    bool_t verbose
  );

void pst_img_graph_compute_cycle_delta_weight
  ( double wtot,
    uint32_t n_neighbours,
    double *deltas,
    double *weights,
    uint32_t ind_e,
    double *d_cycle,
    double *w_cycle
  );

void pst_img_graph_compute_delta_weight_4
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind,
    double *d_default,
    double *w_default
  );

void pst_img_graph_compute_delta_weight_5
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind ,
    double *d_default,
    double *w_default
  );

void pst_img_graph_compute_delta_weight_6
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind ,
    double *d_default,
    double *w_default
  );

void pst_img_graph_integrate_recursive_vertex_remove
  ( pst_img_graph_t *g,
    uint32_t vi,
    double *w_i,
    double wmag,
    bool_t verbose
  );
  
/* IMPLEMENTATIONS */

void pst_img_graph_mark_vertex_removal(pst_img_graph_t *g, uint32_t degree[])
  {
  //   fprintf(stderr,"\n");
    uint32_t DEG_MAX  = 7;
    uint32_t stats[DEG_MAX];
    for (uint32_t deg = 1; deg <= DEG_MAX; deg++)
      { stats[deg-1] = 0;
        for (uint32_t kv = 0; kv < g->NV; kv++)
          { pst_img_graph_vertex_data_t *vd = &(g->vdata[kv]);
            if (vd->vmark == pst_img_graph_mark_DELETED) { continue; }
            uint32_t n_neighbours = (degree != NULL ? degree[kv]: pst_img_graph_vertex_out_degree(g,kv));
            if ((deg  < DEG_MAX) && (n_neighbours != deg) ) { continue; }
            if (n_neighbours >= DEG_MAX )
              { if (vd->vmark == pst_img_graph_mark_NONE)
                 { vd->vmark = pst_img_graph_mark_TO_KEEP; }
                continue;
              }
            if (vd->vmark == pst_img_graph_mark_NONE)
              { stats[deg-1] = stats[deg-1] + 1;
                vd->vmark = pst_img_graph_mark_REMOVED;
                if ( vd->aout != NULL)
                  { haf_arc_t e = vd->aout;
                    do 
                      { uint32_t dest = pst_img_graph_get_arc_origin(g,haf_sym(e));
                        if (dest == -1)
                          { fprintf(stderr,"ALERT - vertex %d is damaged\n", kv);
                            break;
                          }
                        g->vdata[dest].mark = pst_img_graph_mark_TO_KEEP;
                        e = haf_onext(e);
                      } while(e != vd->aout);
                  }
              }
          }
       /* fprintf(stderr,"Removed %ld Vertices with degree %d.\n",stats[deg-1],deg); */
    }
    /* fprintf(stderr,"\n"); */
  }

void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t *g)
  {
    for (uint32_t kv = 0; kv < g->NV; kv++)
      {
        pst_img_graph_vertex_data_t *vd = &(g->vdata[kv]);
        if (kv == -1) continue;
        if (vd->vmark == pst_img_graph_mark_NONE)
          { vd->vmark = pst_img_graph_mark_REMOVED;
            /* fprintf(stderr,"*"); */
            if ( vd->aout != NULL)
              { haf_arc_t e = vd->aout;
                do
                  { uint32_t dest = pst_img_graph_get_arc_origin(g,haf_sym(e));
                    if (dest == -1)
                      { fprintf(stderr,"ALERT - vertex %d is damaged\n",kv);
                        break;
                       }
                    /* fprintf(stderr,"+"); */
                    g->vdata[dest].mark = pst_img_graph_mark_TO_KEEP;
                    e = haf_onext(e);
                  } while(e != vd->aout);
              }
          }
      }
    /* fprintf(stderr,"\n"); */
  }

void pst_img_graph_remove_paralel_edges(pst_img_graph_t *g)
  {
    for (uint32_t eid = 0; eid < g->NE ; eid++)
      { haf_arc_t e = g->aout[eid].aout;
        if (e == NULL) { continue; }
        for (uint32_t k = 0; k < 2; k++)
          { uint32_t org_e = pst_img_graph_get_arc_origin(g,e);
            uint32_t dst_e = pst_img_graph_get_arc_origin(g,haf_sym(e));
            /* uint32_t ind_e = pst_img_graph_get_dir_edge_num(e); */
            haf_arc_t a = haf_onext(e);
            uint32_t org_a = pst_img_graph_get_arc_origin(g,a);
            uint32_t dst_a = pst_img_graph_get_arc_origin(g,haf_sym(a));
            /* uint32_t ind_a = pst_img_graph_get_dir_edge_num(a); */

            /* fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_e,org_e,dst_e); */
            /* fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_a,org_a,dst_a); */
            if ((a != e) && (org_a == org_e) && (dst_e == dst_a))
              { /* fprintf(stderr,"Paralel edge"); */
                /* eliminar o {e} jogando peso e delta em {a} */
                double de = pst_img_graph_get_arc_delta(g,e);
                double da = pst_img_graph_get_arc_delta(g,a);
                double we = pst_img_graph_get_edge_weight(g,e);
                double wa = pst_img_graph_get_edge_weight(g,a);

                assert(we > 0);
                assert(wa > 0);
                double d = (de*we + da*wa)/(wa+we);
                double w = we + wa;


                pst_img_graph_set_arc_delta(g,a,d);
                pst_img_graph_set_edge_weight(g,a,w);
                pst_img_graph_edge_remove(g,e);
                break;
              }
            e = haf_sym(e);
            /* fprintf(stderr,"\n"); */
          }
      }
  }

void pst_img_graph_compute_cycle_delta_weight
  ( double wtot,
    uint32_t n_neighbours,
    double *deltas,
    double *weights,
    uint32_t ind_e,
    double *d_cycle,
    double *w_cycle
  )
  { double lambda = 0.5;

    auto void compute_wd2l(double w0, double w1, double w2,double *w2l);
    
    void compute_wd2l(double w0, double w1, double w2,double *w2l)
      { *w2l  = w0*(w1+w2)/w1; }

    auto uint32_t update_ind(uint32_t i, int32_t increase);
    
    uint32_t update_ind(uint32_t i, int32_t increase)
      { int32_t res = (int32_t)(i + n_neighbours) + increase;
        assert(res >= 0);
        return ((uint32_t)res)%n_neighbours;
      }
      
    uint32_t ind_f = update_ind(ind_e,+1);

    uint32_t ind_org_act = update_ind(ind_e, -2);
    double w0_org = lambda*(weights[ind_e]*weights[ind_org_act])/wtot;

    uint32_t ind_dst_act = update_ind(ind_f,+2);
    double w0_dst = lambda*(weights[ind_f]*weights[ind_dst_act])/wtot;

    for (uint32_t i = 0; i < n_neighbours - 3; i++)
      { uint32_t ind_org_next = update_ind(ind_org_act, -1);
        double w2_org = (weights[ind_e]*weights[ind_org_next])/wtot;
        double w1_org = (weights[ind_org_act]*weights[ind_org_next])/wtot;

        double w2l;
        compute_wd2l(w0_org,w1_org,w2_org,&w2l);
        w0_org = lambda*w2_org + w2l;

        ind_org_act = ind_org_next;

        uint32_t ind_dst_next = update_ind(ind_dst_act,+1);

        double w2_dst = (weights[ind_f]*weights[ind_dst_next])/wtot;
        double w1_dst = (weights[ind_dst_act]*weights[ind_dst_next])/wtot;
        compute_wd2l(w0_dst,w1_dst,w2_dst,&w2l);

        w0_dst = lambda*w2_dst +  w2l;
        ind_dst_act = ind_dst_next; 
      }
    *w_cycle = w0_org + w0_dst ;
    *d_cycle = -deltas[ind_e] + deltas[ind_f];
  }

void pst_img_graph_compute_delta_weight_4
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind,
    double *d_default,
    double *w_default
  )
  { auto uint32_t next_ind(uint32_t i, uint32_t increase);
  
    uint32_t  next_ind(uint32_t i, uint32_t increase)
      { return (i+increase)%4; }
  
    double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];
    double weight = w_edge[ind]*w_edge[next_ind(ind,+1)] + 0.5*(w_edge[ind]*w_edge[next_ind(ind,+2)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+3)]);
    *d_default = delta;
    *w_default = weight/wtot;
  }

void pst_img_graph_compute_delta_weight_5
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind ,
    double *d_default,
    double *w_default
  )
  { double A = 1.1690;
    double B = 0.4425;
    B = 0;
  
    auto uint32_t next_ind(uint32_t i, uint32_t increase);

    uint32_t next_ind(uint32_t i, uint32_t increase)
      { return (i + 5 + increase)%5; }

    double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];

    double weight = w_edge[ind]*w_edge[next_ind(ind,+1)];
    weight+= A*(w_edge[next_ind(ind,+2)]*w_edge[next_ind(ind,+4)] + w_edge[next_ind(ind,+0)]*w_edge[next_ind(ind,+2)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+4)]);
    weight+= B*(w_edge[next_ind(ind,+0)]*w_edge[next_ind(ind,+3)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+3)] );

    if (weight <= 0) { weight = 0; }

    *d_default = delta;
    *w_default = weight/wtot;
  }

void pst_img_graph_compute_delta_weight_6
  ( double wtot,
    double *d_edge,
    double *w_edge,
    uint32_t ind ,
    double *d_default,
    double *w_default
  )
  {
    double A = 2.0;
    double C = 1.5;

    auto uint32_t next_ind(uint32_t i, int32_t increase);
    
    uint32_t  next_ind(uint32_t i, int32_t increase)
      { int32_t res = (int32_t)i + increase;
        assert(res >= 0);
        return (uint32_t)res%6;
      }

    double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];

    double weight = w_edge[ind]*w_edge[next_ind(ind,+1)];
    weight+= A*(w_edge[next_ind(ind, -1)]*w_edge[next_ind(ind,+2)]);
    weight+= C*(w_edge[next_ind(ind, -1)]*w_edge[next_ind(ind,+1)] + w_edge[next_ind(ind, +0)]*w_edge[next_ind(ind, +2)] );

    *d_default = delta;
    *w_default = weight/wtot;
  }
    
void pst_img_graph_integrate_recursive_vertex_remove
  ( pst_img_graph_t *g,
    uint32_t vi,
    double *w_i,
    double wmag,
    bool_t verbose
  )
  { demand(vi < g->NV, "invalid vertex index {vi}");
    pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
    uint32_t n_neighbours = pst_img_graph_vertex_out_degree(g,vi);
    haf_arc_t a0 = vd->aout;
    if ( w_i != NULL) { a0 = pst_img_graph_find_rightmost_arc(g,a0); }
    if (n_neighbours == 0)
      { /* Does nothing ! */ }
    else if (n_neighbours == 1) 
      { pst_img_graph_edge_remove(g,a0); }
    else if (n_neighbours == 2) 
      { haf_arc_t e1 = haf_onext(a0);

        uint32_t v0 = pst_img_graph_get_arc_origin(g,haf_sym(a0));
        uint32_t v1 = pst_img_graph_get_arc_origin(g,haf_sym(e1));

        double d = pst_img_graph_get_arc_delta(g,haf_sym(a0)) + pst_img_graph_get_arc_delta(g,e1);
        double w0 = pst_img_graph_get_edge_weight(g,a0);
        double w1 = pst_img_graph_get_edge_weight(g,e1);
        double w = 1.0/((1/w0) + (1/w1));
        assert(w > 0);

        pst_path_t p0 = pst_img_graph_get_edge_path(g,haf_sym(a0));
        pst_path_t p1 = pst_img_graph_get_edge_path(g,e1);
        pst_path_t p01 = pst_path_concatenate(p0,g->vdata[vi].coords,p1);

        haf_arc_t a01 = pst_img_graph_edge_add(g,v0,v1,d,w,NULL,p01);
        haf_splice(a01,haf_oprev(haf_sym(a0)));
        haf_splice(haf_sym(a01),haf_sym(e1));

        pst_img_graph_edge_remove(g,a0);
        pst_img_graph_edge_remove(g,e1);
      }
    else if (n_neighbours >=3 )
      { pst_img_graph_vertex_remove_general(g,vi,n_neighbours,w_i,wmag,(n_neighbours >=4),verbose); }
    else
      { assert(FALSE); }

     kv = VERTEX_ID_REMOVED;
     vd->aout = NULL;
  }

void pst_img_graph_vertex_remove_general
  (
    pst_img_graph_t *g,
    uint32_t vi,
    uint32_t n_neighbours,
    double *w_i,
    double wmag,
    bool_t merge_diagonals,
    bool_t verbose
  )
  { assert(n_neighbours >=3);
    if (merge_diagonals) { assert(n_neighbours >=4); }
    pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
    haf_arc_t a0 = vd->aout;
    if (verbose) { fprintf(stderr,"  ---old -------------\n"); }
    haf_arc_t e = a0;
    uint32_t count = 0;
    double wtot = 0;
    double w_edge[n_neighbours];
    double d_edge[n_neighbours];
    do
      { w_edge[count] = pst_img_graph_get_edge_weight(g, e);
        d_edge[count] = pst_img_graph_get_arc_delta(g, e);
        if (verbose)
          { fprintf(stderr,"  %2d ", count);
            debug_vertex_remove_edge(g,e);
            fprintf(stderr,"\n");
          }
        assert(w_edge[count] > 0);
        wtot+= w_edge[count];
        e = haf_onext(e);
        count++;
      } while (e != a0);
    assert(count == n_neighbours);
    if (verbose)
      { fprintf(stderr,"  WTOT = %12.6lf  WTOT^2 = %12.6lf \n",wtot,wtot*wtot); 
        fprintf(stderr,"  ---new -------------\n");
      }
    count = 0;
    e = a0;
    do
      { haf_arc_t f = haf_onext(e);
        double d_default,w_default;
        if (n_neighbours == 4)
          { /* fprintf(stderr,"Removed with 4\n"); */
            pst_img_graph_compute_delta_weight_4(wtot,d_edge,w_edge,count,&d_default,&w_default);
          }
        else if (n_neighbours == 5)
          { /* fprintf(stderr,"Removed with 5\n"); */
            pst_img_graph_compute_delta_weight_5(wtot,d_edge,w_edge,count,&d_default,&w_default);
          }
        else if (n_neighbours == 6)
          { /* fprintf(stderr,"Removed with 6\n"); */
            pst_img_graph_compute_delta_weight_6(wtot,d_edge,w_edge,count,&d_default,&w_default);
          }
        else
          { /* fprintf(stderr,"Removed with %ld\n",n_neighbours); */
            pst_img_graph_compute_cycle_delta_weight(wtot,n_neighbours,d_edge,w_edge, count,&d_default,&w_default);
          }
        double w;
        if (w_i != NULL)
          { w = w_i[count]; }
        else { w = w_default; }
        double d = d_default;
        w = wmag*w;
        if (w != 0)
          { uint32_t ve = pst_img_graph_get_arc_origin(g,haf_sym(e));
            uint32_t vf = pst_img_graph_get_arc_origin(g,haf_sym(f));
            pst_path_t p = pst_img_graph_compute_star_wedge_path(g,vi,ve,vf);
            haf_arc_t a = pst_img_graph_edge_add(g,ve,vf,d,w,NULL,p);
            haf_splice(a, haf_lnext(e));
            haf_splice(haf_sym(a), haf_sym(f));

            if (verbose)
              { fprintf(stderr,"  %2d ",count);
                debug_vertex_remove_edge(g,a);
                fprintf(stderr,"  WDEF = %12.6lf WMAG = %12.9lf",w_default, w/w_default);
                fprintf(stderr,"\n");
              }

            if (haf_lnext(haf_lnext(e)) != haf_sym(f))
              { FILE *wr = open_write("debug.txt", TRUE);
                pst_img_graph_write(wr,g);
                fclose(wr);
              }
            assert(haf_lnext(haf_lnext(e)) == haf_sym(f));
            assert(haf_lnext(e) == a);
          }
        else
          { if (verbose)
              { fprintf(stderr," %2d REMOVED !\n", count);
                fprintf(stderr,"  %2d ", count);
                fprintf(stderr,"  WDEF = %12.6lf WMAG = %12.9lf",w_default, w/w_default);
                fprintf(stderr,"\n");
              }
          }
        e = f;
        count++;
      } while (e != a0);
    e = a0;
    do
      { haf_arc_t f = haf_onext(e);
        pst_img_graph_edge_remove(g,e);
        e = (f == e ? NULL: f);
      } while (e != NULL);
  }

void pst_img_graph_shrink(pst_img_graph_t *g,double *w_i,double wmag)
  {
    /*We will first mark the edges for removal*/
    uint32_t *degree = (uint32_t*)notnull(malloc(sizeof(uint32_t)*(g->NV)),"no mem");
    for (uint32_t i = 0; i < g->NV; i++)
      { g->vdata[i].mark = pst_img_graph_mark_NONE;
        degree[i] = pst_img_graph_vertex_out_degree(g,i);
      }
    pst_img_graph_mark_vertex_removal(g,degree);

    for (uint32_t i = 0; i < g->NV; i++)
      { if (g->vdata[i].mark == pst_img_graph_mark_REMOVED)
          { /* fprintf(stderr,"Removing %ld ",i); */
            pst_img_graph_integrate_recursive_vertex_remove(g,i,w_i,wmag,FALSE);
            /* fprintf(stderr,"\n"); */
          } 
      }
    pst_img_graph_remove_paralel_edges(g);
    for (uint32_t i = 0; i < g->NV; i++)
      { g->vdata[i].mark = pst_img_graph_mark_NONE; }
    free(degree);
  }
???
  
  
  fprintf(stderr,"Solving Graph level[%d] with [%d] vertices [%d] valid\n",level, g->NV,g->NV_valid);

  int32_t *ref_tab;
  pst_imgsys_t *S = pst_img_graph_build_system(g, iW, (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &ref_tab);
  
  if (debug)
{
  
    char *filename = jsprintf("%s-%02d-S.txt",out_prefix,level);
    FILE *wr_dump = open_write(filename,FALSE);
    pst_imgsys_write(wr_dump,S);
    fclose(wr_dump);
    free(filename);
    
    pst_img_graph_put_solution_into_image(g,iZ,OZ);
    filename = jsprintf("%s-%02d-iZ.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
    
    pst_img_graph_put_solution_into_image(g,iW,OZ);
    filename = jsprintf("%s-%02d-iW.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
    
    if (RZ != NULL)
{
      pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
      filename = NULL;
      char *filename = jsprintf("%s-%02d-iE.fni",out_prefix,level);
      wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }
  }
  pst_img_graph_solve_system(g,S,iZ,iW,ref_tab,maxIter,convTol,para,szero,verbose);
  

  if (debug)
{
    
    if (RZ != NULL)
{
      pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
      char *filename = jsprintf("%s-%02d-oE.fni",out_prefix,level);
      FILE *wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }

  }
  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  if (debug)
{
    char *filename = jsprintf("%s-%02d-oZ.fni",out_prefix,level);
    FILE *wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
  }
  free(ref_tab);
  pst_imgsys_free(S);
  
}


void pst_img_graph_edge_remove(pst_img_graph_t *g, haf_edge_t e)
  { demand(e != NULL, "invalid null edge {e}");
    
    uint32_t eid = (uint32_t)haf_edge_id(e);
    assert(e == g->hedge[eid]);
    pst_img_graph_edge_data_t *ed = &(g->edata[eid]);

    /* Disconnect the edge from he half-edge mesh structure: */
    haf_arc_t a0 = haf_base_arc(e);
    haf_arc_t a1 = haf_sym(a0);

    haf_arc_t oprev_a0 = haf_oprev(a0);
    haf_arc_t oprev_a1 = haf_oprev(a1);
    
    uint32_t org_id = ed->org[0];
    uint32_t dst_id = ed->org[1];

    if (a0 != oprev_a0) { haf_splice(a0, oprev_a0); }
    if (a1 != oprev_a1) { haf_splice(a1, oprev_a1); }

    assert(haf_dprev(a0) == a0);
    assert(haf_dnext(a0) == a0);
    assert(haf_lnext(a0) == a1);
    assert(haf_lprev(a0) == a1);

    g->vdata[org_id].aout = (oprev_a0 == a0 ? NULL: oprev_a0);
    g->vdata[dst_id].aout = (oprev_a1 == a1 ? NULL: oprev_a1);

    g->hedge[eid] = NULL;
    pst_img_graph_edge_free(g, e);
  }

void pst_img_graph_vertex_remove(pst_img_graph_t *g, uint32_t vi)
  { demand(vi < g->NV, "invalid vertex index {vi}");
    pst_img_graph_vertex_data_t *vd = &(g->vdata[vi]);
    demand(vd->aout == NULL, "vertex is still connected");
    demand(vd->vmark != DELETED, "vertex was already deleted");
    vd->vmark = DELETED;
  }


