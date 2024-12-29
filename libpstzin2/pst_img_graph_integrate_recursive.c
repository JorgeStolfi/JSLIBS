/* See {pst_img_graph.h} */
/* Last edited on 2024-12-24 19:00:05 by stolfi */
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

#include <pst_imgsys.h>
#include <pst_img_graph.h>
#include <pst_img_graph_integrate.h>

#include <pst_img_graph_integrate_recursive.h>

#define VERTEX_ID_REMOVED UINT32_MAX

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

/* IMPLEMENTATIONS */

void pst_img_graph_mark_vertex_removal(pst_img_graph_t *g, uint32_t degree[])
  {
  //   fprintf(stderr,"\n");
    uint32_t DEG_MAX  = 7;
    uint32_t stats[DEG_MAX];
    for (uint32_t deg = 1; deg <= DEG_MAX; deg++)
      { stats[deg-1] = 0;
        for (uint32_t i = 0; i < g->NV; i++)
          { pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
            if (v->id == -1) continue;
            uint32_t n_neighbours = (degree != NULL ? degree[i]: pst_img_graph_vertex_count_neighbours(g,i));
            if ((deg  < DEG_MAX) && (n_neighbours != deg) ) { continue; }
            if (n_neighbours >= DEG_MAX )
              { if (v->mark == MARK_VERTEX_NONE)
                 { v->mark = MARK_VERTEX_PRESERVED; }
                continue;
              }
            if (v->mark == MARK_VERTEX_NONE)
              { stats[deg-1] = stats[deg-1] + 1;
                v->mark = MARK_VERTEX_REMOVED;
                if ( v->edge != NULL)
                  { haf_arc_t e = v->edge;
                    do 
                      { uint32_t dest = pst_img_graph_get_edge_origin(g,haf_sym(e));
                        if (dest == -1)
                          { fprintf(stderr,"ALERT - vertex %d is damaged\n", v->id);
                            break;
                          }
                        g->vertex[dest].mark = MARK_VERTEX_PRESERVED;
                        e = haf_onext(e);
                      } while(e != v->edge);
                  }
              }
          }
       /* fprintf(stderr,"Removed %ld Vertices with degree %d.\n",stats[deg-1],deg); */
    }
    /* fprintf(stderr,"\n"); */
  }

void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t *g)
  {
    for (uint32_t i = 0; i < g->NV; i++)
      {
        pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (v->id == -1) continue;
        if (v->mark == MARK_VERTEX_NONE)
          { v->mark = MARK_VERTEX_REMOVED;
            /* fprintf(stderr,"*"); */
            if ( v->edge != NULL)
              { haf_arc_t e = v->edge;
                do
                  { uint32_t dest = pst_img_graph_get_edge_origin(g,haf_sym(e));
                    if (dest == -1)
                      { fprintf(stderr,"ALERT - vertex %d is damaged\n",v->id);
                        break;
                       }
                    /* fprintf(stderr,"+"); */
                    g->vertex[dest].mark = MARK_VERTEX_PRESERVED;
                    e = haf_onext(e);
                  } while(e != v->edge);
              }
          }
      }
    /* fprintf(stderr,"\n"); */
  }

void pst_img_graph_remove_paralel_edges(pst_img_graph_t *g)
  {
    for (uint32_t i = 0; i < g->NE ; i++)
      { haf_arc_t e = g->edge[i].edge;
        if (e == NULL) { continue; }
        for (uint32_t k = 0; k < 2; k++)
          { uint32_t org_e = pst_img_graph_get_edge_origin(g,e);
            uint32_t dst_e = pst_img_graph_get_edge_origin(g,haf_sym(e));
            /* uint32_t ind_e = pst_img_graph_get_dir_edge_num(e); */
            haf_arc_t a = haf_onext(e);
            uint32_t org_a = pst_img_graph_get_edge_origin(g,a);
            uint32_t dst_a = pst_img_graph_get_edge_origin(g,haf_sym(a));
            /* uint32_t ind_a = pst_img_graph_get_dir_edge_num(a); */

            /* fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_e,org_e,dst_e); */
            /* fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_a,org_a,dst_a); */
            if ((a != e) && (org_a == org_e) && (dst_e == dst_a))
              { /* fprintf(stderr,"Paralel edge"); */
                /* eliminar o {e} jogando peso e delta em {a} */
                double de = pst_img_graph_get_edge_delta(g,e);
                double da = pst_img_graph_get_edge_delta(g,a);
                double we = pst_img_graph_get_edge_weight(g,e);
                double wa = pst_img_graph_get_edge_weight(g,a);

                assert(we > 0);
                assert(wa > 0);
                double d = (de*we + da*wa)/(wa+we);
                double w = we + wa;


                pst_img_graph_set_edge_delta(g,a,d);
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
    pst_img_graph_vertex_data_t *v = &(g->vertex[vi]);
    haf_arc_t e0 = v->edge;
    if (verbose) { fprintf(stderr,"  ---old -------------\n"); }
    haf_arc_t e = e0;
    uint32_t count = 0;
    double wtot = 0;
    double w_edge[n_neighbours];
    double d_edge[n_neighbours];
    do
      { w_edge[count] = pst_img_graph_get_edge_weight(g, e);
        d_edge[count] = pst_img_graph_get_edge_delta(g, e);
        if (verbose)
          { fprintf(stderr,"  %2d ", count);
            debug_vertex_remove_edge(g,e);
            fprintf(stderr,"\n");
          }
        assert(w_edge[count] > 0);
        wtot+= w_edge[count];
        e = haf_onext(e);
        count++;
      } while (e != e0);
    assert(count == n_neighbours);
    if (verbose)
      { fprintf(stderr,"  WTOT = %12.6lf  WTOT^2 = %12.6lf \n",wtot,wtot*wtot); 
        fprintf(stderr,"  ---new -------------\n");
      }
    count = 0;
    e = e0;
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
          { uint32_t ve = pst_img_graph_get_edge_origin(g,haf_sym(e));
            uint32_t vf = pst_img_graph_get_edge_origin(g,haf_sym(f));
            pst_path_t p = pst_img_graph_compute_star_wedge_path(g,vi,ve,vf);
            haf_arc_t a = pst_img_graph_add_edge(g,ve,vf,d,w,NULL,p);
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
                pst_img_graph_print(wr,g);
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
      } while (e != e0);
    e = e0;
    do
      { haf_arc_t f = haf_onext(e);
        pst_img_graph_edge_remove(g,e);
        e = (f == e ? NULL: f);
      } while (e != NULL);
  }

void pst_img_graph_vertex_remove
  ( pst_img_graph_t *g,
    uint32_t vi,
    double *w_i,
    double wmag,
    bool_t verbose
  )
  { uint32_t n_neighbours = pst_img_graph_vertex_count_neighbours(g,vi);
    pst_img_graph_vertex_data_t *v = &(g->vertex[vi]);
    haf_arc_t e0 = v->edge;
    if ( w_i != NULL) { e0 = pst_img_graph_find_leftmost_edge(g,e0); }
    if (n_neighbours == 0)
      { /* Does nothing ! */ }
    else if (n_neighbours == 1) 
      { pst_img_graph_edge_remove(g,e0); }
    else if (n_neighbours == 2) 
      { haf_arc_t e1 = haf_onext(e0);

        uint32_t v0 = pst_img_graph_get_edge_origin(g,haf_sym(e0));
        uint32_t v1 = pst_img_graph_get_edge_origin(g,haf_sym(e1));

        double d = pst_img_graph_get_edge_delta(g,haf_sym(e0)) + pst_img_graph_get_edge_delta(g,e1);
        double w0 = pst_img_graph_get_edge_weight(g,e0);
        double w1 = pst_img_graph_get_edge_weight(g,e1);
        double w = 1.0/((1/w0) + (1/w1));
        assert(w > 0);

        pst_path_t p0 = pst_img_graph_get_edge_path(g,haf_sym(e0));
        pst_path_t p1 = pst_img_graph_get_edge_path(g,e1);
        pst_path_t p01 = pst_path_concatenate(p0,g->vertex[vi].coords,p1);

        haf_arc_t e01 = pst_img_graph_add_edge(g,v0,v1,d,w,NULL,p01);
        haf_splice(e01,haf_oprev(haf_sym(e0)));
        haf_splice(haf_sym(e01),haf_sym(e1));

        pst_img_graph_edge_remove(g,e0);
        pst_img_graph_edge_remove(g,e1);
      }
    else if (n_neighbours >=3 )
      { pst_img_graph_vertex_remove_general(g,vi,n_neighbours,w_i,wmag,(n_neighbours >=4),verbose); }
    else
      { assert(FALSE); }

     v->id = VERTEX_ID_REMOVED;
     v->edge = NULL;
     g->NV_valid--;
  }

void pst_img_graph_shrink(pst_img_graph_t *g,double *w_i,double wmag)
  {
    /*We will first mark the edges for removal*/
    uint32_t *degree = (uint32_t*)notnull(malloc(sizeof(uint32_t)*(g->NV)),"no mem");
    for (uint32_t i = 0; i < g->NV; i++)
      { g->vertex[i].mark = MARK_VERTEX_NONE;
        degree[i] = pst_img_graph_vertex_count_neighbours(g,i);
      }
    pst_img_graph_mark_vertex_removal(g,degree);

    for (uint32_t i = 0; i < g->NV; i++)
      { if (g->vertex[i].mark == MARK_VERTEX_REMOVED)
          { /* fprintf(stderr,"Removing %ld ",i); */
            pst_img_graph_vertex_remove(g,i,w_i,wmag,FALSE);
            /* fprintf(stderr,"\n"); */
          } 
      }
    pst_img_graph_remove_paralel_edges(g);
    for (uint32_t i = 0; i < g->NV; i++)
      { g->vertex[i].mark = MARK_VERTEX_NONE; }
    free(degree);
  }
  
void pst_img_graph_copy_solution_from_shrunk(pst_img_graph_t *jg,double *jZ,pst_img_graph_t *ig, double *iZ)
  {
    assert(ig->NV >= jg->NV);
    uint32_t j = 0;
    for (uint32_t i = 0; i < ig->NV; i++)
      { pst_img_graph_vertex_data_t *vIG = &(ig->vertex[i]);
        if ( vIG->id == -1)
          { iZ[i] = 0;
            continue;
          }
      pst_img_graph_vertex_data_t *vJG = NULL;
      while(j < jg->NV)
        { vJG = &(jg->vertex[j]);
          if ( vJG->id != -1) break;
          j++;
        }

      if ( (j >= jg->NV) || (vJG->id > vIG->id) )
        { iZ[i] = 0; vIG->mark = MARK_VERTEX_REMOVED; 
        }
      else if (vJG->id == vIG->id)
        { iZ[i] = jZ[j]; vIG->mark = MARK_VERTEX_PRESERVED;
          j++;
        }else{
        demand(FALSE,"Vertices out of order!");
      }
    }
  }

void pst_img_graph_estimate_from_shrunk(pst_img_graph_t *jg,double *jZ,pst_img_graph_t *ig, double *iZ)
  {
    assert(ig->NV >= jg->NV);
    pst_img_graph_copy_solution_from_shrunk(jg,jZ,ig,iZ);
    for (uint32_t i = 0; i < ig->NV; i++)
      { pst_img_graph_vertex_data_t *vIG = &(ig->vertex[i]);
        if (vIG->id == -1) continue;
        if (vIG->mark == MARK_VERTEX_REMOVED)
          { haf_arc_t e0 = vIG->edge;
            if (e0 == NULL) continue;
            haf_arc_t e =  e0;
            double sW = 0;
            double sWZ = 0;
            do
              { uint32_t ind_dst = pst_img_graph_get_edge_origin(ig,haf_sym(e));
                assert(ig->vertex[ind_dst].mark == MARK_VERTEX_PRESERVED);
                double d = pst_img_graph_get_edge_delta(ig,e);
                double w = pst_img_graph_get_edge_weight(ig,e);
                sW+= w;
                sWZ+=w*(iZ[ind_dst]-d);
                e = haf_onext(e);
              } while(e0 != e);
            assert(sW > 0);
            iZ[i] = sWZ/sW;
          }
      }
  }

void pst_img_graph_integration_recursive
  ( 
    pst_img_graph_t *g,
    double *iZ,
    double *iW,
    double wmag,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose,
    uint32_t level,
    float_image_t *OZ, /*debug only*/
    float_image_t *RZ,
    char *out_prefix,
    bool_t debug
  )
  {
  fprintf(stderr,"Starting level %d\n",level);
  
  if (debug)
{
    pst_img_graph_put_curl_into_image(g,OZ);
    char *filename = jsprintf("%s-%02d-CL.fni",out_prefix,level);
    FILE *wr_curl = open_write(filename,FALSE);
     float_image_write(wr_curl,OZ);
    fclose(wr_curl);
    free(filename);
  }
  
  if (debug)
{
    pst_img_graph_mark_vertex_removal(g,NULL);
    char *filename = jsprintf("%s-%02d-GR.txt",out_prefix,level);
    FILE *wr_grph = open_write(filename,FALSE);
     pst_img_graph_print(wr_grph,g);
    fclose(wr_grph);
    free(filename);
    filename = jsprintf("%s-%02d-GR.grf",out_prefix,level);
    wr_grph = open_write(filename,FALSE);
     pst_img_graph_write(wr_grph,g);
    fclose(wr_grph);
    free(filename);
  }
  
  if ( g->NV_valid > 2)
{
    fprintf(stderr,"Reducing Graph level[%d] with [%d] vertices [%d] valid\n",level,g->NV,g->NV_valid);
    fprintf(stderr,"Copying...");
    pst_img_graph_t *jg = pst_img_graph_copy(g);
/* fprintf(stderr,"OK.\nShrinking..."); */
    pst_img_graph_shrink(jg,NULL,wmag);
/* fprintf(stderr,"OK.\nNext\n"); */

    double *jZ = rn_alloc(jg->NV);
    double *jW = rn_alloc(jg->NV);
    
/* pst_img_graph_integration_recursive(jg,jZ,jW,wmag,(uint32_t)ceil(maxIter*M_SQRT2),convTol/M_SQRT2,para,szero,verbose,level+1,OZ,RZ,out_prefix,debug); */
    double ratio = (((double)g->NV_valid)/((double)jg->NV_valid));
    double newConvTol = convTol/sqrt(ratio);
    uint32_t newmaxIter = (uint32_t)ceil(sqrt(ratio)*(double)maxIter);
    fprintf(stdout,"%9.6lf\n",ratio);
    fprintf(stderr,"Ratio %9.6lf G_valid %d JG_valid %d New ConvTol %lf NewmaxIter %d \n",ratio,g->NV_valid,jg->NV_valid, newConvTol,newmaxIter);
    pst_img_graph_integration_recursive(jg,jZ,jW,wmag,newmaxIter,newConvTol,para,szero,verbose,level+1,OZ,RZ,out_prefix,debug);
    pst_img_graph_estimate_from_shrunk(jg,jZ,g,iZ);
    /* The vertices must be painted now*/
    free(jZ);
    free(jW);
  }else{
    fprintf(stderr,"End of Recursion - returning\n");
    for (uint32_t i = 0; i < g->NV;i++) iZ[i] = 0; 
  }
  
  
  
  
  fprintf(stderr,"Solving Graph level[%d] with [%d] vertices [%d] valid\n",level,g->NV,g->NV_valid);

  int32_t *ref_tab;
  pst_imgsys_t *S = pst_img_graph_build_integration_system(g, iW, (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &ref_tab);
  
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

