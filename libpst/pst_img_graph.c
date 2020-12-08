/* See {pst_img_graph.h} */
/* Last edited on 2013-10-02 03:23:42 by stolfilocal */
/* Created by Rafael F. V. Saracchini */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <oct.h>
#include <float_image.h>
#include <jsfile.h>
#include <affirm.h>

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

long int pst_img_graph_get_vertex_index_from_image_indices(long int ix, long int iy, long int NX, long int NY)
  { if ((iy < 0) || (iy > NY)) { return -1; }
    if ((ix < 0) || (ix > NX)) { return -1; }
    return iy*NX + ix;
  }

void pst_img_graph_get_vertex_image_indices(r2_t *p,long int NX, long int NY, long int *ix, long int *iy)
  { (*ix) = (long int)floor(p->c[0] + 0.5);
    (*iy) = (long int)floor(p->c[1] + 0.5);
    if ((*ix) >= NX ) { (*ix) = NX -1; }
    if ((*iy) >= NY ) { (*iy) = NY -1; }
    if ((*ix) < 0 ) { (*ix) = 0; }
    if ((*iy) < 0 ) { (*iy) = 0; }
  }

long int  pst_img_graph_vertex_add(pst_img_graph_t *g, long int id,oct_arc_t edge,r2_t coords)
  {  assert(g->n < g->max_n);
    long int ix = g->n;
    pst_vertex_data_t *v = &(g->vertex[ix]);

    v->id = id;
    v->mark = MARK_VERTEX_NONE;
    v->edge = edge;
    v->coords = coords;
    g->n++;
    g->n_valid++;
    return ix;
  }

long int pst_img_graph_get_dir_edge_num(oct_arc_t e)
  { if (e == oct_arc_NULL) { return  -1; }
    return (2*oct_edge_num(oct_edge(e))) + oct_lon_bit(e);
  }

long int pst_img_graph_get_edge_num(oct_arc_t e)
  { if (e == oct_arc_NULL) { return -1; }
    return oct_edge_num(oct_edge(e));
  }

long int pst_img_graph_get_edge_origin(pst_img_graph_t *g, oct_arc_t e)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->org[lbit];
  }

void pst_img_graph_set_edge_origin(pst_img_graph_t *g, oct_arc_t e,long int org)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind < (2*g->max_m)));
    g->edge[ind].data->org[lbit] = org;
  }

double pst_img_graph_get_edge_weight(pst_img_graph_t *g, oct_arc_t e)
  { long int ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->max_m));
    return g->edge[ind].data->weight;
  }

void pst_img_graph_set_edge_weight(pst_img_graph_t *g, oct_arc_t e,double w)
  { long int ind = pst_img_graph_get_edge_num(e);
    assert((ind >= 0) && (ind < g->max_m));
    assert((w == 0) || (w > 1.0e-100));
    g->edge[ind].data->weight = w;
  }

double pst_img_graph_get_edge_delta(pst_img_graph_t *g, oct_arc_t e)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->delta[lbit];
  }

void pst_img_graph_set_edge_delta(pst_img_graph_t *g,oct_arc_t e,double delta)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    g->edge[ind].data->delta[lbit] = delta;
    g->edge[ind].data->delta[!lbit] = -delta;
  }

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t *g, oct_arc_t e)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
    assert((ind_dir >= 0) && (ind_dir < (2*g->max_m)));
    return g->edge[ind].data->path[lbit];
  }

void pst_img_graph_set_edge_path(pst_img_graph_t *g,oct_arc_t e,pst_path_t p)
  { long int ind_dir = pst_img_graph_get_dir_edge_num(e);
    long int ind = pst_img_graph_get_edge_num(e);
    int lbit = ind_dir&1;
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

r2_t pst_path_get_vertex(pst_path_t p, long int i)
  { return ( p.reverse ? p.v[p.n - i -1 ] : p.v[i] );  }

pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1)
  { int n0 = p0.n;
    int n1 = p1.n;
    int n = n0 + n1 + 1;

    r2_t *v = (r2_t*)malloc(sizeof(r2_t)*n);
    long int i;
    for (i = 0; i < n; i++)
      { r2_t c;
        if (i < n0)
          { c = pst_path_get_vertex(p0,i); }
        else if (i == n0)
          { c = coords; }
        else
          { long int j = i - n0 - 1;
            c = pst_path_get_vertex(p1,j);
          }
        v[i] = c;
      }
    return (pst_path_t) { .n=n, .v=v, .reverse=FALSE };
  }

pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t *g,long int vi, long int ve, long int vf)
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
    long int org,long int dst,
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
    long int ind = pst_img_graph_get_edge_num(e);
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
  { long int org_e = pst_img_graph_get_edge_origin(g,e);
    long int dst_e = pst_img_graph_get_edge_origin(g,oct_sym(e));

    oct_arc_t a = oct_oprev(e);
    oct_arc_t b = oct_oprev(oct_sym(e));

    pst_img_graph_set_edge_origin(g,e,-1);
    pst_img_graph_set_edge_origin(g,oct_sym(e),-1);
    pst_img_graph_set_edge_weight(g,e,0);
    pst_img_graph_set_edge_delta(g,e,0);

    pst_path_t p = pst_img_graph_get_edge_path(g,e);
    if (p.v != NULL) { free(p.v); }

    long int ind = pst_img_graph_get_edge_num(e);
    long int ind_e0 = pst_img_graph_get_dir_edge_num(e);
    long int ind_e1 = pst_img_graph_get_dir_edge_num(oct_sym(e));
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

pst_img_graph_t *pst_img_graph_create(long int max_m,long int max_n)
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
    int i;
    for (i = 0; i < max_m;i++) { g->edge[i].edge = oct_arc_NULL; g->edge[i].data = NULL; }
    for (i = 0; i < max_n; i++)
      { g->vertex[i].edge = oct_arc_NULL;
        g->vertex[i].id = -1;
        g->vertex[i].mark = MARK_VERTEX_NONE;
      }

    return g;
  }

void pst_img_graph_print_vertex(FILE *wr, pst_vertex_data_t *v)
  { long int ind_edge_dir = pst_img_graph_get_dir_edge_num(v->edge);
    long int ind_edge = pst_img_graph_get_edge_num(v->edge);
    fprintf(wr,"ID: %ld (%f,%f)  EDGE: %ld (%ld) \n",v->id,v->coords.c[0],v->coords.c[1],ind_edge_dir, ind_edge);
  }

void pst_img_graph_print_edge(FILE *wr,pst_img_graph_t *g, oct_arc_t e)
  { if (e == oct_arc_NULL) 
      { fprintf(wr," EMPTY "); }
    else
      {  long int ind = pst_img_graph_get_edge_num(e);
        double   w   = pst_img_graph_get_edge_weight(g,e);
        double   d   = pst_img_graph_get_edge_delta(g,e);
        double   rd  = pst_img_graph_get_edge_delta(g,oct_sym(e));
        long int org = pst_img_graph_get_edge_origin(g,e);
        long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
        pst_path_t p = pst_img_graph_get_edge_path(g,e);
        char *label = (g->edge[ind].data->label == NULL ? "": g->edge[ind].data->label);
        fprintf(wr,"ID: %ld ORG: %ld DEST: %ld Delta: %lf RDelta %lf W: %lf LABEL %s PATH: %ld ",ind,org,dst,d,rd,w,label,p.n);
        /*Neighbours*/
        fprintf(wr,"ONEXTS ");
        oct_arc_t oe;
        for (oe = oct_onext(e); oe != e; oe = oct_onext(oe))
          { fprintf(wr, "%ld:%ld ",pst_img_graph_get_edge_num(oe),pst_img_graph_get_dir_edge_num(oe)&1); }
      }
    fprintf(wr,"\n");
  }

void pst_img_graph_print(FILE *wr,pst_img_graph_t *g)
{
  long int i;
  for(i = 0; i < g->n; i++)
{
    fprintf(wr,"[%ld] ",i);
    pst_img_graph_print_vertex(wr,&(g->vertex[i]));
  }
  fprintf(wr,"\n");
  for(i = 0; i < g->m; i++)
{
    fprintf(wr,"[%ld] ",i);
    pst_img_graph_print_edge(wr,g,g->edge[i].edge);
  }
  fprintf(wr,"\n");
}

void pst_img_graph_mark_vertex_removal(pst_img_graph_t *g, int degree[])
{
  int i;
//   fprintf(stderr,"\n");
  int DEG_MAX  = 7;
  long int stats[DEG_MAX];
  int deg;
  for(deg = 1; deg <= DEG_MAX; deg++)
{
    stats[deg-1] = 0;
    for(i = 0; i < g->n; i++)
{
      pst_vertex_data_t *v = &(g->vertex[i]);
      
      if(v->id == -1) continue;
      int n_neighbours = (degree != NULL ? degree[i]: pst_img_graph_vertex_count_neighbours(g,i));
      if((deg  < DEG_MAX) && (n_neighbours != deg) ) continue;
      if(n_neighbours >= DEG_MAX )
{
        if(v->mark == MARK_VERTEX_NONE)
{
          v->mark = MARK_VERTEX_PRESERVED;
        }
        continue;
      }
      if(v->mark == MARK_VERTEX_NONE)
{
        stats[deg-1] = stats[deg-1] + 1;
        v->mark = MARK_VERTEX_REMOVED;

        if( v->edge != oct_arc_NULL)
{
          oct_arc_t e = v->edge;
          do{
            
            long int dest = pst_img_graph_get_edge_origin(g,oct_sym(e));
            if(dest == -1)
{
              fprintf(stderr,"ALERT - vertex %ld is damaged\n",v->id);
              
              break;
            }

            g->vertex[dest].mark = MARK_VERTEX_PRESERVED;
            e = oct_onext(e);
          }while(e != v->edge);
        }
      }
    }
//       fprintf(stderr,"Removed %ld Vertices with degree %d.\n",stats[deg-1],deg);
  }

//   fprintf(stderr,"\n");
}

void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t *g)
{
  int i;
//   fprintf(stderr,"\n");
  
    for(i = 0; i < g->n; i++)
{
      pst_vertex_data_t *v = &(g->vertex[i]);
      
      if(v->id == -1) continue;
       if(v->mark == MARK_VERTEX_NONE)
{
        v->mark = MARK_VERTEX_REMOVED;
  //       fprintf(stderr,"*");
        if( v->edge != oct_arc_NULL)
{
          oct_arc_t e = v->edge;
          do{
            
            long int dest = pst_img_graph_get_edge_origin(g,oct_sym(e));
            if(dest == -1)
{
              fprintf(stderr,"ALERT - vertex %ld is damaged\n",v->id);
              
              break;
            }
  //      fprintf(stderr,"+");
            g->vertex[dest].mark = MARK_VERTEX_PRESERVED;
            e = oct_onext(e);
          }while(e != v->edge);
        }
      }
    }
//   fprintf(stderr,"\n");
}




long int pst_img_graph_vertex_count_neighbours(pst_img_graph_t *g, long int vi)
{
  
  pst_vertex_data_t *v = &(g->vertex[vi]);
  
  if(v->edge == oct_arc_NULL) return 0;
  long int count = 1;
  oct_arc_t e = v->edge;
  for(e = oct_onext(v->edge); e!= v->edge; e = oct_onext(e))
{
    count++;
  }
  return count;
}

oct_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t *g,long int vi0, long int vi1)
{
  pst_vertex_data_t *v0 = &(g->vertex[vi0]);
  if( (v0->edge == oct_arc_NULL) || (v0->edge == oct_arc_NULL) )
{  return oct_arc_NULL;}
  oct_arc_t e = v0->edge;
  do
  {
    long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
    if(dst == vi1) return e;
    e = oct_onext(e);
  }while(e != v0->edge);
  
  return oct_arc_NULL;
}

void pst_img_graph_remove_paralel_edges(pst_img_graph_t *g)
{
  long int i;
  for(i = 0; i < g->m ; i++)
{
    oct_arc_t e = g->edge[i].edge;
    
    if(e == oct_arc_NULL) continue;
    int k;
    for(k = 0; k < 2; k++)
{
      
      long int org_e = pst_img_graph_get_edge_origin(g,e);
      long int dst_e = pst_img_graph_get_edge_origin(g,oct_sym(e));
//       long int ind_e = pst_img_graph_get_dir_edge_num(e);
    
      oct_arc_t a = oct_onext(e);
    
      long int org_a = pst_img_graph_get_edge_origin(g,a);
      long int dst_a = pst_img_graph_get_edge_origin(g,oct_sym(a));
//       long int ind_a = pst_img_graph_get_dir_edge_num(a);
    
//       fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_e,org_e,dst_e);
//       fprintf(stderr,"Edge %ld (%ld - %ld) ",ind_a,org_a,dst_a);
      if((a != e) && (org_a == org_e) && (dst_e == dst_a))
{
      
//      fprintf(stderr,"Paralel edge");
        //eliminar o {e} jogando peso e delta em {a}
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
      e = oct_sym(e);
//       fprintf(stderr,"\n");
    }
   
  }
}

// void pst_img_graph_remove_paralel_edges(pst_img_graph_t *g)
//  {
//   long int i;
//   for(i = 0; i < g->m ; i++)
//     {
//     oct_arc_t e = g->edge[i];
//     if(e == oct_arc_NULL) continue;
//     long int dst_e = pst_img_graph_get_edge_origin(g,oct_sym(e));
//     oct_arc_t a = oct_onext(e);
//     long int count_repeated = 0;
//     double sumW = pst_img_graph_get_edge_weight(g,e);
//     double sumD = pst_img_graph_get_edge_delta(g,e)*sumW;
//     while(a != e)
//       {
//       long int dst_a = pst_img_graph_get_edge_origin(g,oct_sym(a));
//       oct_arc_t a_next = oct_onext(a);
//       if(dst_a == dst_e)
//         {
//      double wa = pst_img_graph_get_edge_weight(g,a);
//      double da = pst_img_graph_get_edge_delta(g,a);
//      sumW+=wa;
//      sumD+=wa*da;
//      count_repeated++;
//      pst_img_graph_edge_remove(g,a);
//       }
//       a = a_next;
//     }
//     if(count_repeated > 0)
//       {  
//       pst_img_graph_set_edge_delta(g,a,sumD);
//       pst_img_graph_set_edge_weight(g,a,sumW);
//       
//     }
//   }
// }

oct_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t *g,oct_arc_t e0)
{
  oct_arc_t e = e0;
  assert(e0 != oct_arc_NULL);
  long int org = pst_img_graph_get_edge_origin(g,e0);
  pst_vertex_data_t *v_org = &(g->vertex[org]);
  oct_arc_t le = oct_arc_NULL;
  double cmax = -INF;
  do{
    long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
    pst_vertex_data_t *v_dst = &(g->vertex[dst]);
    r2_t d;
    r2_sub(&(v_dst->coords),&(v_org->coords),&d);
    (void) r2_dir(&d,&d);
    double c = d.c[0];
    if(c > cmax)
{    le = e;     cmax = c;    }
    e = oct_onext(e);
  }while(e != e0);
  
  return le;
}


void debug_vertex_remove_edge(pst_img_graph_t *g, oct_arc_t e)
{
    long int ind = pst_img_graph_get_edge_num(e);
    double   w   = pst_img_graph_get_edge_weight(g,e);
    double   d   = pst_img_graph_get_edge_delta(g,e);
    long int org = pst_img_graph_get_edge_origin(g,e);
    long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
    fprintf(stderr,"ID: %5ld ORG: %5ld DST: %5ld D: %12.6lf W: %12.6lf",ind,org,dst,d,w);
    
}

void pst_img_graph_compute_cycle_delta_weight(
  double wtot,
  long int n_neighbours,
  double *deltas,
  double *weights,
  long int ind_e,
  double *d_cycle,
  double *w_cycle
  )
{
  
  double lambda = 0.5;

  auto void compute_wd2l(double w0, double w1, double w2,double *w2l);
  void compute_wd2l(double w0, double w1, double w2,double *w2l)
    {
      *w2l  = w0*(w1+w2)/w1;
    }
  
  auto long int update_ind(long int i,long int increase);
  long int update_ind(long int i,long int increase)
    {
      return (i + n_neighbours + increase)%n_neighbours;
    }
  long int ind_f = update_ind(ind_e,+1);
  
  long int ind_org_act = update_ind(ind_e,-2);
  double w0_org = lambda*(weights[ind_e]*weights[ind_org_act])/wtot;
  
  long int ind_dst_act = update_ind(ind_f,+2);
  double w0_dst = lambda*(weights[ind_f]*weights[ind_dst_act])/wtot;

  long int i;
  
  for(i = 0; i < n_neighbours - 3; i++)
{
    
    long int ind_org_next = update_ind(ind_org_act,-1);
    
    
    double w2_org = (weights[ind_e]*weights[ind_org_next])/wtot;
    double w1_org = (weights[ind_org_act]*weights[ind_org_next])/wtot;

    double w2l;
    compute_wd2l(w0_org,w1_org,w2_org,&w2l);
    w0_org = lambda*w2_org + w2l;
    
    ind_org_act = ind_org_next;
    

    long int ind_dst_next = update_ind(ind_dst_act,+1);
    
    double w2_dst = (weights[ind_f]*weights[ind_dst_next])/wtot;
    double w1_dst = (weights[ind_dst_act]*weights[ind_dst_next])/wtot;
    compute_wd2l(w0_dst,w1_dst,w2_dst,&w2l);

    w0_dst = lambda*w2_dst +  w2l;
    ind_dst_act = ind_dst_next;
    
  }
  *w_cycle = w0_org + w0_dst ;
  *d_cycle = -deltas[ind_e] + deltas[ind_f];
}

void pst_img_graph_compute_delta_weight_4(double wtot, double *d_edge,double *w_edge,long int ind ,double *d_default,double *w_default)
{
  auto long int next_ind(long int i, int increase);
  long int  next_ind(long int i, int increase)
{
    return (i+increase)%4;
  }
  
  double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];
  double weight = w_edge[ind]*w_edge[next_ind(ind,+1)] + 0.5*(w_edge[ind]*w_edge[next_ind(ind,+2)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+3)]);
  *d_default = delta;
  *w_default = weight/wtot;
}


void pst_img_graph_compute_delta_weight_5(double wtot, double *d_edge,double *w_edge,long int ind ,double *d_default,double *w_default)
{
  double A = 1.1690;
  double B = 0.4425;
  B = 0;
  
  auto long int next_ind(long int i, int increase);
  long int  next_ind(long int i, int increase)
{
    return (i + 5 + increase)%5;
  }
  
  double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];
  
  double weight = w_edge[ind]*w_edge[next_ind(ind,+1)];
  weight+= A*(w_edge[next_ind(ind,+2)]*w_edge[next_ind(ind,+4)] + w_edge[next_ind(ind,+0)]*w_edge[next_ind(ind,+2)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+4)]);
  weight+= B*(w_edge[next_ind(ind,+0)]*w_edge[next_ind(ind,+3)] + w_edge[next_ind(ind,+1)]*w_edge[next_ind(ind,+3)] );
  
  if(weight <= 0)
{
    weight = 0;
  }
    
  *d_default = delta;
  *w_default = weight/wtot;
}

void pst_img_graph_compute_delta_weight_6(double wtot, double *d_edge,double *w_edge,long int ind ,double *d_default,double *w_default)
{
  double A = 2.0;
  double B = 0;
  double C = 1.5;
  double D = 0;
  
  auto long int next_ind(long int i, int increase);
  long int  next_ind(long int i, int increase)
{
    return (i + 6 + increase)%6;
  }
  
  double delta = -d_edge[ind] + d_edge[next_ind(ind,+1)];
  
  double weight = w_edge[ind]*w_edge[next_ind(ind,+1)];
  weight+= A*(w_edge[next_ind(ind,-1)]*w_edge[next_ind(ind,+2)]);
  weight+= C*(w_edge[next_ind(ind,-1)]*w_edge[next_ind(ind,+1)] + w_edge[next_ind(ind,+0)]*w_edge[next_ind(ind,+2)] );
  
  *d_default = delta;
  *w_default = weight/wtot;
}


void pst_img_graph_vertex_remove_general(
  pst_img_graph_t *g,
  long int vi,
  long int n_neighbours,
  double *w_i,
  double wmag,
  bool_t merge_diagonals,
  bool_t verbose)
 {
   
  assert(n_neighbours >=3);
  if(merge_diagonals)
{ assert(n_neighbours >=4); }
      pst_vertex_data_t *v = &(g->vertex[vi]);
      oct_arc_t e0 = v->edge;
    if(verbose)
{
      fprintf(stderr,"  ---old -------------\n");
    }
    oct_arc_t e = e0;
    long int count = 0;
    double wtot = 0;
    double w_edge[n_neighbours];
    double d_edge[n_neighbours];
    do{
      w_edge[count] = pst_img_graph_get_edge_weight(g,e);
      d_edge[count] = pst_img_graph_get_edge_delta(g,e);
      if(verbose)
{
        fprintf(stderr,"  %2ld ",count);
        debug_vertex_remove_edge(g,e);
        fprintf(stderr,"\n");
      }
      assert(w_edge[count] > 0);
      wtot+= w_edge[count];
      e = oct_onext(e);
      count++;
    }while (e != e0);
    assert(count == n_neighbours);
    if(verbose)
{
      fprintf(stderr,"  WTOT = %12.6lf  WTOT^2 = %12.6lf \n",wtot,wtot*wtot); 
      fprintf(stderr,"  ---new -------------\n");
    }
    count = 0;
    e = e0;
    do{
      oct_arc_t f = oct_onext(e);
//       d = d/(wtot*wtot);
//       double w;
//       double we = w_edge[count];
//       double wf = w_edge[(count+1)%n_neighbours];
//       assert((we  > 0) && (wf > 0));
//       double w_default =  (we*wf)/(wtot);
//       if(merge_diagonals)
{
//      double wd = w_edge[(count-1 + n_neighbours)%n_neighbours];
//      double wg = w_edge[(count+2)%n_neighbours];
//      w_default*= (1+ 0.5*(wg/wtot)*((we/wf) +1.0) + 0.5*(wd/wtot)*((wf/we)+1.0)  );
//       }
      double d_default,w_default;
      if(n_neighbours == 4)
{
//      fprintf(stderr,"Removed 4\n");
        pst_img_graph_compute_delta_weight_4(wtot,d_edge,w_edge,count,&d_default,&w_default);
      }else if(n_neighbours == 5)
{
//      fprintf(stderr,"Removed 6\n");
        pst_img_graph_compute_delta_weight_5(wtot,d_edge,w_edge,count,&d_default,&w_default);
      }else if(n_neighbours == 6)
{
//      fprintf(stderr,"Removed 6\n");
        pst_img_graph_compute_delta_weight_6(wtot,d_edge,w_edge,count,&d_default,&w_default);
      }else{
//      fprintf(stderr,"Removed %ld\n",n_neighbours);
        pst_img_graph_compute_cycle_delta_weight(wtot,n_neighbours,d_edge,w_edge, count,&d_default,&w_default);
      }
      double w;
      if(w_i != NULL)
{
        w = w_i[count];
      }else {
        w = w_default;
      }
      double d = d_default;
      w = wmag*w;
      
      if(w != 0)
{
      
        long int ve = pst_img_graph_get_edge_origin(g,oct_sym(e));
        long int vf = pst_img_graph_get_edge_origin(g,oct_sym(f));
      
        pst_path_t p = pst_img_graph_compute_star_wedge_path(g,vi,ve,vf);
  
        oct_arc_t a = pst_img_graph_edge_insert(g,ve,vf,d,w,NULL,p);
        oct_splice(a, oct_lnext(e));
        oct_splice(oct_sym(a), oct_sym(f));
        
        if(verbose) {
          fprintf(stderr,"  %2ld ",count);
          debug_vertex_remove_edge(g,a);
          fprintf(stderr,"  WDEF = %12.6lf WMAG = %12.9lf",w_default, w/w_default);
          fprintf(stderr,"\n");
        }
        
        if (oct_lnext(oct_lnext(e)) != oct_sym(f))
{
          FILE *wr = open_write("debug.txt",TRUE);
          pst_img_graph_print(wr,g);
          fclose(wr);
        }
        assert(oct_lnext(oct_lnext(e)) == oct_sym(f));
        assert(oct_lnext(e) == a);
      }else{
        if(verbose) {
          fprintf(stderr," %2ld REMOVED !\n",count);
          fprintf(stderr,"  %2ld ",count);
          fprintf(stderr,"  WDEF = %12.6lf WMAG = %12.9lf",w_default, w/w_default);
          fprintf(stderr,"\n");
        }
      }
      e = f;
      count++;
    }while (e != e0);
    
    e = e0;
    do{
      oct_arc_t f = oct_onext(e);
      pst_img_graph_edge_remove(g,e);
      e = (f == e ? oct_arc_NULL: f);
    }while (e != oct_arc_NULL);
}


void pst_img_graph_vertex_remove(pst_img_graph_t *g, long int vi,double *w_i,double wmag,bool_t verbose)
{
  
  
  
  long int n_neighbours = pst_img_graph_vertex_count_neighbours(g,vi);
  pst_vertex_data_t *v = &(g->vertex[vi]);
  oct_arc_t e0 = v->edge;
  if( w_i != NULL)
{
    e0 = pst_img_graph_find_leftmost_edge(g,e0);
  }
  
  if(n_neighbours == 0)
{
    /* Does nothing ! */
  }
  else if(n_neighbours == 1) {
     pst_img_graph_edge_remove(g,e0);
  }
  else if(n_neighbours == 2) {
    
    oct_arc_t e1 = oct_onext(e0);
    
    long int v0 = pst_img_graph_get_edge_origin(g,oct_sym(e0));
    long int v1 = pst_img_graph_get_edge_origin(g,oct_sym(e1));
   
    double d = pst_img_graph_get_edge_delta(g,oct_sym(e0)) + pst_img_graph_get_edge_delta(g,e1);
    double w0 = pst_img_graph_get_edge_weight(g,e0);
    double w1 = pst_img_graph_get_edge_weight(g,e1);
    double w = 1.0/((1/w0) + (1/w1));
    assert(w > 0);
    
    pst_path_t p0 = pst_img_graph_get_edge_path(g,oct_sym(e0));
    pst_path_t p1 = pst_img_graph_get_edge_path(g,e1);
    pst_path_t p01 = pst_path_concatenate(p0,g->vertex[vi].coords,p1);
    
    oct_arc_t e01 = pst_img_graph_edge_insert(g,v0,v1,d,w,NULL,p01);
    oct_splice(e01,oct_oprev(oct_sym(e0)));
    oct_splice(oct_sym(e01),oct_sym(e1));
    
    pst_img_graph_edge_remove(g,e0);
    pst_img_graph_edge_remove(g,e1);
    
  }
  else if (n_neighbours >=3 )
{
     pst_img_graph_vertex_remove_general(g,vi,n_neighbours,w_i,wmag,(n_neighbours >=4),verbose);
  }else{
    assert(FALSE);
  }
  
   v->id = -1;
   v->edge = oct_arc_NULL;
   g->n_valid--;
}

void pst_img_graph_shrink(pst_img_graph_t *g,double *w_i,double wmag)
{
  /*We will first mark the edges for removal*/
  int *degree = (int*)notnull(malloc(sizeof(int)*(g->n)),"no mem");
  long int i;
  for(i = 0; i < g->n; i++)
{
    g->vertex[i].mark = MARK_VERTEX_NONE;
    degree[i] =pst_img_graph_vertex_count_neighbours(g,i);
  }
  pst_img_graph_mark_vertex_removal(g,degree);
  
  for(i = 0; i < g->n; i++)
{
    if(g->vertex[i].mark == MARK_VERTEX_REMOVED)
{
//       fprintf(stderr,"Removing %ld ",i);
      pst_img_graph_vertex_remove(g,i,w_i,wmag,FALSE);
//       fprintf(stderr,"\n");
      
    } 
  }
  
  pst_img_graph_remove_paralel_edges(g);
  for(i = 0; i < g->n; i++)
{
    g->vertex[i].mark = MARK_VERTEX_NONE;
  }
  free(degree);
}


void pst_img_graph_check_consistency(pst_img_graph_t *g)
{
  long int i;
  bool_t test = TRUE;
  for(i = 0; i < g->n; i++)
{
    pst_vertex_data_t *v = &(g->vertex[i]);
    if(v->edge != oct_arc_NULL)
{
      if(v->id == -1)
{
        fprintf(stderr,"Dead-alive vertex [%ld](%ld) with edge not NULL\n",i,v->id);
        test = FALSE;
      }else{
        long int ind_edge = pst_img_graph_get_edge_num(v->edge);
        if(g->edge[ind_edge].edge == oct_arc_NULL)
{
          fprintf(stderr,"False link at vertex [%ld](%ld) pointing to [%ld]\n",i,v->id,ind_edge);
          test = FALSE;
        }else{
          long int orig_e = pst_img_graph_get_edge_origin(g,v->edge);
          long int orig_ed = pst_img_graph_get_edge_origin(g,g->edge[ind_edge].edge);
          long int orig_ed_s = pst_img_graph_get_edge_origin(g,oct_sym(g->edge[ind_edge].edge));
          if( (orig_e != orig_ed) && (orig_e != orig_ed_s))
{
            fprintf(stderr,"Inconsitency of org-dest at vertex [%ld](%ld) pointing to [%ld]\n",i,v->id,ind_edge);
            test = FALSE;
          }
        }
      }
    
      
    }
    
    
  }
  if( test)  fprintf(stderr,"Result - OK\n");
  else fprintf(stderr,"Result - FAIL\n");
}

void pst_img_graph_reorganize(pst_img_graph_t *g)
{

}

imgsys_t *pst_img_graph_build_integration_system(pst_img_graph_t *g,double *iW,long int NX_Z,long int NY_Z,long int* *ref_tab)
{
  
//   fprintf(stderr,"Oi\n");
  long int *ind_ix = (long int*)malloc(sizeof(long int)*g->n);
  imgsys_equation_t *eq = (imgsys_equation_t*)notnull(malloc((g->n)*sizeof(imgsys_equation_t)), "no mem");
  int i;
  int N = 0;
  for(i = 0; i < g->n ; i++)
{
     
    imgsys_equation_t *eqk = &(eq[N]);
    
    int nt = 0; /* Number of terms in equation. */
    eqk->rhs = 0.0;
    eqk->ix[nt] = i;  eqk->cf[nt] = 1.00; nt++;
    ind_ix[i] = N;
    pst_vertex_data_t *v = &(g->vertex[i]);
    
    
    
    
    int num_neighbours = (v->id == -1 ? 0:pst_img_graph_vertex_count_neighbours(g,i));
    oct_arc_t e0 = v->edge;
    if( num_neighbours > 0 )
{
      oct_arc_t e = e0;

      do{
        double d = pst_img_graph_get_edge_delta(g,e);
        double w = pst_img_graph_get_edge_weight(g,e);
        long int dest = pst_img_graph_get_edge_origin(g,oct_sym(e));
        assert(dest != -1);
        long int id = g->vertex[dest].id;
        assert(id != -1);
        eqk->ix[nt] = dest;
        assert(w > 0) ;
        eqk->cf[nt] = -w;
        eqk->rhs += -w*d;
        nt++;
        e = oct_onext(e);
      }while(e != e0);
      assert(nt <= MAXCOEFS);
      double wtot = 0;
      int j;
      for(j = 1; j < nt; j++)
{ assert(eqk->cf[j] < 0); wtot+= -eqk->cf[j]; } 
      assert(wtot > 0);
      iW[i] = wtot;
      for(j = 1; j < nt; j++)
{ eqk->cf[j]/=wtot;  }
      eqk->rhs /= wtot;
       eqk->nt = nt;
      N++;
    }else{
      ind_ix[i] = -1;
      iW[i] = 0;
    }
    
  }
  
   long int k;
   for (k = 0; k < N; k++)
   { imgsys_equation_t *eqk = &(eq[k]);
     int nt = eqk->nt;
     int mt = 0;
     int i;
     for(i = 0; i < nt; i++)
      { /* Get the temporay index {xyi}: */
        int xyi = eqk->ix[i];
        /* Get the definitive index {ki}: */
        int ki = ind_ix[xyi];
        if (ki >= 0)
        { /* Append the term to the equation: */
          int j = mt;
          eqk->ix[j] = ki;
          eqk->cf[j] = eqk->cf[i];
          assert(!isnan(eqk->cf[j]));
          mt++;
        }
      }
      eqk->nt = mt;
  }

  /* Build the inverse tables {col[0..N-1],row[0..N-1]}: */
   int *col = (int *)notnull(malloc(N*sizeof(int)), "no mem");
   int *row = (int *)notnull(malloc(N*sizeof(int)), "no mem");
   long int count_idx = 0;
   for(i =0; i < g->n; i++)
{
     if(ind_ix[i] >= 0)
{
       long int x,y;
       pst_img_graph_get_vertex_image_indices(&(g->vertex[i].coords),NX_Z,NY_Z,&x,&y);
       col[count_idx] = x;
       row[count_idx] = y;
       count_idx++;
     }
   }
   assert(count_idx == N);
   
//   long int *ix = (long int *)notnull(malloc((NXY_Z)*sizeof(long int)), "no mem");
//   { 
//     long int xy;
//     for(xy = 0; xy < NXY_Z; xy++)
{
//       ix[xy] = -1;
//     }
//     long int count_idx = 0;
//     for (xy = 0; xy <g->n; xy++) 
//     { 
//       if( ind_ix[xy] >= 0)
//       {
//      long int k = ind_ix[xy];
//      long int x,y;
//      pst_img_graph_get_vertex_image_indices(&(g->vertex[xy]),NX_Z,NY_Z,&x,&y);
//      
//      col[count_idx] = x;
//      row[count_idx] = y;
//      ix[g->vertex[xy].id] =k;
//      count_idx++;
//       }
//     }
//     assert(count_idx == N);
//     
//   }
  
    /* Now package the equations as a system: */
//     imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, ix, col, row);
    *ref_tab = ind_ix;
    imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, NULL, col,row);  
    return S;
}

void pst_img_graph_solve_system(
    pst_img_graph_t *g,
    imgsys_t *S,
    double *iZ, 
    double *iW,
    long int *ref_tab,
    long int maxIter, 
    double convTol, 
    int para, 
    int szero, 
    bool_t verbose 
  )
  {
   
   if( g->n_valid == S->N )
{
     fprintf(stderr,"IGUAL\n");
   }else{ fprintf(stderr,"NAO IGUAL\n"); }
   
    double *Z = (double*)malloc(sizeof(double)*S->N);
    
    long int count_vt = 0;
    long int i;
    for(i = 0; i < g->n; i++)
{
      if(ref_tab[i] >= 0)
{
        Z[count_vt] = iZ[i];
        count_vt++;
      }
    }
    assert(count_vt == S->N);
    
     long int *queue = pst_img_graph_sort_equations(S,ref_tab,iW,g->n );
    pst_imgsys_solve(S, Z,queue, maxIter, convTol, para, szero, verbose, 0, NULL);
    count_vt = 0;
    for(i = 0; i < g->n; i++)
{
      if( ref_tab[i] >= 0)
{
        iZ[i] = Z[count_vt];
        count_vt++;
      }
    }
    assert(count_vt == S->N);
    free(Z);
//     free(queue);
  }
  
  
  
void pst_img_graph_copy_solution_from_shrunk(pst_img_graph_t *jg,double *jZ,pst_img_graph_t *ig, double *iZ)
{
    long int i,j;
  assert(ig->n >= jg->n);
  
  j = 0;
  for(i = 0; i < ig->n; i++)
{
    
    pst_vertex_data_t *vIG = &(ig->vertex[i]);
    if( vIG->id == -1) {
      iZ[i] = 0;
      continue;
    }
    pst_vertex_data_t *vJG = NULL;
    while(j < jg->n)
{
      vJG = &(jg->vertex[j]);
      if( vJG->id != -1) break;
      j++;
    }
    
    if( (j >= jg->n) || (vJG->id > vIG->id) )
{
      iZ[i] = 0; vIG->mark = MARK_VERTEX_REMOVED; 
    }else if ( vJG->id == vIG->id ) {
      iZ[i] = jZ[j]; vIG->mark = MARK_VERTEX_PRESERVED;
      j++;
    }else{
      demand(FALSE,"Vertices out of order!");
    }

  }
}

void pst_img_graph_estimate_from_shrunk(pst_img_graph_t *jg,double *jZ,pst_img_graph_t *ig, double *iZ)
{
  long int i;
  assert(ig->n >= jg->n);
  
  pst_img_graph_copy_solution_from_shrunk(jg,jZ,ig,iZ);
  
  for(i = 0; i < ig->n; i++)
{
    pst_vertex_data_t *vIG = &(ig->vertex[i]);
    if(vIG->id == -1) continue;
    if(vIG->mark == MARK_VERTEX_REMOVED)
{ 
      oct_arc_t e0 = vIG->edge;
      if(e0 == oct_arc_NULL) continue;
      oct_arc_t e =  e0;
      double sW = 0;
      double sWZ = 0;
      do{
        long int ind_dst = pst_img_graph_get_edge_origin(ig,oct_sym(e));
        assert(ig->vertex[ind_dst].mark == MARK_VERTEX_PRESERVED);
        double d = pst_img_graph_get_edge_delta(ig,e);
        double w = pst_img_graph_get_edge_weight(ig,e);
        sW+= w;
        sWZ+=w*(iZ[ind_dst]-d);
        e = oct_onext(e);
      }while(e0 != e);
      assert(sW > 0);
      iZ[i] = sWZ/sW;
    }
  }
  
    
}


pst_img_graph_t *pst_img_graph_copy(pst_img_graph_t *g)
{
  /*First, count only valids*/
  
  auto oct_arc_t linha( oct_arc_t ed); 

  
  
  long int i;
  long int *eq_vector_vt = (long int*)malloc(sizeof(long int)*(g->n));
  long int valid_n = 0;
  for(i = 0; i < g->n;i++)
{
    if(g->vertex[i].id != -1)
{
      eq_vector_vt[i] = valid_n;
      valid_n++;
    }else{
      eq_vector_vt[i] = -1;
    }
  }
  
  long int *eq_vector_ed = (long int*)malloc(sizeof(long int)*(g->m));
  oct_arc_t *ref_tab = (oct_arc_t*)malloc(sizeof(oct_arc_t)*(g->m)*2);
  long int valid_m  = 0;
  for(i = 0; i < g->m;i++)
{
    if(g->edge[i].edge != oct_arc_NULL)
{
      eq_vector_ed[i] = valid_m;
      valid_m++;
    }else{
      eq_vector_ed[i] = -1;
    }
  }
  
  pst_img_graph_t  *ng = pst_img_graph_create(valid_m,valid_n);
  /*Easy part, copy the raw data*/
  for(i = 0; i < g->n;i++)
{
    pst_vertex_data_t *v = &(g->vertex[i]);
    if(eq_vector_vt[i] != -1)
{
      pst_img_graph_vertex_add(ng,v->id,oct_arc_NULL, v->coords );

    }
  }
  
  
  for(i = 0; i < g->m;i++)
{
    oct_arc_t e = g->edge[i].edge;
    if(eq_vector_ed[i] != -1)
{
      long int ind_e = pst_img_graph_get_edge_num(e);
      
      long int org =  pst_img_graph_get_edge_origin(g,e);
      long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
       pst_path_t p = pst_img_graph_get_edge_path(g,e);
      double delta = pst_img_graph_get_edge_delta(g,e);
      double weight = pst_img_graph_get_edge_weight(g,e);
      
      char *label = g->edge[ind_e].data->label;
      char *nl = (label == NULL ? NULL : (char*)malloc(sizeof(char)*(strlen(label)+1)));
      if(label != NULL ) strcpy(nl,label);
      
       pst_path_t np = p;
       if(p.v != NULL ) {
        np.v = (r2_t*)malloc(sizeof(r2_t)*(p.n));
        memcpy(np.v,p.v,sizeof(r2_t)*(p.n));
       }
      /*We dont do much with it now...*/
      long int new_org = eq_vector_vt[org];
      long int new_dst = eq_vector_vt[dst];
      oct_arc_t new_e =  pst_img_graph_edge_insert(ng,new_org,new_dst,delta,weight,nl,np);
      long int ind_e_dir = pst_img_graph_get_dir_edge_num(e);
      long int ind_e_sym = pst_img_graph_get_dir_edge_num(oct_sym(e));
      ref_tab[ind_e_dir] = new_e;
      ref_tab[ind_e_sym] = oct_sym(new_e);
      
     
    }
  }
 
  oct_arc_t linha( oct_arc_t ed)
{
    long int ind =  pst_img_graph_get_dir_edge_num(ed);
    return ref_tab[ind];
  }
 
 
  /*Hard part, assemble the graph data*/
  for(i = 0; i < g->m;i++)
{
    oct_arc_t e = g->edge[i].edge;
    long int ind_el = eq_vector_ed[i];
    if( ind_el != -1)
{
        
        oct_splice(oct_oprev(linha(e)),linha(oct_oprev(e)));
        oct_splice(oct_oprev(oct_sym(linha(e))),linha(oct_oprev(oct_sym(e))));
    }
  }
  
//   fprintf(stderr,"Checking consistency...");
//   for(i = 0; i < g->m; i++)
{
//     oct_arc_t e = g->edge[i].edge;
//     if( e != oct_arc_NULL)
{
//       long int ind = pst_img_graph_get_dir_edge_num(e);
//       oct_arc_t e_eqv = ref_tab[ind];
//       
//       assert(oct_onext(e_eqv) == linha(oct_onext(e)) );
//       
//       assert(oct_onext(oct_sym(e_eqv)) == linha(oct_onext(oct_sym(e))) );    
//     
//     }
//   }
//   fprintf(stderr,"OK\n");
  
  
  free(eq_vector_vt);
  free(eq_vector_ed);
  free(ref_tab);
  return ng;
}



void pst_img_graph_write_vertex(FILE *wr, pst_vertex_data_t *v)
{
  long int ind_edge = (v->edge == oct_arc_NULL ? -1 : pst_img_graph_get_dir_edge_num(v->edge));
  fprintf(wr,"%09ld %d %9.6f %9.6f %09ld\n",v->id,v->mark,v->coords.c[0],v->coords.c[1],ind_edge);
}

void pst_img_graph_write_path(FILE *wr, pst_path_t p)
{
  fprintf(wr,"%ld %d ",p.n,p.reverse);
  long int i;
  for(i = 0; i < p.n; i++)
{
    fprintf(wr,"%9.6f %9.6f ",p.v[i].c[0],p.v[i].c[1]);
  }
}

void pst_img_graph_write_edge(FILE *wr, pst_img_graph_t *g, oct_arc_t e)
{
  if(e == oct_arc_NULL)
{ 
    fprintf(wr,"-1\n");
    return;
  }
  long int ind = pst_img_graph_get_edge_num(e);
  long int ind_dir = pst_img_graph_get_dir_edge_num(e);
  long int org = pst_img_graph_get_edge_origin(g,e);
  long int dst = pst_img_graph_get_edge_origin(g,oct_sym(e));
  double delta = pst_img_graph_get_edge_delta(g,e);
  double weight = pst_img_graph_get_edge_weight(g,e);
  char *label = g->edge[ind].data->label;
  pst_path_t p = pst_img_graph_get_edge_path(g,e);
  
  fprintf(wr,"%09ld %09ld %09ld %09ld %9.6e %9.6e ",ind,ind_dir,org,dst,delta,weight);


  oct_arc_t e_next = oct_onext(e);
  
  
  long int ind_e_next = pst_img_graph_get_edge_num(e_next);
  long int ind_dir_e_next = pst_img_graph_get_dir_edge_num(e_next);
  int long_bit_enext = ind_dir_e_next&1;
    
  oct_arc_t e_next_sym = oct_onext(oct_sym(e));
  
  long int ind_e_next_sym = pst_img_graph_get_edge_num(e_next_sym);
  long int ind_dir_e_next_sym = pst_img_graph_get_dir_edge_num(e_next_sym);
  int long_bit_enext_sym = ind_dir_e_next_sym&1;
  
  fprintf(wr,"%ld %d %ld %d ",ind_e_next,long_bit_enext,ind_e_next_sym,long_bit_enext_sym);
  
  
  if(label == NULL)
{
    fprintf(wr,"0 ");
  }else{
    fprintf(wr,"%d ",strlen(label));
    fprintf(wr,"%s ",label);
  }
  pst_img_graph_write_path(wr,p);
  fprintf(wr,"\n");
}


void pst_img_graph_write(FILE *wr, pst_img_graph_t *g)
{
  /*First, write static data*/
  fprintf(wr,"%ld %ld %ld %ld \n",g->n_valid, g->m_valid,g->n, g->m);
  long int i;
  for(i = 0; i < g->n; i++)
{
    pst_img_graph_write_vertex(wr,&(g->vertex[i]));
  }
  for(i = 0; i < g->m; i++)
{
    pst_img_graph_write_edge(wr,g,g->edge[i].edge);
  }
  fprintf(wr,"\n");
}

void pst_img_graph_read_vertex(FILE *wr, pst_img_graph_t  *g)
{
  long int id;
  int mark;
  r2_t coords;
  long int ind_edge;
  demand(
    fscanf(wr,"%ld %d %lf %lf %ld",&id,&mark,&(coords.c[0]),&(coords.c[1]),&ind_edge) == 5,
    "Cannot read vertex"
  );
  /*Who is responsible to fill the  real edge's value is the pst_img_read_edges*/
  long int ix = pst_img_graph_vertex_add(g,id,oct_arc_NULL,coords);
  g->vertex[ix].mark = mark;
  
}

pst_path_t pst_img_graph_read_path(FILE *wr)
{
  pst_path_t p;
  int reverse;
  demand(fscanf(wr,"%ld %d",&(p.n),&(reverse)) == 2, "Cannot read path");
  p.reverse = (reverse == 1);
  p.v = NULL;
  if(p.n == 0) return p;
  p.v = (r2_t*) malloc(sizeof(r2_t)*(p.n));
  long int i;
  for(i = 0; i < p.n; i++)
{
    demand(fscanf(wr,"%lf %lf",&(p.v[i].c[0]),&(p.v[i].c[1])) == 2, "Cannot read path elements");
  }
  return p;
}

void pst_img_graph_read_edge(FILE *wr, pst_img_graph_t *g, long int list_onext[], int list_long_bit[])
{
  
  long int ind;
  long int ind_dir;
  long int org ;
  long int dst;
  double delta;
  double weight;

  demand(fscanf(wr,"%ld",&ind) == 1, "Cannot read edge index");
  if(ind != -1)
{
  
    demand(
      fscanf(wr,"%ld %ld %ld %lf %lf",&ind_dir,&org,&dst,&delta,&weight) == 5,
      "Cannot read edge data"
      );
      
    demand(
      fscanf(wr,"%ld %d %ld %d",&(list_onext[0]),&(list_long_bit[0]),&(list_onext[1]),&(list_long_bit[1])) == 4,
      "Cannot read edge connectivity"
    );
    
    long int label_len;
    demand( fscanf(wr,"%ld",&(label_len)) == 1, "Cannot read edge label lenght");
    char *label = (label_len == 0  ? NULL : (char*)malloc(sizeof(char)*(label_len +1)));
    long int i;
    fgetc(wr); /*skip the first blank space*/
    for(i = 0; i < label_len; i++)
{
      fscanf(wr,"%c",&(label[i]));
    }
    if(label != NULL) label[i] = '\0';
    
    pst_path_t p = pst_img_graph_read_path(wr);
    pst_img_graph_edge_insert(g,org,dst,delta,weight,label,p);
  }
  else{
    g->edge[g->m].edge = oct_arc_NULL;
    g->m++;
  }
  
    
  
}



pst_img_graph_t *pst_img_graph_read(FILE *wr)
{
  long int n_valid,m_valid, n,m;
  demand( fscanf(wr,"%ld %ld %ld %ld",&n_valid,&m_valid,&n,&m) == 4,"Cannot read graph header properly");
  assert( (n_valid >= 0) && (m_valid >= 0) && (n >= 0) && (m >=0) );
  
  pst_img_graph_t *g = pst_img_graph_create(m,n);
  
  long int i;
  for(i = 0; i < n ; i++)
{
    pst_img_graph_read_vertex(wr,g);
  }
  long int *list_next = (long int*)malloc(sizeof(long int)*2*m);
  int *list_long_bit = (int*)malloc(sizeof(int)*2*m);
  for(i = 0; i < m ; i++)
{
    pst_img_graph_read_edge(wr,g,&(list_next[2*i]),&(list_long_bit[2*i]));
  }
  
  auto oct_arc_t recover_edge(long int index, int long_bit);
  oct_arc_t  recover_edge(long int index, int long_bit)
{
    assert((index >= 0 ) && (index < g->m));
    oct_arc_t e = g->edge[index].edge;
    assert(e != oct_arc_NULL);
    if(long_bit == 1) e = oct_sym(e);
    return e;
  }
  
  /*Hard part*/
  for(i = 0; i < g->m;i++)
{
    oct_arc_t e = g->edge[i].edge;
    if(e != oct_arc_NULL)
{
        
        oct_arc_t e_next = recover_edge(list_next[2*i],list_long_bit[2*i]);
        oct_arc_t e_next_sym = recover_edge(list_next[(2*i)+1],list_long_bit[(2*i)+1]);
        
        oct_splice( oct_oprev(e_next),e);
        oct_splice( oct_oprev(e_next_sym),oct_sym(e));
    }
  }
  
  free(list_next);
  free(list_long_bit);
  
  return g;
}

long int pst_img_graph_find_nearest_vertex(pst_img_graph_t *g, r2_t p)
{
 long int i;
 double min_dist_sqr = +INF;
 long int min_i = -1;
 for(i = 0; i < g->n; i++)
{
   pst_vertex_data_t *v = &(g->vertex[i]);
   if(v->id == -1) continue;
   double d= r2_dist_sqr(&p,&(v->coords));
   if(d < min_dist_sqr)
{ min_dist_sqr = d; min_i = i; }
 }
  
  return min_i ;
}

double pst_img_graph_compute_left_face_curl(pst_img_graph_t *g, oct_arc_t e0)
{
  oct_arc_t e = e0;
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
  long int count= 0;
  do{
    long int ind_org = pst_img_graph_get_edge_origin(g,e);
    assert(ind_org != -1);
    pst_vertex_data_t *v = &(g->vertex[ind_org]);
    assert(v->id != -1);
    r2_add(&p,&(v->coords),&p);
    count++;
    e = oct_lnext(e);
  }while(e != e0);
  r2_scale(1.0/(double)count,&p,&p);
  return p;
}

void pst_img_graph_put_curl_into_image(pst_img_graph_t *g,float_image_t *OZ)
{
  long int i;
  for(i = 0; i < g->m; i++)
{
    oct_arc_t e = g->edge[i].edge;
    if(e != oct_arc_NULL)
{
      int j;
      for(j = 0; j < 2; j++)
{
        double c = pst_img_graph_compute_left_face_curl(g,e);
        r2_t p = pst_img_graph_left_face_baricenter(g,e);
        long int x,y;
        pst_img_graph_get_vertex_image_indices(&p,OZ->sz[1],OZ->sz[2],&x,&y);
        float_image_set_sample(OZ,0,x,y,c);
        e = oct_sym(e);
      }
    }
  }
}

void pst_img_graph_put_solution_into_image(pst_img_graph_t *g, double *iZ, float_image_t *OZ)
{
  long int i;
  float_image_fill_channel(OZ,0,0.0);
  for(i = 0; i < g->n; i++)
{
    pst_vertex_data_t *v = &(g->vertex[i]);
    if(v->id == -1) continue;
    double z = iZ[i];
    long int x,y;
    pst_img_graph_get_vertex_image_indices(&(v->coords),OZ->sz[1],OZ->sz[2],&x,&y);
    float_image_set_sample(OZ,0,x,y,z);
  }
}

void pst_img_graph_put_error_into_image(pst_img_graph_t *g,double *iZ,double*iW, float_image_t *RZ,float_image_t *OZ)
{
  long int i;
  float_image_fill_channel(OZ,0,0.0);
  
  double sumWD= 0;
  double sumW = 1.0e-300;
  
  for(i = 0; i < g->n; i++)
{
    pst_vertex_data_t *v = &(g->vertex[i]);
    
    if(v->id == -1) continue;
    double iz = iZ[i];
    double iw = iW[i];
    if(iw == 0) continue;
    long int x,y;
    pst_img_graph_get_vertex_image_indices(&(v->coords),OZ->sz[1],OZ->sz[2],&x,&y);
    double rz = float_image_get_sample(RZ,0,x,y);
    double dz = iz - rz;
    sumW+= iw;
    sumWD+= iw*dz;
  }
  
  double avgdz = sumWD/sumW;
  
  for(i = 0; i < g->n; i++)
{
    pst_vertex_data_t *v = &(g->vertex[i]);
    if(v->id == -1) continue;
    double iw = iW[i];
    if(iw == 0) continue;
    double iz = iZ[i];
    long int x,y;
    pst_img_graph_get_vertex_image_indices(&(v->coords),OZ->sz[1],OZ->sz[2],&x,&y);
    double rz = float_image_get_sample(RZ,0,x,y);
    double dz = iz - rz;
    float_image_set_sample(OZ,0,x,y,dz - avgdz);
  }
  
}

void pst_img_graph_integration_recursive( 
  pst_img_graph_t *g,
  double *iZ,
  double *iW,
  double wmag,
  long int maxIter,
  double convTol, 
  int para, 
  int szero, 
  bool_t verbose,
  int level,
  float_image_t *OZ, /*debug only*/
  float_image_t *RZ,
  char *out_prefix,
  bool_t debug
)
{
  fprintf(stderr,"Starting level %d\n",level);
  
  if(debug)
{
    char *filename = NULL;
    pst_img_graph_put_curl_into_image(g,OZ);
    asprintf(&filename,"%s-%02d-CL.fni",out_prefix,level);
    FILE *wr_curl = open_write(filename,FALSE);
     float_image_write(wr_curl,OZ);
    fclose(wr_curl);
    free(filename);
  }
  
  if(debug)
{
    pst_img_graph_mark_vertex_removal(g,NULL);
    char *filename = NULL;
    asprintf(&filename,"%s-%02d-GR.txt",out_prefix,level);
    FILE *wr_grph = open_write(filename,FALSE);
     pst_img_graph_print(wr_grph,g);
    fclose(wr_grph);
    free(filename);
    filename = NULL;
    asprintf(&filename,"%s-%02d-GR.grf",out_prefix,level);
    wr_grph = open_write(filename,FALSE);
     pst_img_graph_write(wr_grph,g);
    fclose(wr_grph);
    free(filename);
  }
  
  if( g->n_valid > 2)
{
    fprintf(stderr,"Reducing Graph level[%d] with [%ld] vertices [%ld] valid\n",level,g->n,g->n_valid);
    fprintf(stderr,"Copying...");
    pst_img_graph_t *jg = pst_img_graph_copy(g);
//     fprintf(stderr,"OK.\nShrinking...");
    pst_img_graph_shrink(jg,NULL,wmag);
//     fprintf(stderr,"OK.\nNext\n");

    double *jZ = (double*)notnull(malloc(sizeof(double)*jg->n),"bug");
    double *jW = (double*)notnull(malloc(sizeof(double)*jg->n),"bug");
    
//     pst_img_graph_integration_recursive(jg,jZ,jW,wmag,(int)ceil(maxIter*M_SQRT2),convTol/M_SQRT2,para,szero,verbose,level+1,OZ,RZ,out_prefix,debug);
    double ratio = (((double)g->n_valid)/((double)jg->n_valid));
    double newConvTol = convTol/sqrt(ratio);
    int newmaxIter = ceil(sqrt(ratio)*maxIter);
    fprintf(stdout,"%9.6lf\n",ratio);
    fprintf(stderr,"Ratio %9.6lf G_valid %ld JG_valid %ld New ConvTol %lf NewmaxIter %d \n",ratio,g->n_valid,jg->n_valid, newConvTol,newmaxIter);
    pst_img_graph_integration_recursive(jg,jZ,jW,wmag,newmaxIter,newConvTol,para,szero,verbose,level+1,OZ,RZ,out_prefix,debug);
    pst_img_graph_estimate_from_shrunk(jg,jZ,g,iZ);
    /* The vertices must be painted now*/
    free(jZ);
    free(jW);
  }else{
    fprintf(stderr,"End of Recursion - returning\n");
    long int i;
    for(i = 0; i < g->n;i++) iZ[i] = 0; 
  }
  
  
  
  
  fprintf(stderr,"Solving Graph level[%d] with [%ld] vertices [%ld] valid\n",level,g->n,g->n_valid);

  long int *ref_tab;
  imgsys_t *S = pst_img_graph_build_integration_system(g,iW,OZ->sz[1],OZ->sz[2],&ref_tab);
  
  if(debug)
{
  
    char *filename = NULL;
    asprintf(&filename,"%s-%02d-S.txt",out_prefix,level);
    FILE *wr_dump = open_write(filename,FALSE);
    pst_imgsys_write(wr_dump,S);
    fclose(wr_dump);
    free(filename);
    
    pst_img_graph_put_solution_into_image(g,iZ,OZ);
    filename = NULL;
    asprintf(&filename,"%s-%02d-iZ.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
    
    pst_img_graph_put_solution_into_image(g,iW,OZ);
    filename = NULL;
    asprintf(&filename,"%s-%02d-iW.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
    
    if(RZ != NULL)
{
      pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
      filename = NULL;
      asprintf(&filename,"%s-%02d-iE.fni",out_prefix,level);
      wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }
  }
  pst_img_graph_solve_system(g,S,iZ,iW,ref_tab,maxIter,convTol,para,szero,verbose);
  

  if(debug)
{
    
    if(RZ != NULL)
{
      pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
      char *filename = NULL;
      asprintf(&filename,"%s-%02d-oE.fni",out_prefix,level);
      FILE *wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }

  }
  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  if(debug)
{
    char *filename = NULL;
    asprintf(&filename,"%s-%02d-oZ.fni",out_prefix,level);
    FILE *wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
    free(filename);
  }
  free(ref_tab);
  pst_imgsys_free(S);
  
}


void pst_img_graph_integration( 
  pst_img_graph_t *g,
  double *iZ,
  double *iW,
  long int maxIter,
  double convTol, 
  int para, 
  int szero, 
  bool_t verbose,
  int level,
  float_image_t *OZ, /*debug only*/
  float_image_t *RZ,
  char *out_prefix
)
{
  
  char *filename_grph = NULL;
  asprintf(&filename_grph,"%s-Graph-%02d.txt",out_prefix,level);
  FILE *wr_grph = open_write(filename_grph,FALSE);
   pst_img_graph_print(wr_grph,g);
  fclose(wr_grph);
  long int i;
  
  for(i = 0; i < g->n;i++) iZ[i] = 0; 
  fprintf(stderr,"Solving Graph with [%ld] vertices [%ld] valid\n",g->n,g->n_valid);

  long int *ref_tab;
  imgsys_t *S = pst_img_graph_build_integration_system(g,iW,OZ->sz[1],OZ->sz[2],&ref_tab);
  
  
  char *filename = NULL;
  asprintf(&filename,"%s-%02d-S.txt",out_prefix,level);
  FILE *wr_dump = open_write(filename,FALSE);
  pst_imgsys_write(wr_dump,S);
  fclose(wr_dump);
  
  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  filename = NULL;
  asprintf(&filename,"%s-%02d-iZ.fni",out_prefix,level);
  wr_dump = open_write(filename,FALSE);
  float_image_write(wr_dump,OZ);
  fclose(wr_dump);
  
  if(RZ != NULL)
{
    pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
    filename = NULL;
    asprintf(&filename,"%s-%02d-iE.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
  }
  
  pst_img_graph_solve_system(g,S,iZ,iW,ref_tab,maxIter,convTol,para,szero,verbose);
  

  
  if(RZ != NULL)
{
    pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
    filename = NULL;
    asprintf(&filename,"%s-%02d-oE.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
  }

  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  filename = NULL;
  asprintf(&filename,"%s-%02d-oZ.fni",out_prefix,level);
  wr_dump = open_write(filename,FALSE);
  float_image_write(wr_dump,OZ);
  fclose(wr_dump);

  free(ref_tab);
  pst_imgsys_free(S);
  
}

void pst_img_graph_free(pst_img_graph_t *g)
{
  long int i;
  for(i = 0; i < g->m; i++)
{
    oct_arc_t e = g->edge[i].edge;
    if(e != oct_arc_NULL) pst_img_graph_edge_remove(g,e);
    if(g->edge[i].data != NULL)
{
      pst_edge_data_free(g->edge[i].data);
      
    }
  }
  free(g->edge);
  free(g->vertex);
  
  free(g);
}


long int *pst_img_graph_sort_equations(imgsys_t *S,long int *ref_tab, double *iW ,long int g_n )
{
 
  #define MAXDEGREE (2*((MAXCOEFS)-1))
  int N = S->N;
  assert(MAXCOEFS <= 20);
  /* Build the graph of the equations */
  long int *graph = (long int *)notnull(malloc(MAXDEGREE*N*sizeof(long int)), "no mem");
  /* Value of out_degree[k] is the number of neighbours of pixel[k] with weight less than OW[k]  */
  long int *out_degree =  (long int *)notnull(malloc(N*sizeof(long int)), "no mem");
  /* Value of in_degree[k] is the number of unprocessed neighbours of pixel[k] with weight greater than OW[k]  */
  long int *in_degree =  (long int *)notnull(malloc(N*sizeof(long int)), "no mem");
  
  
  long int *inv_ref_tab = (long int*)malloc(sizeof(long int)*S->N);
  long int count_vt = 0;
  int k;
  for(k = 0; k < g_n; k++)
{
    if(ref_tab[k] != -1)
{
      inv_ref_tab[count_vt] = k;
      count_vt++;
    }
  }
  assert(count_vt == S->N);
  
  auto bool_t arrow( long int k1, long int k2);
  bool_t arrow( long int k1, long int k2)
{
    long int refk1 = inv_ref_tab[k1];
    long int refk2 = inv_ref_tab[k2];
    double w1 = iW[refk1];
    double w2 = iW[refk2];
    return w1 < w2;
  }
  
  
  
  for(k = 0; k < N; k++)
{ in_degree[k] = 0; out_degree[k] = 0; }
  for(k = 0; k < N; k++)
{
    imgsys_equation_t *eqk = &(S->eq[k]);
    assert( eqk->ix[0] == k);
    int j;
    for(j = 1; j < eqk->nt; j++)
    {
      long int i = eqk->ix[j];
      if(arrow(k,i))
{ in_degree[i]++; graph[MAXDEGREE*k + out_degree[k]] = i; out_degree[k]++;}
      if(arrow(i,k))
{ in_degree[k]++; graph[MAXDEGREE*i + out_degree[i]] = k; out_degree[i]++;}
    }
    
  }
  
  
    
  long int queue_free = 0;
  long int queue_start = 0;
  long int *queue = (long int *)notnull(malloc(N*sizeof(long int)), "no mem");
  
  auto void eq_queue_insert(long int index);
  void eq_queue_insert(long int index)
{
    assert(queue_free < N);
    queue[queue_free] = index;
    queue_free++;
  }

  auto long int eq_queue_remove( void );
  long int eq_queue_remove( void )
{
    assert(queue_start < queue_free );
    long int i = queue[queue_start];
    queue_start++;
    return i;
    
  }  
  
  for(k = 0; k < N; k++)
  {if( in_degree[k] == 0)
{ eq_queue_insert(k); }}
  
  while(queue_start < queue_free )
  {
    long int k = eq_queue_remove();
    int j;
    for(j = 0; j < out_degree[k]; j++)
    {
      long int i = graph[MAXDEGREE*k + j];
      assert(in_degree[i] > 0);
      in_degree[i]--;
      if(in_degree[i] == 0)
{ eq_queue_insert(i); }
    }
  }
  
  assert(queue_free == N);
  
  free(graph);
  free(in_degree);
  free(out_degree);
  free(inv_ref_tab);
  
  return queue;
 
}
