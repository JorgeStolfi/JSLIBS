/* See {pst_img_graph.h} */
/* Last edited on 2024-12-24 18:59:41 by stolfi */
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

pst_imgsys_t *pst_img_graph_build_integration_system
  ( pst_img_graph_t *g,
    double *iW,
    int32_t NX_Z,
    int32_t NY_Z,
    int32_t **ref_tab
  )
  { /* fprintf(stderr,"Oi\n"); */
    int32_t *ind_ix = talloc(g->NV, int32_t);
    pst_imgsys_equation_t *eq = talloc(g->NV, pst_imgsys_equation_t);
    uint32_t N = 0;
    for (uint32_t i = 0; i < g->NV ; i++)
      { pst_imgsys_equation_t *eqk = &(eq[N]);
        uint32_t nt = 0; /* Number of terms in equation. */
        eqk->rhs = 0.0;
        eqk->ix[nt] = i;  eqk->cf[nt] = 1.00; nt++;
        ind_ix[i] = (int32_t)N;
        pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        uint32_t num_neighbours = (v->id == -1 ? 0:pst_img_graph_vertex_count_neighbours(g,i));
        haf_arc_t e0 = v->edge;
        if (num_neighbours > 0)
          { haf_arc_t e = e0;
            do
              { double d = pst_img_graph_get_edge_delta(g,e);
                double w = pst_img_graph_get_edge_weight(g,e);
                uint32_t dest = pst_img_graph_get_edge_origin(g,haf_sym(e));
                assert(dest != -1);
                uint32_t id = g->vertex[dest].id;
                assert(id != -1);
                eqk->ix[nt] = dest;
                assert(w > 0) ;
                eqk->cf[nt] = -w;
                eqk->rhs += -w*d;
                nt++;
                e = haf_onext(e);
              } while(e != e0);
            assert(nt <= MAXCOEFS);
            double wtot = 0;
            for (uint32_t j = 1; j < nt; j++)
              { assert(eqk->cf[j] < 0); wtot+= -eqk->cf[j]; } 
            assert(wtot > 0);
            iW[i] = wtot;
            for (uint32_t j = 1; j < nt; j++)
              { eqk->cf[j]/=wtot;  }
            eqk->rhs /= wtot;
            eqk->nt = nt;
            N++;
         }
       else
         { ind_ix[i] = -1;
           iW[i] = 0;
         }
      }
    for (uint32_t k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(eq[k]);
        uint32_t nt = eqk->nt;
        uint32_t mt = 0;
        for (uint32_t i = 0; i < nt; i++)
          { /* Get the temporay index {xyi}: */
            uint32_t xyi = eqk->ix[i];
            /* Get the definitive index {ki}: */
            int32_t ki = ind_ix[xyi];
            if (ki >= 0)
              { /* Append the term to the equation: */
                int32_t j = (int32_t)mt;
                eqk->ix[j] = (uint32_t)ki;
                eqk->cf[j] = eqk->cf[i];
                assert(!isnan(eqk->cf[j]));
                mt++;
              }
          }
        eqk->nt = mt;
      }

    /* Build the inverse tables {col[0..N-1],row[0..N-1]}: */
    uint32_t *col = talloc(N, uint32_t);
    uint32_t *row = talloc(N, uint32_t);
    uint32_t count_idx = 0;
    for (uint32_t i =0; i < g->NV; i++)
      { if (ind_ix[i] >= 0)
          { int32_t x, y;
            pst_img_graph_get_vertex_image_indices(&(g->vertex[i].coords), NX_Z, NY_Z, &x, &y);
            col[count_idx] = (uint32_t)x;
            row[count_idx] = (uint32_t)y;
            count_idx++;
          }
     }
     assert(count_idx == N);
   
    /* Now package the equations as a system: */
    /* pst_imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, ix, col, row); */
    *ref_tab = ind_ix;
    pst_imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, NULL, col,row);  
    return S;
  }

void pst_img_graph_solve_system
  ( pst_img_graph_t *g,
    pst_imgsys_t *S,
    double *iZ, 
    double *iW,
    int32_t *ref_tab,
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose 
  )
  { if ( g->NV_valid == S->N )
      { fprintf(stderr,"IGUAL\n"); }
    else
      { fprintf(stderr,"NAO IGUAL\n"); }
   
    double *Z = rn_alloc(S->N);
    
    int32_t count_vt = 0;
    for (uint32_t i = 0; i < g->NV; i++)
      { if (ref_tab[i] >= 0)
          { Z[count_vt] = iZ[i];
            count_vt++;
          }
       }
    assert(count_vt == S->N);
    
    uint32_t *queue = pst_img_graph_sort_equations(S,ref_tab,iW,g->NV );
    pst_imgsys_solve(S, Z, queue, maxIter, convTol, para, szero, verbose, 0, NULL);
    count_vt = 0;
    for (uint32_t i = 0; i < g->NV; i++)
      { if ( ref_tab[i] >= 0)
          { iZ[i] = Z[count_vt];
            count_vt++;
          }
       }
    assert(count_vt == S->N);
    free(Z);
    /* free(queue); */
  }

void pst_img_graph_put_solution_into_image(pst_img_graph_t *g, double *iZ, float_image_t *OZ)
  {
    float_image_fill_channel(OZ, 0, 0.0);
    for (uint32_t i = 0; i < g->NV; i++)
      { pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (v->id == -1) continue;
        double z = iZ[i];
        int32_t x,y;
        pst_img_graph_get_vertex_image_indices(&(v->coords), (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &x, &y);
        float_image_set_sample(OZ, 0, x, y, (float)z);
      }
  }

void pst_img_graph_put_error_into_image(pst_img_graph_t *g, double *iZ, double*iW, float_image_t *RZ, float_image_t *OZ)
  {
    float_image_fill_channel(OZ,0,0.0);

    double sumWD= 0;
    double sumW = 1.0e-300;

    for (uint32_t i = 0; i < g->NV; i++)
      { pst_img_graph_vertex_data_t *v = &(g->vertex[i]);

        if (v->id == -1) continue;
        double iz = iZ[i];
        double iw = iW[i];
        if (iw == 0) continue;
        int32_t x,y;
        pst_img_graph_get_vertex_image_indices(&(v->coords), (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &x, &y);
        double rz = float_image_get_sample(RZ, 0, x, y);
        double dz = iz - rz;
        sumW+= iw;
        sumWD+= iw*dz;
      }

    double avgdz = sumWD/sumW;

    for (uint32_t i = 0; i < g->NV; i++)
      {
        pst_img_graph_vertex_data_t *v = &(g->vertex[i]);
        if (v->id == -1) continue;
        double iw = iW[i];
        if (iw == 0) continue;
        double iz = iZ[i];
        int32_t x,y;
        pst_img_graph_get_vertex_image_indices(&(v->coords), (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &x,&y);
        double rz = float_image_get_sample(RZ,0,x,y);
        double dz = iz - rz;
        float_image_set_sample(OZ, 0, x, y, (float)(dz - avgdz));
      }
  }

void pst_img_graph_integration( 
  pst_img_graph_t *g,
  double *iZ,
  double *iW,
  uint32_t maxIter,
  double convTol, 
  bool_t para, 
  bool_t szero, 
  bool_t verbose,
  uint32_t level,
  float_image_t *OZ, /*debug only*/
  float_image_t *RZ,
  char *out_prefix
)
{ char *filename_grph = jsprintf("%s-Graph-%02d.txt",out_prefix,level);
  FILE *wr_grph = open_write(filename_grph,FALSE);
   pst_img_graph_print(wr_grph,g);
  fclose(wr_grph);
  for (uint32_t i = 0; i < g->NV;i++) iZ[i] = 0; 
  fprintf(stderr,"Solving Graph with [%d] vertices [%d] valid\n", g->NV, g->NV_valid);

  int32_t *ref_tab;
  pst_imgsys_t *S = pst_img_graph_build_integration_system(g, iW, (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &ref_tab);
  
  
  char *filename = jsprintf("%s-%02d-S.txt",out_prefix,level);
  FILE *wr_dump = open_write(filename,FALSE);
  free(filename);
  pst_imgsys_write(wr_dump,S);
  fclose(wr_dump);
  
  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  filename = jsprintf("%s-%02d-iZ.fni",out_prefix,level);
  wr_dump = open_write(filename,FALSE);
  free(filename);
  float_image_write(wr_dump,OZ);
  fclose(wr_dump);
  
  if (RZ != NULL)
{
    pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
    filename = jsprintf("%s-%02d-iE.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    free(filename);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
  }
  
  pst_img_graph_solve_system(g,S,iZ,iW,ref_tab,maxIter,convTol,para,szero,verbose);
  

  
  if (RZ != NULL)
{
    pst_img_graph_put_error_into_image(g,iZ,iW,RZ,OZ);
    filename = jsprintf("%s-%02d-oE.fni",out_prefix,level);
    wr_dump = open_write(filename,FALSE);
    free(filename);
    float_image_write(wr_dump,OZ);
    fclose(wr_dump);
  }

  pst_img_graph_put_solution_into_image(g,iZ,OZ);
  filename = jsprintf("%s-%02d-oZ.fni",out_prefix,level);
  wr_dump = open_write(filename,FALSE);
  free(filename);
  float_image_write(wr_dump,OZ);
  fclose(wr_dump);

  free(ref_tab);
  pst_imgsys_free(S);
  
}

uint32_t *pst_img_graph_sort_equations
  ( pst_imgsys_t *S,
    int32_t *ref_tab,
    double *iW ,
    uint32_t NV 
  )
{
 
  #define MAXDEGREE (2*((MAXCOEFS)-1))
  uint32_t N = S->N;
  assert(MAXCOEFS <= 20);
  /* Build the graph of the equations */
  uint32_t *graph = talloc(MAXDEGREE*N, uint32_t);
  /* Value of out_degree[k] is the number of neighbours of pixel[k] with weight less than OW[k]  */
  uint32_t *out_degree =  talloc(N, uint32_t);
  /* Value of in_degree[k] is the number of unprocessed neighbours of pixel[k] with weight greater than OW[k]  */
  uint32_t *in_degree =   talloc(N, uint32_t);
  
  
  int32_t *inv_ref_tab = (int32_t*)malloc(sizeof(int32_t)*S->N);
  int32_t count_vt = 0;
  for (int32_t k = 0; k < NV; k++)
{
    if (ref_tab[k] != -1)
{
      inv_ref_tab[count_vt] = k;
      count_vt++;
    }
  }
  assert(count_vt == S->N);
  
  auto bool_t arrow( uint32_t k1, uint32_t k2);
  bool_t arrow( uint32_t k1, uint32_t k2)
{
    int32_t refk1 = inv_ref_tab[k1];
    int32_t refk2 = inv_ref_tab[k2];
    double w1 = iW[refk1];
    double w2 = iW[refk2];
    return w1 < w2;
  }
  
  
  
  for (uint32_t k = 0; k < N; k++)
{ in_degree[k] = 0; out_degree[k] = 0; }
  for (uint32_t k = 0; k < N; k++)
{
    pst_imgsys_equation_t *eqk = &(S->eq[k]);
    assert( eqk->ix[0] == k);
    for (uint32_t j = 1; j < eqk->nt; j++)
    {
      uint32_t i = eqk->ix[j];
      if (arrow(k,i))
{ in_degree[i]++; graph[MAXDEGREE*k + out_degree[k]] = i; out_degree[k]++;}
      if (arrow(i,k))
{ in_degree[k]++; graph[MAXDEGREE*i + out_degree[i]] = k; out_degree[i]++;}
    }
    
  }
  
  
    
  uint32_t queue_free = 0;
  uint32_t queue_start = 0;
  uint32_t *queue = talloc(N, uint32_t);
  
  auto void eq_queue_insert(uint32_t index);
  void eq_queue_insert(uint32_t index)
{
    assert(queue_free < N);
    queue[queue_free] = index;
    queue_free++;
  }

  auto uint32_t eq_queue_remove( void );
  uint32_t eq_queue_remove( void )
{
    assert(queue_start < queue_free );
    uint32_t i = queue[queue_start];
    queue_start++;
    return i;
    
  }  
  
  for (uint32_t k = 0; k < N; k++)
  {if ( in_degree[k] == 0)
{ eq_queue_insert(k); }}
  
  while(queue_start < queue_free )
  {
    uint32_t k = eq_queue_remove();
    for (uint32_t j = 0; j < out_degree[k]; j++)
    {
      uint32_t i = graph[MAXDEGREE*k + j];
      assert(in_degree[i] > 0);
      in_degree[i]--;
      if (in_degree[i] == 0)
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
