/* See {tosl_mesh_make_ico.h} */
/* Last edited on 2024-11-20 05:21:47 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_mesh_make_ico.h>

tosl_mesh_t *tosl_mesh_make_ico(double R)
  {
    /* Compute the two basic (quantized, even) coordinates {Rq, Sq}: */
    R = fmax(R, 10.0); 
    double golden = (sqrt(5) - 1.0)/2.0;
    double S = golden*R;     
    tosl_coord_t Rq = (tosl_coord_t)floor(R + 0.5); Rq = 2*(Rq/2);
    tosl_coord_t Sq = (tosl_coord_t)floor(S + 0.5); Sq = 2*(Sq/2);
    
    /* Allocate the mesh structure: */
    int32_t NV = 12;
    int32_t NE = 30;
    tosl_mesh_t *mesh = tosl_mesh_new(NE, NV);

    fprintf(stderr, "  creating the %d vertices...\n", NV);
    for (int32_t i = 0; i < 3; i++)
      { int32_t j = (i + 1) % 3;
        for (int32_t dj = -1; dj <= +1; dj += 2)
          for (int32_t di = -1; di <= +1; di += 2)
            { tosl_point_t v = (tosl_point_t){{ 0, 0, 0 }};
              v.c[i] = di*Rq;
              v.c[j] = dj*Sq;
              char *lab = jsprintf("v%d.%d%d", i, (di+1)/2, (dj+1)/2); 
              tosl_mesh_add_vert(&v, lab, mesh);
            }
      }
    assert(mesh->NV == NV);
        
    fprintf(stderr, "  adding the 30 edges and linking the 12 axial triangles...\n");
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t r = 0; r <= 1; r++)
          { /* Axial edge: */
            tosl_vert_id_t kv0 = 4*i + r;
            tosl_vert_id_t kv1 = kv0 + 2;
            tosl_arc_id_t ja = tosl_mesh_add_edge(kv0, kv1, 0, mesh);
            /* Diagonal edges: */
            tosl_vert_id_t kv2 = (kv0 + r + 8) % 12;
            tosl_vert_id_t kv3 = kv2 + 1;
            tosl_arc_id_t ja02 = tosl_mesh_add_edge(kv0, kv2, 0, mesh); 
            tosl_arc_id_t ja03 = tosl_mesh_add_edge(kv0, kv3, 0, mesh);
            tosl_arc_id_t ja12 = tosl_mesh_add_edge(kv1, kv2, 0, mesh);
            tosl_arc_id_t ja13 = tosl_mesh_add_edge(kv1, kv3, 0, mesh);
            /* Triangles with the axial edge: */
            tosl_mesh_link_triang((r+1)%2, ja, ja12, tosl_sym(ja02), mesh);
            tosl_mesh_link_triang((r+0)%2, ja, ja13, tosl_sym(ja03), mesh);
          }
      }
    assert(mesh->NE == NE);
    
    /* Debug: */
    /* if (debug) { tosl_mesh_print(stderr, mesh); } */
    
    fprintf(stderr, "  linking the octant triangles...\n");
    tosl_mesh_link_triang(1,  2, 42, 22, mesh);
    tosl_mesh_link_triang(0, 12, 46, 24, mesh);
    tosl_mesh_link_triang(0,  6, 44, 32, mesh);
    tosl_mesh_link_triang(1, 16, 48, 34, mesh);
    tosl_mesh_link_triang(0,  4, 52, 26, mesh);
    tosl_mesh_link_triang(1, 14, 56, 28, mesh);
    tosl_mesh_link_triang(1,  8, 54, 36, mesh);
    tosl_mesh_link_triang(0, 18, 58, 38, mesh);

    return mesh;                   
  }
