/* See {pst_gr_from_slope_map.h} */
/* Last edited on 2025-03-15 14:12:55 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <jsfile.h>
#include <affirm.h>

#include <pst_slope_map.h>

#include <pst_gr.h>
#include <pst_gr_path.h>

#include <pst_gr_from_slope_map.h>

void pst_gr_get_diagonal_edge_data_from_maps
  ( float_image_t* IG,
    int32_t x, int32_t y,
    int32_t dirx,int32_t diry,
    double *d_P, double *w_P
  );
  
#define NONE pst_gr_NONE

/* IMPLEMENTATIONS */

pst_gr_t* pst_gr_from_slope_map(float_image_t* IG, bool_t add_diags)
  { 
    demand(! add_diags, "{add_diags} not implemented yet");
    demand(IG != NULL, "gradient map must not be {NULL}");
    
    int32_t NC_G, NX_G, NY_G; /* Cols and rows of {IG} and {IW}. */
    float_image_get_size(IG, &NC_G, &NX_G, &NY_G);
    demand((NC_G == 2) || (NC_G == 3), "gradient map must have 2 or 3 channels");

    /* Cols and rows of height map: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;

    /* The vertices are the height map pixels: */
    int32_t NV = NX_Z*NY_Z;
    
    /* The edges connect horz or vert adjacent adjacent height map pixels: */
    int32_t NE = NX_Z*NY_G + NX_G*NY_Z;
    
    pst_gr_t* gr = pst_gr_new((uint32_t)NV, (uint32_t)NE, (uint32_t)NX_Z, (uint32_t)NY_Z);

    auto void add_edge(int32_t orgx, int32_t orgy, int32_t dstx, int32_t dsty);
      /* Adds to the graph {gr} an edge betewwn vertex {org=[orgx,orgy]}
        and vertex {dst=[dstx,dsty]}. where {orgx,dstx} are in
        {0..NX_Z-1} and {orgy,dsty} are in {0..NY_Z-1}. */

    auto void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d_P, double *w_P);
      /* Sets {*d_P} to the height difference of the graph edge from
        vertex {u=[ux,uy]} to vertex {v=[vx,y]}, as estimated from the
        gradient map {IG}. Also sets {*w_P]} to the corresponding weight,
        estimated from the weight map {IW}; if {IW} is {NULL}, {*w_P]} is
        set to 1. Note that {u} and {v} correspond to pixels of the
        height map, hence to corners of the gradient and weight map
        grids. */

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { r2_t coords = (r2_t){{ x, y }};
            uint32_t vi = pst_gr_add_vertex(gr, x, y, coords);
            assert(vi == pst_gr_get_vertex_from_map_indices(gr, x, y));
          }
      }
                
    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { if (y > 0)
              { /* Add vert edge from {x,y} to {x,y-1}: */
                add_edge(x,y, x,y-1);
              }
            if (x > 0)
              { /* Add horz edge from {x,y} to {x-1,y}: */
                add_edge(x,y, x-1,y);
              }
          }
      }

    return gr;

    void add_edge(int32_t orgx, int32_t orgy, int32_t dstx, int32_t dsty)
      { pst_gr_vertex_t org = pst_gr_get_vertex_from_map_indices(gr, orgx,orgy);
        demand((org != NONE) && (org >= 0) && (org < NV), "invalid origin {orgx,orgy}");
        pst_gr_vertex_t dst = pst_gr_get_vertex_from_map_indices(gr, dstx,dsty);
        demand((dst != NONE) && (dst >= 0) && (dst < NV), "invalid destination {dstx,dsty}");
        double we, de;
        get_edge_parameters(orgx,orgy, dstx,dsty, &de, &we);
        if (we > 0)
          { pst_gr_arc_t ai = pst_gr_add_edge(gr, (uint32_t)org,(uint32_t)dst, de,we, pst_gr_path_NULL, TRUE);
            assert(ai != NONE);
          }
      }

    void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d_P, double *w_P)
      { pst_slope_map_get_edge_data(IG, xu, yu, xv-xu, yv-yu, d_P, w_P); }
      
  }
