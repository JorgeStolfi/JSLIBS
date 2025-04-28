/* See {psgr_from_slope_map.h} */
/* Last edited on 2025-04-27 03:22:41 by stolfi */
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
#include <i2.h>

#include <pst_slope_map.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_path.h>

#include <psgr_from_slope_map.h>

void psgr_get_diagonal_edge_data_from_maps
  ( float_image_t* IG,
    int32_t x, int32_t y,
    int32_t dirx,int32_t diry,
    double *d_P, double *w_P
  );
  
#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE
  /* Shorter names. */  

/* IMPLEMENTATIONS */

psgr_t* psgr_from_slope_map(float_image_t* IG, bool_t extrapolate, bool_t add_diags)
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
    
    psgr_t* gr = psgr_new((uint32_t)NV, (uint32_t)NE, (uint32_t)NX_Z, (uint32_t)NY_Z);

    auto void add_edge(i2_t *ipo, i2_t *ipd);
      /* Adds to the graph {gr} an edge from vertex {org} with pixel indices {ipo}
        to vertex {dst} with pixel indices {ipd}.  The pixel indices must be in
        {{0..NX_Z-1}..{0..NY_Z-1}}. */

    auto void get_edge_parameters(i2_t *ipu, i2_t *ipv,  double *d_P, double *w_P);
      /* Sets {*d_P} to the height difference of the graph edge from
        vertex {u} of {gr} with pixel indices {ipu} to vertex {v} with
        pixel indices {ipv}, as estimated from the gradient map {IG}.
        Also sets {*w_P]} to the corresponding weight, estimated from
        the weight map {IW}; if {IW} is {NULL}, {*w_P]} is set to 1.
        Note that {u} and {v} correspond to pixels of the height map,
        hence to corners of the gradient and weight map grids. */

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { i2_t ip = (i2_t){{ x, y }};
            psgr_vertex_t vi = psgr_add_vertex(gr, &ip);
            assert(vi.ix == psgr_vertex_from_pixel(gr, &ip).ix);
          }
      }
                
    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { i2_t ip = (i2_t){{ x, y }};
            if (y > 0)
              { /* Add vert edge from {x,y} to {x,y-1}: */
                i2_t ip_vert =  (i2_t){{ x, y-1 }};
                add_edge(&ip, &ip_vert);
              }
            if (x > 0)
              { /* Add horz edge from {x,y} to {x-1,y}: */
                i2_t ip_horz =  (i2_t){{ x, y-1 }};
                add_edge(&ip, &ip_horz);
              }
          }
      }

    return gr;

    void add_edge(i2_t *ipo, i2_t *ipd)
      { psgr_vertex_t org = psgr_vertex_from_pixel(gr, ipo);
        demand((org.ix != VNONE.ix) && (org.ix >= 0) && (org.ix < NV), "invalid origin {orgx,orgy}");
        psgr_vertex_t dst = psgr_vertex_from_pixel(gr, ipd);
        demand((dst.ix != VNONE.ix) && (dst.ix >= 0) && (dst.ix < NV), "invalid destination {dstx,dsty}");
        double we, de;
        get_edge_parameters(ipo, ipd, &de, &we);
        if (we > 0)
          { psgr_arc_t ai = psgr_add_edge(gr, org,dst, de,we, psgr_path_NULL);
            assert(ai.ix != ANONE.ix);
          }
      }

    void get_edge_parameters(i2_t *ipu, i2_t *ipv,  double *d_P, double *w_P)
      { int32_t xu = ipu->c[0], yu = ipu->c[1];
        int32_t xv = ipv->c[0], yv = ipv->c[1];
        pst_slope_map_get_edge_data(IG, xu,yu, xv-xu, yv-yu, extrapolate, d_P, w_P);
      }
      
  }
