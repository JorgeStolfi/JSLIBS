/* See {pst_img_graph_from_maps.h} */
/* Last edited on 2025-02-11 10:11:16 by stolfi */
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

#include <pst_slope_map.h>
#include <pst_imgsys.h>
#include <pst_img_graph.h>

#include <pst_img_graph_from_maps.h>

void pst_img_graph_get_diagonal_edge_data_from_maps
  ( float_image_t* IG,
    float_image_t* IW,
    int32_t x, int32_t y,
    int32_t dirx,int32_t diry,
    double *d, double*w
  );

/* IMPLEMENTATIONS */

int32_t pst_img_graph_from_maps_get_vertex_index_from_height_map_indices(int32_t ix, int32_t iy, int32_t NX_Z, int32_t NY_Z)
  { if ((iy < 0) || (iy >= NY_Z)) { return -1; }
    if ((ix < 0) || (ix >= NX_Z)) { return -1; }
    int32_t xy = ix + iy*NX_Z;
    return xy;
  }

void pst_img_graph_from_maps_get_height_map_indices_from_vertex_index(int32_t xy, int32_t NX_Z, int32_t NY_Z, int32_t *ix, int32_t *iy)
  { int32_t NXY_Z = NX_Z*NY_Z;
    if ((xy < 0) || (xy >= NXY_Z)) 
      { (*ix) = -1; (*iy) = -1; }
    else
      { (*ix) = xy % NX_Z;
        (*iy) = xy / NX_Z;
      }
  }

pst_img_graph_t* pst_img_graph_from_gradient_and_weight_maps(float_image_t* IG,  bool_t add_diags)
  { 
    demand(! add_diags, "{add_diags} not implemented yet");
    demand(IG != NULL, "gradient map must not be {NULL}");
    
    int32_t NC_G, NX_G, NY_G; /* Cols and rows of {IG} and {IW}. */
    float_image_get_size(IG, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "gradient map must have 3 channels");

    /* Cols and rows of height map: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;

    /* The vertices are the height map pixels: */
    int32_t NV = NX_Z*NY_Z;
    
    /* The edges connect horz or vert adjacent adjacent height map pixels: */
    int32_t NE = NX_Z*NY_G + NX_G*NY_Z;
    
    pst_img_graph_t* g = pst_img_graph_new((uint32_t)NV, (uint32_t)NE);

    auto haf_arc_t add_edge(int32_t orgx, int32_t orgy, int32_t dstx, int32_t dsty);
      /* Adds to the graph {g} an edge betewwn vertex {org=(orgx,orgy)}
        and vertex {dst=(dstx,dsty)}. where {orgx,dstx} are in
        {0..NX_Z-1} and {orgy,dsty} are in {0..NY_Z-1}. Creates the
        corresponding {haf_edge_t} returns the {haf_arc_t} from {org} to
        {dst}. Does NOT splice the arc into the half-edge structure. */

    auto void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d, double *w);
      /* Sets {*d} to the height difference of the graph edge from
        vertex {u=(ux,uy)} to vertex {v=(vx,y)}, as estimated from the
        gradient map {IG}. Also sets {*w} to the corresponding weight,
        estimated from the weight map {IW}; if {IW} is {NULL}, {*w} is
        set to 1. Note that {u} and {v} correspond to pixels of the
        height map, hence to corners of the gradient and weight map
        grids. */

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { r2_t coords = (r2_t){{ x, y }};
            uint32_t kv = pst_img_graph_add_vertex(g, x, y, NULL, coords);
            assert(kv == pst_img_graph_from_maps_get_vertex_index_from_height_map_indices(x,y, NX_Z,NY_Z));
          }
      }

    /* The west-pointing arc out of vertex {x} in prev/current row is {west[x]}: */
    haf_arc_t *west = talloc(NX_Z, haf_arc_t);
    
    /* The south-pointing arc out of vertex {x} in prev/current row is {south[x]}: */
    haf_arc_t *south = talloc(NX_Z, haf_arc_t);
    for (int32_t x = 0; x < NX_Z; x++) { west[x] = south[x] = NULL; }
    
    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { if ((y > 0) && (x > 0) && (x < NX_Z-1))
              { /* Add vert edge from {x,y} to {x,y-1}: */
                haf_arc_t ev = add_edge(x,y, x,y-1);
                haf_arc_t ew = west[x];
                if (ew != NULL) { haf_splice(ev,ew); }
                south[x] = ev;
              }
            if ((x > 0) && (y > 0) && (y < NY_Z-1))
              { /* Add horz edge from {x,y} to {x-1.y}: */
                haf_arc_t eh = add_edge(x,y, x-1,y);
                haf_arc_t es = south[x];
                if (es != NULL) { haf_splice(eh,es); }
                haf_arc_t et = south[x-1];
                if (et != NULL) { haf_splice(haf_sym(eh), et); }
                west[x] = eh;
              }
          }
      }

    return g;

    haf_arc_t add_edge(int32_t orgx,int32_t orgy, int32_t dstx, int32_t dsty)
      { int32_t org = pst_img_graph_from_maps_get_vertex_index_from_height_map_indices(orgx,orgy, NX_Z,NY_Z);
        demand((org >= 0) && (org < NV), "invalid origin {orgx,orgy}");
        int32_t dst = pst_img_graph_from_maps_get_vertex_index_from_height_map_indices(dstx,dsty, NX_Z,NY_Z);
        demand((dst >= 0) && (dst < NV), "invalid origin {orgx,orgy}");
        char* elabel = jsprintf("from %ld=[%ld,%ld] to %ld=[%ld,%ld]", org,orgx,orgy, dst,dstx,dsty);
        double we, de;
        get_edge_parameters(orgx,orgy, dstx,dsty, &de, &we);
        haf_arc_t e = pst_img_graph_edge_add(g, (uint32_t)org,(uint32_t)dst, de,we, elabel, pst_path_create_empty());
        return e;
      }

    void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d, double *w)
      { pst_slope_map_get_edge_data(IG, xu, yu, xv, yv, d, w); }
      
  }
  
