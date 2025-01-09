/* See {pst_img_graph_from_maps.h} */
/* Last edited on 2025-01-06 18:47:54 by stolfi */
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

#include <pst_img_graph_from_maps.h>

void pst_img_graph_get_diagonal_edge_data_from_maps
  ( float_image_t* IG,
    float_image_t* IW,
    int32_t x, int32_t y,
    int32_t dirx,int32_t diry,
    double *d, double*w
  );

/* IMPLEMENTATIONS */

int32_t pst_img_graph_get_vertex_index_from_height_map_indices(int32_t ix, int32_t iy, int32_t NX_Z, int32_t NY_Z)
  { if ((iy < 0) || (iy >= NY_Z)) { return -1; }
    if ((ix < 0) || (ix >= NX_Z)) { return -1; }
    return ix + iy*NX_Z;
  }

void pst_img_graph_get_height_map_indices_from_vertex_index(int32_t vid, int32_t NX_Z, int32_t NY_Z, int32_t *ix, int32_t *iy)
  { int32_t NV = NX_Z*NY_Z;
    if ((vid < 0) || (vid >= NV)) { (*ix) = -1; (*iy) = -1; return; }
    (*ix) = vid % NX_Z;
    (*iy) = vid / NX_Z;
  }

pst_img_graph_t* pst_img_graph_from_gradient_and_weight_maps
  ( float_image_t* IG,
    float_image_t* IW,
    bool_t add_diags
  )
  { 
    demand(! add_diags, "{add_diags} not implemented yet");
    demand(IG != NULL, "gradient map must not be {NULL}");
    
    int32_t NC_G, NX_G, NY_G; /* Cols and rows of {IG} and {IW}. */
    float_image_get_size(IG, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 2, "gradient map must have 2 channels");
    if (IW != NULL) { float_image_check_size(IW, 1, NX_G, NY_G); }

    /* Cols and rows of height map: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;

    /* The vertices are the height map pixels: */
    int32_t NV = NX_Z*NY_Z;
    
    /* The edges connect horz or vert adjacent adjacent height map pixels: */
    int32_t NE = NX_Z*NY_G + NX_G*NY_Z;
    
    pst_img_graph_t* g = pst_img_graph_new(NV, NE);

    auto haf_arc_t add_edge(int32_t orgx, int32_t orgy, int32_t dstx, int32_t dsty);
      /* Adds to the graph {g} an edge betewwn vertex {org=(orgx,orgy)} and vertex {dst=(dstx,dsty)}.
        where {orgx,dstx} are in {0..NX_Z-1} and {orgy,dsty} are
        in {0..NY_Z-1}.  Creates the corresponding {haf_edge_t} returns the {haf_arc_t} from
        {org} to {dst}.  Does NOT splice the arc into the half-edge structure.  */

    auto void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d, double *w);
      /* Sets {*d} to the height difference of the graph edge from vertex {u=(ux,uy)}
        to vertex {v=(vx,y)}, as estimated from the gradient map {IG}.  Also sets {*w} to the
        corresponding weight, estimated from the weight map {IW}; if {IW} is {NULL},
        {*w} is set to 1.  Note that {u} and {v} correspond to pixels of the height map,
        hence to corners of the gradient and weight map grids. */

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { int32_t id = pst_img_graph_get_vertex_index_from_height_map_indices(x,y,NX_Z,NY_Z);
            r2_t coords = (r2_t){{x,y}};
            pst_img_graph_add_vertex(g, id, NULL, coords);
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
      { int32_t org = pst_img_graph_get_vertex_index_from_height_map_indices(orgx,orgy, NX_Z,NY_Z);
        int32_t dst = pst_img_graph_get_vertex_index_from_height_map_indices(dstx,dsty, NX_Z,NY_Z);
        char* elabel = jsprintf("from %ld (%ld,%ld) to %ld(%ld,%ld)", org,orgx,orgy, dst,dstx,dsty);
        double we, de;
        get_edge_parameters(orgx,orgy, dstx,dsty, &de, &we);
        haf_arc_t e = pst_img_graph_add_edge(g, org,dst, de,we, elabel, pst_path_create_empty());
        return e;
      }

    void get_edge_parameters(int32_t xu, int32_t yu, int32_t xv, int32_t yv,  double *d, double *w)
      {
        int32_t dx = xv - xu;
        int32_t dy = yv - yu;
        if (dx == 0) 
          { assert(abs(dy) == 1);
            pst_img_graph_get_axial_edge_data_from_maps(IG,IW, xu,yu, Y_AXIS, dy, d,w);
          }
        else if (dy == 0)
          { assert(abs(dx) == 1);
            pst_img_graph_get_axial_edge_data_from_maps(IG,IW, xu,yu, X_AXIS, dx, d,w);
          }
        else
          { assert(abs(dx) == 1);
            assert(abs(dy) == 1);
            pst_img_graph_get_diagonal_edge_data_from_maps(IG,IW, xu,yu, dx,dy, d,w);
          }
        /* Is this possible? else{ *d = 0; *w =  1/r2_dist_sqr(&u,&v);  } */

        return;

        void get_nearest_grid_corner(r2_t pt, int32_t NX, int32_t NY, int32_t* x, int32_t* y)
            {
          assert( (pt.c[0] >= 0) && (pt.c[0] <= NX ));
          assert( (pt.c[1] >= 0) && (pt.c[1] <= NY ));
          *x = pt.c[0];
          *y = pt.c[1];
        }

      }
      
  }
  
void pst_img_graph_compute_edge_parameters
  ( int32_t NX, int32_t NY,
    int32_t gradientFunction,
    int32_t weightFunction,
    r2_t u, r2_t v,
    double *d, double *w
  )
  {
    double ww;
    r2_t p;
    r2_mix(0.5,&u,0.5,&v,&p);

    switch(weightFunction)
        {
      case WGHT_FUNC_NONE :
        break;
      case WGHT_FUNC_CONST :
        *w = 2.5;
        break;
      case WGHT_FUNC_RAMP :
        *w = (0.25/NX)*p.c[0] + (0.75/NY)*p.c[1] ;
        break;
      case WGHT_FUNC_RANDOM :
        ww = (sin(17*p.c[0]*p.c[0]) + cos(29*p.c[1]*p.c[1]))*47;
        (*w) = ww - floor(ww);
        break;
      default:
        demand(FALSE,"Invalid Weight function");
    };

    r2_t uv;
    r2_sub(&v,&u,&uv);
    double dx,dy,dd;
    switch(gradientFunction)
        {
      case GRAD_FUNC_NONE :
        break;
      case GRAD_FUNC_CONST :
        *d = 2.5;
        break;
      case GRAD_FUNC_RAMP :
        dx = 2*((0.25/NX)*p.c[0] + (0.75/NY)*p.c[1]) -1 ;
        dy = 2*((0.75/NX)*p.c[0] + (0.25/NY)*p.c[1]) -1 ;
        *d = dx*uv.c[0] + dy*uv.c[1];
        break;
      case GRAD_FUNC_RANDOM :
        dd = (sin(17*p.c[0]*p.c[0]) + cos(29*p.c[1]*p.c[1]))*47;
        dd = 2*(dd - floor(dd)) -1;
        dd *= r2_dist(&u,&v);
        (*d)= dd; 
        break;
      default:
        demand(FALSE,"Invalid Grad function");
    };
  }

void pst_img_graph_get_diagonal_edge_data_from_maps(
  float_image_t* IG,
  float_image_t* IW,
  int32_t x, int32_t y,
  int32_t dirx,int32_t diry,
  double *d, double*w
  )
      {
  
  if ((IG == NULL) && (IW == NULL) )
      {
    *d = 0;
    *w = 0.5;
    return;
  }
  
  assert( (abs(dirx) == 1) && (abs(diry) == 1) );
  int32_t NX =  (IG == NULL ?  IW->sz[1] : IG->sz[1]);
  int32_t NY =  (IG == NULL ?  IW->sz[2] : IG->sz[2]);
  
  int32_t ix = (dirx > 0 ? x:x-1);
  int32_t iy = (diry > 0 ? y:y-1);
  
  assert( (ix >=0) && (ix < NX));
  assert( (iy >=0) && (iy < NY));
  
  double dx = dirx*(IG == NULL ? 0 : float_image_get_sample(IG,0,ix,iy));
  double dy = diry*(IG == NULL ? 0 : float_image_get_sample(IG,1,ix,iy));
  *w = 0.5*(IW == NULL ? 1.0 : float_image_get_sample(IW,0,ix,iy));
 
  *d = dx+dy;
  if (*w == 0) { *d = 0; }
}
