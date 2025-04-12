/* See {pst_graph_from_maps.h}. */
/* Last edited on 2025-03-15 21:30:48 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <float_image.h>

#include <pst_interpolate.h>

#include <pst_graph.h>
#include <pst_graph_from_maps.h>

void pst_graph_interpolate_two_samples
  (  float_image_t* I, float_image_t* W,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *v, double* w
   );

/* IMPLEMENTATIONS */

int32_t pst_graph_compute_vertex_index(int32_t ix, int32_t iy, int32_t NX, int32_t NY)
  {
    if ((iy < 0 ) || (iy > NY)) { return -1; }
    if ((ix < 0 ) || (ix > NX)) { return -1; }
    return ix + iy*NX;
  }

void pst_graph_restore_vertex_index(int32_t id, int32_t NX, int32_t NY, int32_t *ix, int32_t *iy)
  {
    demand((NX >= 1) && (NY >= 1), "invalid image size");
    if ((id < 0 )  || (id >= NX*NY)) 
      { *ix = *iy = -1; }
    else
      { *ix = id%(int32_t)NX; *iy = id/(int32_t)NX; }
  }

void pst_graph_interpolate_two_samples
  (  float_image_t *I, float_image_t *W,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *v, double *w
   )
   {
     int32_t NX = (int32_t)I->sz[1]; 
     if (W != NULL) { assert(W->sz[1] == NX); }
     int32_t NY = (int32_t)I->sz[2]; 
     if (W != NULL) { assert(W->sz[2] == NY); }
     
     double v0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : float_image_get_sample(I, (int32_t)c, x0, y0));
     double v1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : float_image_get_sample(I, (int32_t)c, x1, y1));
     double w0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x0, y0)));
     double w1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x1, y1)));
   
     /* First we get the interpolation of one diagonal*/
     *v = (w0*v0+w1*v1)/2;
     *w = (w0+w1)/2;
   }

pst_graph_t *pst_graph_create_from_gradient_and_weight_maps(float_image_t *IG, float_image_t *IW)
  {
    assert(IG->sz[0] == 2);
    int32_t NX = (int32_t)IG->sz[1];
    int32_t NY = (int32_t)IG->sz[2];

    int32_t NX_Z = NX+1;
    int32_t NY_Z = NY+1;
    
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert((NX ==  IW->sz[1]) && (NY == IW->sz[2]));
      }

    uint32_t NV_max = (uint32_t)(NX_Z*NY_Z);
    uint32_t NE_max = (uint32_t)(2*(NX_Z*NY_Z) - NX_Z - NY_Z);

    pst_graph_t *g = pst_graph_new(NV_max, NE_max);

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { int32_t ioo = pst_graph_compute_vertex_index(x, y, NX_Z, NY_Z);
            assert((ioo >= 0) && (ioo < NV_max));
            pst_graph_add_vertex(g, (uint32_t)ioo, x, y);

            /* Edge from {(x,y)} to {(x-1,y)}: */
            double dxm, wxm;
            pst_map_interpolate_samples(IG, IW, 0, 2, x-1, y-1, x-1, y+0, &dxm, &wxm);
            if (wxm > 0) 
              { int32_t imo = pst_graph_compute_vertex_index(x-1, y, NX_Z, NY_Z);
                assert((imo >= 0) && (imo < NV_max));
                pst_graph_add_edge(g, (uint32_t)ioo, (uint32_t)imo, dxm, wxm, 0);
              }

            /* Edge from {(x,y)} to {(x,y-1)}: */
            double dym, wym;
            pst_map_interpolate_samples(IG, IW, 1, 2, x-1, y-1, x+0, y-1, &dym, &wym);
            if (wym > 0) 
              { int32_t iom = pst_graph_compute_vertex_index(x, y-1, NX_Z, NY_Z);
                assert((iom >= 0) && (iom < NV_max));
                pst_graph_add_edge(g, (uint32_t)ioo, (uint32_t)iom, dym, wym, 1);
              }
          }
      }

    pst_graph_update_neighbours(g);
    return g;
  }
