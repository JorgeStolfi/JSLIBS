/* See pst_slope_map.h */
/* Last edited on 2025-01-07 22:07:47 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <jsfile.h>
#include <r2.h>
#include <jswsize.h>
#include <rn.h>

#include <pst_basic.h>
#include <pst_interpolate.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_slope_map.h>

r2_t pst_slope_map_get_pixel(float_image_t *G, int32_t x, int32_t y)
  { demand(G->sz[0] == 2, "wrong slope map depth");
    r2_t grd;
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { grd.c[axis] = (*p); p += G->st[0]; }
    return grd;
  }

void pst_slope_map_set_pixel(float_image_t *G, int32_t x, int32_t y, r2_t *grd)
  { demand(G->sz[0] == 2, "wrong slope map depth");
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { (*p) = (float)(grd->c[axis]); p += G->st[0]; }
  }

void  pst_slope_map_get_axial_edge_data
  ( float_image_t* G,
    float_image_t* W,
    int32_t x, int32_t y,
    int32_t axis, int32_t dir,
    double *dP, double *wP
  )
  {
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 2, "slope map must have 2 channels");
    if (W == NULL) { float_image_check_size(W, 1, NX_G, NY_G); }
    demand((dir == -1) || (dir == +1), "invalid direction {dir}");
    
    demand((x >= 0) && (x <= NX_G), "invalid grid corner {X} (1)");
    demand((y >= 0) && (y <= NY_G), "invalid grid corner {Y} (1)");
    
    int32_t x1 = x + (axis == 0 ? dir : 0);
    int32_t y1 = y + (axis == 1 ? dir : 0);
    
    
    demand((x1 >= 0) && (x1 <= NX_G), "invalid grid corner {X} (2)");
    demand((y1 >= 0) && (y1 <= NY_G), "invalid grid corner {Y} (2)");
  
    double d,w;
    if (axis == 0)
      { if (dir == +1){ pst_interpolate_four_samples(G,W, 0, x+0,y+0, x+0,y-1, &d,&w); }
        if (dir == -1){ pst_interpolate_four_samples(G,W, 0, x-1,y+0, x-1,y-1, &d,&w); }
      }
    else if (axis == 1)
      { if (dir == +1){ pst_interpolate_four_samples(G,W, 1, x-1,y+0, x+0,y+0, &d,&w); }
        if (dir == -1){ pst_interpolate_four_samples(G,W, 1, x-1,y-1, x+0,y-1, &d,&w); }
      }
    else
      { demand(FALSE,"invalid  AXIS"); }

    /* Reverse the sign of the difference depending on {dir}: */
    if (dir == -1){ d = -d; }

    /* If the weight is zero, the difference is irrelevant, so clear it: */
    if (w == 0){ d = 0; }

    (*dP) = d; (*wP) = w;
  }

void  pst_slope_map_get_edge_data
  ( float_image_t* G,
    float_image_t* W,
    int32_t x, int32_t y,
    int32_t ux, int32_t uy,
    double *dP, double *wP
  )
  {
    demand((ux != 0) || (uy != 0), "invalid increments {ux,uy}");

    auto void get_two_edge_path_data(int32_t axis, double *dpP, double *wpP);
      /* Estimates the difference alog the two-edge path from {x,y} to {x',y'}
        that begins along the specified {axis}. */
    
    double d, w;
    if (ux == 0)
      { pst_slope_map_get_axial_edge_data(G, W, x, y, 1, uy, &d, &w); }
    else if (uy == 0)
      { pst_slope_map_get_axial_edge_data(G, W, x, y, 0, ux, &d, &w); }
    else
      { assert((ux == -1) || (ux == +1));
        assert((uy == -1) || (uy == +1));
        int32_t NC_G, NX_G, NY_G;
        float_image_get_size(G, &NC_G, &NX_G, &NY_G);
        demand(NC_G == 2, "slope map must have 2 channels");
        if (W == NULL) { float_image_check_size(W, 1, NX_G, NY_G); }
        
        /* Walk horizontally then vertically: */
        double dp0, wp0;
        get_two_edge_path_data(0, &dp0, &wp0);
        
        /* Walk vertically then horizontally: */
        double dp1, wp1;
        get_two_edge_path_data(1, &dp1, &wp1);
        
        /* Combine the two paths: */
        w = wp0 + wp1;
        if (w > 0)
          { d = (wp0*dp0 + wp1*dp1)/w; }
        else
          { /* Indeterminate, might as well: */ d = 0.0; }
          
          
     }
   (*dP) = d; (*wP) = w;
   return;
   
   void get_two_edge_path_data(int32_t axis, double *dpP, double *wpP)
     { 
       /* Get differences and weights of the two edges: */
       double d1, w1, d2, w2;
       if (axis == 0)
         { pst_slope_map_get_axial_edge_data(G, W, x, y, 0, ux, &d1, &w1);
           pst_slope_map_get_axial_edge_data(G, W, x+ux, y, 1, uy, &d2, &w2);
         }
       else if (axis == 1)
         { pst_slope_map_get_axial_edge_data(G, W, x, y, 1, uy, &d1, &w1);
           pst_slope_map_get_axial_edge_data(G, W, x, y+uy, 0, ux, &d2, &w2);
         }
       else
         { assert(FALSE); }
         
       /* Combine them: */
       double d = d1 + d2;
       double w = w1*w2;
       if (w == 0) { /* Might as well: */ d = 0; }
       (*dpP) = d; (*wpP) = w;
     }
  }

void pst_slope_and_weight_map_shrink
  ( float_image_t *IG, 
    float_image_t *IW, 
    float_image_t **SG, 
    float_image_t **SW
  )
  { 
    int32_t NC = (int32_t)IG->sz[0];
    int32_t NXI = (int32_t)IG->sz[1]; int32_t NXJ = (NXI+1)/2;
    int32_t NYI = (int32_t)IG->sz[2]; int32_t NYJ = (NYI+1)/2;
    float_image_t *JG = float_image_new((int32_t)NC, (int32_t)NXJ, (int32_t)NYJ);
    float_image_t *OW = NULL;
    
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert(IW->sz[1] == NXI);
        assert(IW->sz[2] == NYI);
	OW = float_image_new(1, (int32_t)NXJ, (int32_t)NYJ);
      }

    for (int32_t jy = 0; jy < NYJ; jy++)
      { for (int32_t jx = 0; jx < NXJ; jx++)
          { int32_t ix = 2*jx, iy = 2*jy;
	    float wmin = INF;
	    for (int32_t c = 0; c < 2; c++)
	    {
	      double da,wa;
	      pst_interpolate_two_samples(IG, IW, (int32_t)c, ix,iy,  ix+1,iy+1, &da,&wa);
	      double db,wb;
 	      pst_interpolate_two_samples(IG, IW, (int32_t)c, ix,iy+1, ix+1,iy,  &db,&wb);
	      
	      float d = (float)( wa+wb > 0 ? (wa*da + wb*db)/(wa + wb) : (da+db)/2.0);
	      float_image_set_sample(JG, (int32_t)c, jx, jy, d);
	      float w = (float)(4/(1/wa + 1/wb));
	      if (w < wmin) wmin = w;
	    }
   
	    if (OW != NULL) { float_image_set_sample(OW, 0, jx, jy, wmin);}
          
          }
      }
    *SG = JG;
    *SW = OW;
  }
