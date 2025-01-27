/* See pst_slope_map.h */
/* Last edited on 2025-01-23 13:20:12 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <jsfile.h>
#include <r2.h>
#include <jswsize.h>
#include <jsrandom.h>
#include <rn.h>

#include <pst_basic.h>
#include <pst_interpolate.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_slope_map.h>

r2_t pst_slope_map_get_gradient(float_image_t *G, int32_t x, int32_t y)
  { demand((G->sz[0] == 2) || (G->sz[0] == 3), "slope map must have 2 or 3 channels");
    r2_t grd;
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { grd.c[axis] = (*p); p += G->st[0]; }
    return grd;
  }

float pst_slope_map_get_weight(float_image_t *G, int32_t x, int32_t y)
  { demand((G->sz[0] == 2) || (G->sz[0] == 3), "slope map must have 2 or 3 channels");
    float w = (G->sz[0] == 3 ? float_image_get_sample(G, 2, x, y) : 1.0f);
    return w;
  }

void pst_slope_map_set_gradient(float_image_t *G, int32_t x, int32_t y, r2_t *grd)
  { demand((G->sz[0] == 2) || (G->sz[0] == 3), "slope map must have 2 or 3 channels");
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { (*p) = (float)(grd->c[axis]); p += G->st[0]; }
  }

void pst_slope_map_set_weight(float_image_t *G, int32_t x, int32_t y, float w)
  { demand(G->sz[0] == 3, "slope map must have 3 channels");
    float_image_set_sample(G, 2, x, y, w);
  }

void pst_slope_map_get_axial_edge_data
  ( float_image_t* G,
    int32_t x, int32_t y,
    int32_t axis, int32_t dir,
    double *dP, double *wP
  )
  {
    int32_t NC, NX, NY;
    float_image_get_size(G, &NC, &NX, &NY);
    demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");
    demand((dir == -1) || (dir == +1), "invalid direction {dir}");
    
    demand((x >= 0) && (x <= NX), "invalid grid corner {X} (1)");
    demand((y >= 0) && (y <= NY), "invalid grid corner {Y} (1)");
    
    int32_t x1 = x + (axis == 0 ? dir : 0);
    int32_t y1 = y + (axis == 1 ? dir : 0);
    
    
    demand((x1 >= 0) && (x1 <= NX), "invalid grid corner {X} (2)");
    demand((y1 >= 0) && (y1 <= NY), "invalid grid corner {Y} (2)");
  
    double d,w;
    if (axis == 0)
      { if (dir == +1){ pst_slope_map_interpolate_four_samples(G, 0, x+0,y+0, x+0,y-1, &d,&w); }
        if (dir == -1){ pst_slope_map_interpolate_four_samples(G, 0, x-1,y+0, x-1,y-1, &d,&w); }
      }
    else if (axis == 1)
      { if (dir == +1){ pst_slope_map_interpolate_four_samples(G, 1, x-1,y+0, x+0,y+0, &d,&w); }
        if (dir == -1){ pst_slope_map_interpolate_four_samples(G, 1, x-1,y-1, x+0,y-1, &d,&w); }
      }
    else
      { demand(FALSE,"invalid  AXIS"); }

    /* Reverse the sign of the difference depending on {dir}: */
    if (dir == -1){ d = -d; }

    /* If the weight is zero, the difference is irrelevant, so clear it: */
    if (w == 0){ d = 0; }

    (*dP) = d; (*wP) = w;
  }

void pst_slope_map_interpolate_two_samples
  (  float_image_t *G,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *vRP, double *wRP
   )
   {
     int32_t NC, NX, NY;
     float_image_get_size(G, &NC, &NX, &NY);
     demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");
     demand((c == 0) || (c == 1), "invalid channel index {c}");
    
     /* Extract the sample values and weights: */
     double v0, v1;
     double w0, w1;
     
     if ((x0 >= 0) && (y0 >= 0) && (x0 < NX) && (y0 < NY))
       { v0 = float_image_get_sample(G,c,x0,y0);
         w0 = (NC == 3 ? float_image_get_sample(G,2,x0,y0) : 1.0); 
       }
     else
       { v0 = 0; w0 = 0; }
       
     if ((x1 >= 0) && (y1 >= 0) && (x1 < NX) && (y1 < NY))
       { v1 = float_image_get_sample(G,c,x1,y1);
         w1 = (NC == 3 ? float_image_get_sample(G,2,x1,y1) : 1.0); 
       }
     else
       { v1 = 0; w1 = 0; }
       
     /* Replicate the other value if weight is zero: */
     if (w0 == 0) { v0 = v1; }
     if (w1 == 0) { v1 = v0; }
       
     /* Interpolate: */
     pst_interpolate_two_values(v0, w0, v1, w1, vRP, wRP); 
   }


void pst_slope_map_interpolate_four_samples
  (  float_image_t *G,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *vRP, double *wRP
   )
   {
     int32_t NC, NX, NY;
     float_image_get_size(G, &NC, &NX, &NY);
     demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");
     demand((c == 0) || (c == 1), "invalid channel index {c}");
     
     int32_t dx = x1 - x0, dy = y1 - y0;
     demand((abs(dx) <= 1) && (abs(dy) <= 1) && (abs(dx) + abs(dy) == 1), "pixels are not adjacent");
  
     /* Compute indices of auxiliary pixels: */
     int32_t xm = 2*x0 - x1, ym = 2*y0 - y1;
     int32_t xp = 2*x1 - x0, yp = 2*y1 - y0;
     
     /* Extract the sample values and weights: */
     double vm, v0, v1, vp;
     double wm, w0, w1, wp;
     
     if ((xm >= 0) && (ym >= 0) && (xm < NX) && (ym < NY))
       { vm = float_image_get_sample(G,c,xm,ym);
         wm = (NC == 3 ? float_image_get_sample(G,2,xm,ym) : 1.0); 
       }
     else
       { vm = 0; wm = 0; }
     
     if ((x0 >= 0) && (y0 >= 0) && (x0 < NX) && (y0 < NY))
       { v0 = float_image_get_sample(G,c,x0,y0);
         w0 = (NC == 3 ? float_image_get_sample(G,2,x0,y0) : 1.0); 
       }
     else
       { v0 = 0; w0 = 0; }
       
     if ((x1 >= 0) && (y1 >= 0) && (x1 < NX) && (y1 < NY))
       { v1 = float_image_get_sample(G,c,x1,y1);
         w1 = (NC == 3 ? float_image_get_sample(G,2,x1,y1) : 1.0); 
       }
     else
       { v1 = 0; w1 = 0; }
       
     if ((xp >= 0) && (yp >= 0) && (xp < NX) && (yp < NY))
       { vp = float_image_get_sample(G,c,xp,yp);
         wp = (NC == 3 ? float_image_get_sample(G,2,xp,yp) : 1.0); 
       }
     else
       { vp = 0; wp = 0; }
	    
     /* Interpolate: */
     pst_interpolate_four_values(vm,wm, v0,w0, v1,w1, vp,wp, vRP,wRP);
   }

void  pst_slope_map_get_edge_data
  ( float_image_t* G,
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
      { pst_slope_map_get_axial_edge_data(G, x, y, 1, uy, &d, &w); }
    else if (uy == 0)
      { pst_slope_map_get_axial_edge_data(G, x, y, 0, ux, &d, &w); }
    else
      { assert((ux == -1) || (ux == +1));
        assert((uy == -1) || (uy == +1));
        int32_t NC, NX, NY;
        float_image_get_size(G, &NC, &NX, &NY);
        demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");
        
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
         { pst_slope_map_get_axial_edge_data(G, x,    y, 0, ux, &d1, &w1);
           pst_slope_map_get_axial_edge_data(G, x+ux, y, 1, uy, &d2, &w2);
         }
       else if (axis == 1)
         { pst_slope_map_get_axial_edge_data(G, x, y,    1, uy, &d1, &w1);
           pst_slope_map_get_axial_edge_data(G, x, y+uy, 0, ux, &d2, &w2);
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

float_image_t *pst_slope_map_shrink(float_image_t *IG)
  { 
    int32_t NC, NXI, NYI;
    float_image_get_size(IG, &NC, &NXI, &NYI);
    demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");

    int32_t NXJ = (NXI+1)/2;
    int32_t NYJ = (NYI+1)/2;
    float_image_t *JG = float_image_new(NC, (int32_t)NXJ, (int32_t)NYJ);

    for (int32_t yJ = 0; yJ < NYJ; yJ++)
      { for (int32_t xJ = 0; xJ < NXJ; xJ++)
          { double sum_wv0 = 0;
            double sum_wv1 = 0;
            double sum_w = 0; 
            float min_w = +INF;
            for (int32_t dx = 0; dx <= 1; dx++)
              { for (int32_t dy = 0; dy <= 1; dy++)
                  { int32_t xI = 2*xJ + dx;
                    int32_t yI = 2*yJ + dy;
                    if ((xI < NXI) && (yI < NYI))
                      { float vI[NC];
                        float_image_get_pixel(IG, xI, yI, vI);
                        float w = (NC == 3 ? vI[2] : 1.0f);
                        sum_wv0 += w*vI[0];
                        sum_wv1 += w*vI[1];
                        sum_w += w;
                        min_w = fminf(min_w, w);
                      }
                  }
              }
            float vJ[3];
            vJ[0] = (float)(sum_w == 0 ? 0.0 : sum_wv0/sum_w);
            vJ[1] = (float)(sum_w == 0 ? 0.0 : sum_wv1/sum_w);
            if (NC == 3) { vJ[2] = min_w; }
            float_image_set_pixel(JG, xJ, yJ, vJ);
          }
      }
    return JG;
  }

void pst_slope_map_perturb(float_image_t *G, double sigma, uint32_t seed)
  { 
    int32_t NC, NX, NY;
    float_image_get_size(G, &NC, &NX, &NY);
    demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");
    
    if (seed != 0) { srand(seed); srandom(seed); }
    
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { for (int32_t c = 0; c < 2; c++) 
              { float *smp = float_image_get_sample_address(G, c, x, y);
                (*smp) = (float)((*smp) + sigma*dgaussrand());
              }
          }
      }
  }
