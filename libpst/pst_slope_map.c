/* See pst_slope_map.h */
/* Last edited on 2025-03-01 02:57:56 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <jsrandom.h>

#include <pst_basic.h>
#include <pst_interpolate.h>
#include <pst_cell_map_shrink.h>

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

void pst_slope_map_interpolate_four_samples
  (  float_image_t *G,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *vR_P, double *wR_P
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
     
     auto void fetch(int32_t x, int32_t y, double *v_P,double *w_P);
       /* Fetches the value {v = G[c,X,Y]} and the weight {w=G[2,X,Y]},
         returning them in {*v_P} and {*w_P}. However, if pixel {x,y}
         does not exist, or {v} is not finite, of {w} is zero, sets
         {v=NAN} and {w=0}. */
     
     /* Extract the sample values and weights: */
     double vm, wm; fetch(xm, ym, &vm, &wm);
     double v0, w0; fetch(x0, y0, &v0, &w0);
     double v1, w1; fetch(x1, y1, &v1, &w1);
     double vp, wp; fetch(xp, yp, &vp, &wp);
     pst_interpolate_four_values(vm, wm, v0, w0, v1, w1, vp, wp, vR_P, wR_P);

     return;
     
     void fetch(int32_t x, int32_t y, double *v_P,double *w_P)
       { double v = NAN, w = 0.0;
         if ((x >= 0) && (y >= 0) && (x < NX) && (y < NY))
           { w = (NC == 3 ? float_image_get_sample(G,2,x,y) : 1.0); 
             assert(isfinite(w) && (w >= 0));
             if (w > 0)
               { v = float_image_get_sample(G,c,x,y);
                 if (! isfinite(v)) { v = NAN; w = 0; }
               }
           }
         (*v_P) = v;
         (*w_P) = w;
       }
   }

void  pst_slope_map_get_edge_data
  ( float_image_t* G,
    int32_t x, int32_t y,
    int32_t ux, int32_t uy,
    double *dP, double *wP
  )
  {
    demand((ux != 0) || (uy != 0), "invalid increments {ux,uy}");
    
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
        
        /* Get the gradient {gx,gy} and weight {w} in cell with corners {(x,y),(x',y')}: */
        int32_t xg = (ux > 0 ? x : x-1);
        int32_t yg = (uy > 0 ? y : y-1);
        
        double gx, gy;
        if ((xg >= 0) && (xg < NX) && (yg >= 0) && (yg < NY))
          { float v[NC];
            float_image_get_pixel(G, xg, yg, v);
            gx = v[0]; gy = v[1]; w = (NC >= 3 ? v[2] : 1.0);
            demand(isfinite(w) && (w >= 0), "invalid gradient weight");
            if ((! isfinite(gx)) || (! isfinite(gy))) { w = 0.0; }
          }
        else
          { w = 0.0; }
        
        /* Compute the height difference {d} from {(x,y)} and {(x',y')}: */
        if (w == 0)
          { d = NAN; }
        else
          { d = ux*gx + uy*gy;
            if (! isfinite(d)) { d = NAN; w = 0.0; }
          }
     }
   (*dP) = d; (*wP) = w;
   return;
  }

float_image_t *pst_slope_map_shrink(float_image_t *IG, double scale)
  { 
    int32_t NC, NXI, NYI;
    float_image_get_size(IG, &NC, &NXI, &NYI);
    demand((NC == 2) || (NC == 3), "slope map must have 2 or 3 channels");

    /* Shrunk map dimensions: */
    int32_t NXJ = (NXI+1)/2;
    int32_t NYJ = (NYI+1)/2;

    /* Allocate the reduced map: */
    float_image_t *JG = pst_cell_map_shrink(IG, 2, NXJ, NYJ, scale);
    return JG;
  }
  
float_image_t* pst_slope_map_merge(float_image_t *GA, float_image_t *GB)
  {
    int32_t NC_A, NX, NY;
    float_image_get_size(GA, &NC_A, &NX, &NY);
    demand((NC_A == 2) || (NC_A == 3), "slope map {GA} must have 2 or 3 channels");
    
    int32_t NC_B;
    float_image_get_size(GB, &NC_B, NULL, NULL);
    demand((NC_B == 2) || (NC_B == 3), "slope map {GB} must have 2 or 3 channels");

    float_image_check_size(GB, -1, NX, NY, "map size mismatch {GA,GB}");
    
    float_image_t *G = float_image_new(3, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { r2_t dzA = pst_slope_map_get_gradient(GA, x, y);
            r2_t dzB = pst_slope_map_get_gradient(GB, x, y);
            double wA = pst_slope_map_get_weight(GA, x, y); demand(isfinite(wA) && (wA >= 0), "invalid {A} weight");
            double wB = pst_slope_map_get_weight(GB, x, y); demand(isfinite(wB) && (wB >= 0), "invalid {B} weight");
            r2_t dz; double w;
            if ((wA == 0) || (wB == 0) || (! r2_is_finite(&dzA)) || (! r2_is_finite(&dzB)))
              { w = 0.0; dz = (r2_t){{ NAN, NAN }}; }
            else
              { w = 2.0/(1.0/wA + 1.0/wB);
                if ((! isfinite(w)) || (w < 1.0e-20))
                  { w = 0.0; dz = (r2_t){{ NAN, NAN }}; }
                else
                  { double wtot = wA + wB;
                    r2_mix(wA/wtot, &dzA, wB/wtot, &dzB, &dz);
                  }
              }
            pst_slope_map_set_gradient(G, x, y, &dz);
            pst_slope_map_set_weight(G, x, y, (float)w);
          }
      }
      
    return G;
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
