/* See pst_normal_map.h */
/* Last edited on 2016-03-16 16:09:49 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <jsrandom.h>
#include <affirm.h>

#include <pst_normal_map.h>
#include <pst_slope_map.h>

r3_t pst_normal_map_get_pixel(float_image_t *NRM, int x, int y)
  { demand(NRM->sz[0] == 3, "wrong normal map depth");
    r3_t nrm;
    float *p = float_image_get_sample_address(NRM, 0, x, y);
    int axis;
    for (axis = 0; axis < 3; axis++)
      { nrm.c[axis] = (*p); p += NRM->st[0]; }
    return nrm;
  }
  
void pst_normal_map_set_pixel(float_image_t *NRM, int x, int y, r3_t *nrm)
  { demand(NRM->sz[0] == 3, "wrong normal map depth");
    float *p = float_image_get_sample_address(NRM, 0, x, y);
    int axis;
    for (axis = 0; axis < 3; axis++)
      { (*p) = (float)nrm->c[axis]; p += NRM->st[0]; }
  }

r2_t pst_normal_map_scene_pt_from_image_pt(r2_t *xy, r3x3_t *xym_to_uvm)
  { r3_t xym = (r3_t) {{ xy->c[0], xy->c[1], 1.0 }};
    r3_t uvm;
    r3x3_map_row(&xym, xym_to_uvm, &uvm);
    r2_t uv = (r2_t){{ uvm.c[0]/uvm.c[2], uvm.c[1]/uvm.c[2] }};
    return uv;
  }

r3_t pst_normal_map_image_dir_from_scene_dir(r3_t *uvw, r3x3_t *uvw_to_xyz)
  { r3_t xyz;
    if (r3_L_inf_norm(uvw) == 0.0)
      { /* Normal is undefined - return null normal: */
        xyz = (r3_t) {{ 0.0, 0.0, 0.0 }};
      }
    else
      { /* Map normal from {U,V,W} to {X,Y,Z} system: */
        r3x3_map_row(uvw, uvw_to_xyz, &xyz);
        /* Normalize: */
        (void)r3_dir(&xyz, &xyz);
      }
    return xyz;
  }

r3_t pst_normal_map_eval
  ( pst_normal_map_proc_t nrmf, /* Normal-computing funtion. */
    double x,                   /* X-coordinate of projected point in image system. */
    double y,                   /* Y-coordinate of projected point in image system. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {xy} coords to model {uv} coords. */
    r3x3_t *uvw_to_xyz          /* Linear map of {uvw} coords to normal {xyz} coords. */
  )
  { /* Convert point from {X,Y} coords to {U,V} coords: */
    r2_t pt_xy = (r2_t){{ x, y }};
    r2_t pt_uv = pst_normal_map_scene_pt_from_image_pt(&pt_xy, xym_to_uvm);

    /* Compute normal in {U,V,W} system: */
    r3_t nrm_uvw = nrmf(&pt_uv);
    
    /* Convert normal direction back to image {X,Y,Z} coords. */
    r3_t nrm_xyz =  pst_normal_map_image_dir_from_scene_dir(&nrm_uvw, uvw_to_xyz);
    
    return nrm_xyz;
  }

void pst_normal_map_from_proc
  ( pst_normal_map_proc_t nrmf, /* Normal-computing funtion. */
    int NS,                     /* Order of subsampling grid within each pixel. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {xy} coords to model {uv} coords. */
    r3x3_t *uvw_to_xyz,         /* Linear map of {uvw} coords to normal {xyz} coords. */
    float_image_t *NRM          /* (OUT) Computed normal map. */
  )
  { /* Get/check image dimensions: */
    demand(NRM->sz[0] == 3, "bad normal map");
    int NX = (int)(NRM->sz[1]);      /* Number of columns in output image. */
    int NY = (int)(NRM->sz[2]);      /* Number of rows in output image. */
    int x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { /* Compute average normal inside pixel in column {x}, row {y}: */
            r3_t nrm = pst_normal_map_pixel_avg(nrmf, x, y, NS, xym_to_uvm, uvw_to_xyz);
            /* Store in normal map: */
            pst_normal_map_set_pixel(NRM, x, y, &nrm);
          }
      }
  }
  
r3_t pst_normal_map_pixel_avg
  ( pst_normal_map_proc_t nrmf, /* Normal-computing function. */
    int x, int y,               /* Pixel indices (coords of lower left corner). */
    int NS,                     /* Order of sub-sampling grid in pixel. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {xy} coords to model {uv} coords. */
    r3x3_t *uvw_to_xyz          /* Linear map of {uvw} coords to normal {xyz} coords. */
  )
  { double step = 1.0/((double)NS);
    r2_t grd_uv_sum = (r2_t){{ 0, 0 }}; /* Slope totals in {U,V,W} system. */
    int sx, sy;
    for (sy = 0; sy < NS; sy++)
      { for (sx = 0; sx < NS; sx++) 
          { 
            /* Coordinates of sample point: */
            r2_t pt_xy;
            pt_xy.c[0] = x + (sx + 0.5)*step;
            pt_xy.c[1] = y + (sy + 0.5)*step;
            
            /* Convert point from {X,Y} coords to {U,V} coords: */
            r2_t pt_uv = pst_normal_map_scene_pt_from_image_pt(&pt_xy, xym_to_uvm);

            /* Compute normal in {U,V,W} system: */
            r3_t nrm_uvw = nrmf(&pt_uv);
            
            /* If any sample is undefined, the average is undefined: */
            if (r3_L_inf_norm(&nrm_uvw) == 0.0) { return (r3_t){{ 0,0,0 }}; }

            /* Convert normal direction to slopes: */
            r2_t grd_uv = pst_normal_map_slope_from_normal(&nrm_uvw, +INF);
            grd_uv_sum.c[0] += grd_uv.c[0];
            grd_uv_sum.c[1] += grd_uv.c[1];
         }
      }
      
    /* Compute average: */
    r2_scale(1.0/(NS*NS), &grd_uv_sum, &grd_uv_sum);
    
    /* Convert gradient to normal: */
    r3_t nrm_uvw = pst_normal_map_normal_from_slope(&grd_uv_sum);

    /* Convert normal direction back to image {X,Y,Z} coords. */
    r3_t nrm_xyz =  pst_normal_map_image_dir_from_scene_dir(&nrm_uvw, uvw_to_xyz);

    return nrm_xyz;
  }
  
void pst_normal_map_perturb(float_image_t *NRM, double noise)
  { assert(NRM->sz[0] == 3);
    int NX = (int)(NRM->sz[1]);
    int NY = (int)(NRM->sz[2]);
    int x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            /* Perturb normal if so requested: */
            if ((noise > 0.0) && (r3_L_inf_norm(&nrm) > 0.0)) 
              { pst_perturb_normal(&nrm, noise);
                pst_normal_map_set_pixel(NRM, x, y, &nrm);
              }
          }
      }
  }

void pst_perturb_normal(r3_t *nrm, double amt)
  { /* Generate a Gaussian random vector with mean 0, deviation {1} in each coordinate: */
    r3_t rnd; 
    r3_throw_normal(&rnd);
    /* Scale {rnd} by {1/sqrt(2)} so that a small {amt} gives
      an rms displacement of {amt}, then mix it into {nrm}
      with weights {amt}:{1-amt}: */
    double crnd = amt/M_SQRT2;
    double cnrm = 1 - amt;
    r3_mix(crnd, &rnd, cnrm, nrm, nrm);
    /* Re-normalize: */
    double len = r3_dir(nrm, nrm);
    /* If we were so unlucky as to get a near-zero mix, try again: */
    while (len < 0.00001) { r3_throw_normal(nrm); len = r3_dir(nrm, nrm); }
  }

float_image_t *pst_normal_map_to_slope_map(float_image_t *NRM, double maxSlope)
  {
    demand(NRM->sz[0] == 3, "normal map should have three channels"); 
    int NX = (int)(NRM->sz[1]); 
    int NY = (int)(NRM->sz[2]);
    float_image_t *GRD = float_image_new(2, NX, NY);
    
    int x,y;
    for (y =0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_pixel(NRM, x, y);
            r2_t grd = pst_normal_map_slope_from_normal(&nrm, maxSlope);
            pst_slope_map_set_pixel(GRD, x, y, &grd);
          }
      }
    return GRD;
  }
  
r2_t pst_normal_map_slope_from_normal(r3_t *nrm, double maxSlope)
  { double nx = nrm->c[0], ny = nrm->c[1], nz = nrm->c[2];
    double nr = hypot(nx, ny);
    double nzmin = nr / maxSlope;
    if (nz < nzmin) { nz = nzmin; }
    double dZdX = -nx / nz;
    double dZdY = -ny / nz;
    return (r2_t) {{ dZdX, dZdY }};
  }

float_image_t *pst_normal_map_from_slope_map(float_image_t *GRD)
  {
    demand(GRD->sz[0] == 2, "slope map should have two channels"); 
    int NX = (int)(GRD->sz[1]); 
    int NY = (int)(GRD->sz[2]);
    float_image_t *NRM = float_image_new(3, NX, NY);
    
    int x,y;
    for (y =0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r2_t grd = pst_slope_map_get_pixel(GRD, x, y);
            r3_t nrm = pst_normal_map_normal_from_slope(&grd);
            pst_normal_map_set_pixel(NRM, x, y, &nrm);
          }
      }
    return NRM;
  }

r3_t pst_normal_map_normal_from_slope(r2_t *grd)
  { double dZdX = grd->c[0]; 
    double dZdY = grd->c[1]; 
    double m = sqrt(1.0 + dZdX*dZdX + dZdY*dZdY);
    double nx = dZdX/m; 
    double ny = dZdY/m; 
    double nz = 1.0/m; 
    return (r3_t){{ nx, ny, nz }};
  }
