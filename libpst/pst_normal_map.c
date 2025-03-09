/* See pst_normal_map.h */
/* Last edited on 2025-03-04 18:36:06 by stolfi */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <jsrandom.h>
#include <affirm.h>

#include <pst_normal_map.h>
#include <pst_slope_map.h>

void pst_normal_map_pixel_avg
  ( pst_normal_func_t nrmf,   /* Normal-computing function. */
    int32_t x, int32_t y,     /* Pixel indices (coords of lower left corner). */
    int32_t NS,               /* Order of sub-sampling grid in pixel. */
    r3x3_t *img2_to_scn2,     /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img,       /* Linear map of vectors from 3D scene to 3D image coords. */
    r3_t *nrm_P,
    float *w_P
  );
  /* Computes the outwards normal vector of a surface,
    averaged over the patch that projects onto the pixel in column
    {x}, row {y} --- that is, over the square with diagonal corners
    {(x,y)} and {(x+1,y+1)} in the image coordinate system. 
    
    The average is taken over a grid of {NS×NS} samples inside the
    pixel. Normals are averaged by converting them to slopes
    {dW/dU,dW/dV} in the scene's {sx,sy,sz} system, computing the
    arithmetic mean of the slopes, and converting the result back to a
    normal vector. However,iff the computed normal is {(0,0,0)} or
    {(NAN,NAN,NAN)} for any sample point, the pixel average is
    {(NAN,NAN,NAN)}.
    
    The normal function {nrmf} is evaluated with {pst_normal_map_eval},
    using the perpective matrix {img2_to_scn2} to convert the 2D image
    coords {ix,iy} to 2D scene coords {sx,sy}, and then the linear
    matrix {scn_to_img} to convert tthe normal vector from scene's
    coords {sx,sy,sz} to image coords {ix,iy,iz}. */

/* IMPLEMENTATIONS */

r3_t pst_normal_map_get_vector(float_image_t *NRM, int32_t x, int32_t y)
  { demand((NRM->sz[0] == 3) || (NRM->sz[0] == 4), "wrong normal map depth");
    r3_t nrm;
    float *p = float_image_get_sample_address(NRM, 0, x, y);
    for (int32_t axis = 0; axis < 3; axis++)
      { nrm.c[axis] = (*p); p += NRM->st[0]; }
    return nrm;
  }
  
void pst_normal_map_set_vector(float_image_t *NRM, int32_t x, int32_t y, r3_t *nrm)
  { demand((NRM->sz[0] == 3) || (NRM->sz[0] == 4), "wrong normal map depth");
    float *p = float_image_get_sample_address(NRM, 0, x, y);
    int32_t axis;
    for (axis = 0; axis < 3; axis++)
      { (*p) = (float)nrm->c[axis]; p += NRM->st[0]; }
  }

float pst_normal_map_get_weight(float_image_t *NRM, int32_t x, int32_t y)
  { demand((NRM->sz[0] == 3) || (NRM->sz[0] == 4), "wrong normal map depth");
    float w = (NRM->sz[0] == 4 ? float_image_get_sample(NRM, 3, x, y) : 1.0f);
    return w;
  }

void pst_normal_map_set_weight(float_image_t *NRM, int32_t x, int32_t y, float w)
  { demand(NRM->sz[0] == 4, "normal map has no weight channel");
    float_image_set_sample(NRM, 3, x, y, w);
  }

r2_t pst_normal_map_scene_pt_from_image_pt(r2_t *imgp, r3x3_t *img2_to_scn2)
  { r3_t imgh = (r3_t) {{ imgp->c[0], imgp->c[1], 1.0 }};
    r3_t scnh;
    r3x3_map_row(&imgh, img2_to_scn2, &scnh);
    r2_t pscn = (r2_t){{ scnh.c[0]/scnh.c[2], scnh.c[1]/scnh.c[2] }};
    return pscn;
  }

r3_t pst_normal_map_image_dir_from_scene_dir(r3_t *scnv, r3x3_t *scn_to_img)
  { r3_t imgv;
    if (r3_L_inf_norm(scnv) == 0.0)
      { /* Normal is undefined - return null normal: */
        imgv = (r3_t) {{ 0.0, 0.0, 0.0 }};
      }
    else
      { /* Map normal from {U,V,W} to {X,Y,Z} system: */
        r3x3_map_row(scnv, scn_to_img, &imgv);
        /* Normalize: */
        (void)r3_dir(&imgv, &imgv);
      }
    return imgv;
  }

r3_t pst_normal_map_eval
  ( pst_normal_func_t nrmf,  /* Normal-computing funtion. */
    double ix,               /* X-coordinate of projected point in image system. */
    double iy,               /* Y-coordinate of projected point in image system. */
    r3x3_t *img2_to_scn2,    /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img       /* Linear mapof vectors from 3D scene to 3D image coords. */
  )
  { /* Convert point from 2D image coords to 2D scene coords: */
    r2_t imgp = (r2_t){{ ix, iy }};
    r2_t scnp = pst_normal_map_scene_pt_from_image_pt(&imgp, img2_to_scn2);

    /* Compute normal in 3D scene system: */
    r3_t scnv = nrmf(&scnp);
    
    /* Convert normal direction back to 3D image coords. */
    r3_t imgv =  pst_normal_map_image_dir_from_scene_dir(&scnv, scn_to_img);
    
    return imgv;
  }

void pst_normal_map_from_proc
  ( pst_normal_func_t nrmf,  /* Normal-computing funtion. */
    int32_t NS,              /* Order of subsampling grid within each pixel. */
    r3x3_t *img2_to_scn2,    /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img,      /* Linear map of vectors from 3D scene to 3D image coords. */
    float_image_t *NRM       /* (OUT) Computed normal map. */
  )
  { /* Get/check image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(NRM, &NC, &NX, &NY);
    demand((NC == 3) || (NC == 4), "bad normal map");

    int32_t x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { /* Compute average normal inside pixel in column {x}, row {y}: */
            r3_t nrm; float w;
            pst_normal_map_pixel_avg(nrmf, x, y, NS, img2_to_scn2, scn_to_img, &nrm, &w);
            /* Store in normal map: */
            pst_normal_map_set_vector(NRM, x, y, &nrm);
            if (NC == 4) { pst_normal_map_set_weight(NRM, x, y, w); }
          }
      }
  }
  
void pst_normal_map_pixel_avg
  ( pst_normal_func_t nrmf,   /* Normal-computing function. */
    int32_t x, int32_t y,     /* Pixel indices (coords of lower left corner). */
    int32_t NS,               /* Order of sub-sampling grid in pixel. */
    r3x3_t *img2_to_scn2,     /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img,       /* Linear map of vectors from 3D scene to 3D image coords. */
    r3_t *nrm_P,
    float *w_P
  )
  { double step = 1.0/((double)NS);
    r2_t scng_sum = (r2_t){{ 0, 0 }}; /* Slope totals in scene system. */
    float w = 1.0;
    int32_t sx, sy;
    for (sy = 0; (sy < NS) && (w > 0); sy++)
      { for (sx = 0; (sx < NS) && (w > 0); sx++) 
          { /* Coordinates of sample point: */
            r2_t imgp = (r2_t) {{ x + (sx + 0.5)*step, y + (sy + 0.5)*step }};
            
            /* Convert point from {X,Y} coords to {U,V} coords: */
            r2_t scnp = pst_normal_map_scene_pt_from_image_pt(&imgp, img2_to_scn2);

            /* Compute normal in {U,V,W} system: */
            r3_t scnv = nrmf(&scnp);
            
            /* If any sample is undefined, the average is undefined: */
            double mag = r3_L_inf_norm(&scnv);
            if (isfinite(scnv.c[0]) && isfinite(mag) && (mag > 0.0))
              { /* Convert normal direction to slopes, and accumulat: */
                r2_t scng = pst_slope_from_normal(&scnv, +INF);
                scng_sum.c[0] += scng.c[0];
                scng_sum.c[1] += scng.c[1];
              }
            else
              { w = 0.0; }
         }
      }
      
    if (w == 0.0)
      { (*nrm_P) = (r3_t){{ NAN, NAN, NAN }}; }
    else
      { /* Compute average: */
        r2_scale(1.0/(NS*NS), &scng_sum, &scng_sum);

        /* Convert gradient to normal: */
        r3_t scnv = pst_normal_from_slope(&scng_sum);

        /* Convert normal direction back to image {X,Y,Z} coords. */
        r3_t imgv =  pst_normal_map_image_dir_from_scene_dir(&scnv, scn_to_img);
        
        (*nrm_P) = imgv;
      }
    (*w_P) = w;
  }
  
void pst_normal_map_perturb(float_image_t *NRM, double noise)
  { /* Get/check image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(NRM, &NC, &NX, &NY);
    demand((NC == 3) || (NC == 4), "bad normal map");

    if (noise == 0) { return; }
    int32_t x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
            /* Perturb normal unless undefined */
            double mag = r3_L_inf_norm(&nrm);
            if (isfinite(nrm.c[0]) && isfinite(mag) && (mag > 0.0))
              { pst_perturb_normal(&nrm, noise);
                pst_normal_map_set_vector(NRM, x, y, &nrm);
              }
            else
              { if (NC == 4) { pst_normal_map_set_weight(NRM, x, y, 0); } }
          }
      }
  }

float_image_t *pst_normal_map_to_slope_map(float_image_t *NRM, double maxSlope)
  {
    /* Get/check image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(NRM, &NC, &NX, &NY);
    demand((NC == 3) || (NC == 4), "bad normal map");
    
    float_image_t *GRD = float_image_new(3, NX, NY);
    
    int32_t x,y;
    for (y =0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
            r2_t grd = pst_slope_from_normal(&nrm, maxSlope);
            pst_slope_map_set_gradient(GRD, x, y, &grd);
            float w = pst_normal_map_get_weight(NRM, x, y);
            double mag = r3_L_inf_norm(&nrm);
            if ((! isfinite(grd.c[0])) || (! isfinite(mag)) || mag == 0) { w = 0; }
            pst_slope_map_set_weight(GRD, x, y, w);
          }
      }
    return GRD;
  }

float_image_t *pst_normal_map_from_slope_map(float_image_t *GRD)
  {
    int32_t NC, NX, NY;
    float_image_get_size(GRD, &NC, &NX, &NY);
    demand((NC == 2) || (NC == 3), "bad slope map depth");

    float_image_t *NRM = float_image_new(4, NX, NY);
    int32_t x,y;
    for (y =0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { r2_t grd = pst_slope_map_get_gradient(GRD, x, y);
            r3_t nrm = pst_normal_from_slope(&grd);
            pst_normal_map_set_vector(NRM, x, y, &nrm);
            float w = pst_slope_map_get_weight(GRD, x, y);
            double mag = r3_L_inf_norm(&nrm);
            if ((! isfinite(nrm.c[0])) || (! isfinite(mag)) || mag == 0) { w = 0; }
            pst_normal_map_set_weight(NRM, x, y, w);
          }
      }
    return NRM;
  }
