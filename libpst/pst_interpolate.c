/* See {pst_interpolate.h}  */

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2024-12-22 12:37:27 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <affirm.h>
#include <float_image.h>
#include <pst_interpolate.h>

void pst_interpolate_two_values
  ( double v0, double w0,
    double v1, double w1,
    double *vR, double *wR
  )
  { *vR = (v0+v1)/2;
    *wR = 4.0/(1/w0 + 1/w1);
  }

void pst_interpolate_two_samples
  (  float_image_t *I, float_image_t *W,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *vR, double *wR
   )
   {
     if((I == NULL) && (W == NULL)) { *vR = 0; *wR = 1;  return;  }
     
     int32_t NX = (int32_t)(I == NULL ? W->sz[1] : I->sz[1]);
     int32_t NY = (int32_t)(I == NULL ? W->sz[2] : I->sz[2]);
     
     if((I != NULL) && (W != NULL)) 
       { demand((I->sz[1] == W->sz[1]) && (I->sz[2] == W->sz[2]), "image/mask size mismatch"); }
     if (W != NULL)
       { demand(W->sz[0] == 1, "invalid weight mask channels"); }
     
     /* Extract the sample values and weights: */
     double v0, v1;
     double w0, w1;
     
     if ((x0 >= 0) && (y0 >= 0) && (x0 < NX) && (y0 < NY))
       { v0 = (I == NULL ? 0 : float_image_get_sample(I,c,x0,y0));
         w0 = (W == NULL ? 1 : float_image_get_sample(W,0,x0,y0)); 
       }
     else
       { v0 = 0; w0 = 0; }
       
     if ((x1 >= 0) && (y1 >= 0) && (x1 < NX) && (y1 < NY))
       { v1 = (I == NULL ? 0 : float_image_get_sample(I,c,x1,y1));
         w1 = (W == NULL ? 1 : float_image_get_sample(W,0,x1,y1)); 
       }
     else
       { v1 = 0; w1 = 0; }
       
     /* Replicate the other value if weight is zero: */
     if (w0 == 0) { v0 = v1; }
     if (w1 == 0) { v1 = v0; }
       
     /* Interpolate: */
     pst_interpolate_two_values(v0, w0, v1, w1, vR, wR); 
   }

void pst_interpolate_four_values
  (  double vm, double wm,
     double v0, double w0,
     double v1, double w1,
     double vp, double wp,
     double *vR, double *wR
  )
  {
    double dist_factor = 0.25;
    
    /* Extrapolate {vm,v0 --> va}: */
    double va = (3*v0 - vm)/2.0;
    double wa = dist_factor*4.0/(9.0/w0 + (1.0/(w0*wm)));
    
    /* Interpolate {v0,v1 --> vb}: */
    double vb = (v0+v1)/2.0;
    double wb = 4.0/(1.0/w0 + 1.0/w1);
    
    /* Extrapolate {v1,vp --> vc}: */
    double vc = (3*v1 - vp)/2.0;
    double wc = dist_factor*4.0/(9.0/w1 + (1.0/(w1*wp)));
  
    /* Weighted average of the estimates: */
    if ((wa == 0) && (wb == 0) && (wc == 0))
      { *wR = 0; *vR = vb; }
    else
      { *wR = wa + wb + wc;
        *vR = (wa*va + wb*vb + wc*vc)/(wa + wb + wc);
      }
  }

void pst_interpolate_four_samples
  (  float_image_t *I, float_image_t *W,
     int32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *vR, double *wR
   )
   {
     if ((I == NULL) && (W == NULL)) { *vR = 0; *wR = 1;  return; }
     
     int32_t NX = (int32_t)(I == NULL ? W->sz[1] : I->sz[1]);
     int32_t NY = (int32_t)(I == NULL ? W->sz[2] : I->sz[2]);
     
     if((I != NULL) && (W != NULL)) 
       { demand((I->sz[1] == W->sz[1]) && (I->sz[2] == W->sz[2]), "image/mask size mismatch"); }
     if (W != NULL)
       { demand(W->sz[0] == 1, "invalid weight mask channels"); }
  
     /* Compute indices of auxiliary pixels: */
     int32_t xm = 2*x0 - x1, ym = 2*y0 - y1;
     int32_t xp = 2*x1 - x0, yp = 2*y1 - y0;
     
     /* Extract the sample values and weights: */
     double vm, v0, v1, vp;
     double wm, w0, w1, wp;
     
     if ((xm >= 0) && (ym >= 0) && (xm < NX) && (ym < NY))
       { vm = (I == NULL ? 0 : float_image_get_sample(I,c,xm,ym));
         wm = (W == NULL ? 1 : float_image_get_sample(W,0,xm,ym)); 
       }
     else
       { vm = 0; wm = 0; }
     
     if ((x0 >= 0) && (y0 >= 0) && (x0 < NX) && (y0 < NY))
       { v0 = (I == NULL ? 0 : float_image_get_sample(I,c,x0,y0));
         w0 = (W == NULL ? 1 : float_image_get_sample(W,0,x0,y0)); 
       }
     else
       { v0 = 0; w0 = 0; }
       
     if ((x1 >= 0) && (y1 >= 0) && (x1 < NX) && (y1 < NY))
       { v1 = (I == NULL ? 0 : float_image_get_sample(I,c,x1,y1));
         w1 = (W == NULL ? 1 : float_image_get_sample(W,0,x1,y1)); 
       }
     else
       { v1 = 0; w1 = 0; }
       
     if ((xp >= 0) && (yp >= 0) && (xp < NX) && (yp < NY))
       { vp = (I == NULL ? 0 : float_image_get_sample(I,c,xp,yp));
         wp = (W == NULL ? 1 : float_image_get_sample(W,0,xp,yp)); 
       }
     else
       { vp = 0; wp = 0; }
	    
     /* Interpolate: */
     pst_interpolate_four_values(vm,wm, v0,w0, v1,w1, vp,wp, vR,wR);
   }
