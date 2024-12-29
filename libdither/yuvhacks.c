/* See {yuvhacks.h} */

/* Last edited on 2024-12-25 12:44:37 by stolfi */

#include <yuvhacks.h>

typedef struct yuv_matrix
  { /* Fractional coefficients. */
    double RY, GY, BY;
    double RGU;
    double RV, GV, BV;
    double y0;  /* Luminance offset (rel. range [0 _ 1]). */
  } yuv_matrix;
  
static yuv_matrix YUVM;

void set_yuv_matrix(void)
  { /* The RGB to YUY transformation consists of multiplying the
      RGB coordinates (in the range [0 _ 1]) by the matrix

       [ 306,  601, 117]
       [ 291, -291,   0] / 1024
       [-142,  -87, 229]

      The result ranges in [0..1024], [-291..+291], [-229..+229]
      over 1024. 
      
      Then we scale U and V by (1 + Y0)/(Y + Y0), where Y0 is
      a fixed bias.
    */
    YUVM.RY  =   306.0/1024.0;
    YUVM.GY  =   601.0/1024.0;
    YUVM.BY  =   117.0/1024.0;
    YUVM.RGU =   136.0/1024.0;
    YUVM.RV  =  -169.0/1024.0;
    YUVM.GV  =    39.0/1024.0;
    YUVM.BV  =   130.0/1024.0;
    
    YUVM.y0  =    32.0/1024.0;
  }

void rgb_to_yuv(int32_t R, int32_t G, int32_t B, uint32_t maxval, float *yp, float *up, float *vp)
  { 
    /* These quantities are scaled by {maxval}: */
    register double Y, U, V, Y0, YOFF;
    Y = YUVM.RY  * (double)R + YUVM.GY * (double)G + YUVM.BY * (double)B;
    U = YUVM.RGU * (double)(R - G);
    V = YUVM.RV  * (double)R + YUVM.GV * (double)G + YUVM.BV * (double)B;
    Y0 = YUVM.y0 * (double)maxval;
    YOFF = Y + Y0;
    if (YOFF <= Y0/2) YOFF = Y0/2;
    /* Projective correction: */
    double m = (1.0 + YUVM.y0)/YOFF;
    (*yp) = (float)(Y/((double)maxval));
    (*up) = (float)(U*m);
    (*vp) = (float)(V*m);
  }

float rgb_to_y(int32_t R, int32_t G, int32_t B, uint32_t maxval)
  { register double Y;
    Y = YUVM.RY * (double)R + YUVM.GY * (double)G + YUVM.BY * (double)B;
    return (float)(Y/((double)maxval));
  }

