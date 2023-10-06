#define PROG_NAME "test_image_window_op"
#define PROG_DESC "test of {image_window_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-09-24 00:12:13 by stolfi */ 
/* Created on 2012-01-25 by J. Stolfi, UNICAMP */

#define test_image_window_op_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <bool.h>
#include <r2.h>
#include <frgb_ops.h>
#include <frgb.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_test.h>
#include <image_window_op.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

typedef double tiwo_func_t (int32_t x, int32_t y);
  /* Type of a function that computes a window sample value given the
    indices {x,y} in {-1,00,+1}. */

int32_t main(int32_t argn, char **argv);

#define NUM_TEST_IMAGES 29
/* Number of distinct test images. */

void tiwo_fill_test_image(int32_t it, float_image_t *timg);
  /* Fills {timg} with test image {it}. */

void tiwo_get_function_and_mapping(int32_t it, int32_t *ipP, int32_t *igP);
  /* Returns the basic function {*ipP} and the domain mapping {*igP} 
    corresponding to test image {it}. */
   
double tiwo_compute_test_function(double x, double y, int32_t ip);
  /* Evaluates test funtion {ip} at a point with relative image coordinates {x,y}.
    Assumes that origin {(0,0)} is the center of the image. */

void tiwo_apply_domain_mapping(double X, double Y, int32_t ig, double *xP, double *yP);
  /* Converts relative image coordinates {X,Y} to function arguments {*xP,*yP}
    according to the domain mapping with index {ig}. */

void tiwo_test_operator(image_window_op_t op, int32_t NX, int32_t NY);
  /* Applies the operator to all input functions and geometric mappings,
    converted to {NX} by [NY} images, pastes the test and
    result images (plain and squared) into {oimg}. */
  
void tiwo_do_single_test
  ( float_image_t *timg,
    image_window_op_t op,
    bool_t smoothed,
    bool_t squared,
    float_image_t *rimg
  );

void tiwo_paste_image
  ( float_image_t *iimg, 
    double lo, 
    double hi, 
    int32_t xs, 
    int32_t ys, 
    bool_t sgn, 
    float_image_t *oimg
  );
  /* Pastes into color image {oimg} the greyscale test image {iimg}, with lower left pixel
    at position {(x,y)}. Assumes that the samples of {iimg} range in {[lo_hi]}.
    
    If {sgn} is false, the image {iimg} is converted to gray, mapping {lo}
    to bçack and {hi} to white.
    
    Ifs {sgn} is true, the image {iimg} is converted to scale, mapping negative
    values in {[lo _ 0)} from black to blue, 0 to gray, and
    positive values in {(0 _ hi]} from red to white.*/

void tiwo_paste_image_signed(float_image_t *iimg, int32_t c, int32_t x, int32_t y, float_image_t *oimg);
  /* Pastes into color image {oimg} the greyscale operator result image {iimg}, with lower left pixel
    at position {(x,y)}. The image is converted to gray by color mapping that
    preserves brightness but uses blues for negatives, red for positives. */

void tiwo_write_image(const char *name, float_image_t *img);
  /* Writes image {img} to file "out/{io}-{name}.pgm". */

void tiwo_check_op_range(image_window_op_t op);
  /* Estimates the actual range of all versions of operator {op}
    (smoothed or not, squared or not) by evaluating it on 
    all possible windows with sample values spanning {[0_1]}.
    Compares it with the range returned by {image_window_op_get_range}. */

/* CHECKING VALUES ON SPECIFIC FUNCTIONS */ 

void tiwo_test_unsmoothed_diff_ops(void);
  /* Tests the unsmoothed operators {f,fx,fy,fxx,fxy,fxxy,fxyy,fxxyy}
    on windows generated from the quadtatic polynomial {P}. */
  
void tiwo_test_smoothed_diff_ops(void);
  /* Tests the smoothed operators {f,fx,fy,fxx,fxy,fxxy,fxyy,fxxyy}
    on windows generated from the quadratic polynomial {Q}. */
 
void tiwo_test_diff_ops_on_quadratic
  ( tiwo_func_t W,
    char *Wname,
    bool_t smoothed,
    bool_t squared,
    double f_exp,
    double fx_exp,
    double fy_exp,
    double fxx_exp,
    double fxy_exp,
    double fyy_exp,
    double fxxy_exp,
    double fxyy_exp,
    double fxxyy_exp
  );
  /* Calls {image_window_op_apply} to evaluate the nine differential operators on a window
    whose samples are produced by the function {W}. Checks whether
    their values match the given expected values {f_exp}, {fx_exp}, etc.
    Also tests the operators {laplacian,orthicity,elongation}.
    
    The {squared} and {smoothed} flags are passed to {image_window_op_apply}.
    If {squared} is true, compares the results with {f_exp^2}, {fx_exp^2}, etc. */
  
void tiwo_test_linfit_ops(void);
  /* Tests the operators {linf,linfx,linfy,linvar,lindev}
    on windows generated from an affine function {L}
    plus noise that is orthogonal to it. */
 
void tiwo_test_linfit_ops_on_quadratic
  ( tiwo_func_t W,
    char *Wname,
    bool_t squared,
    double linf_exp,
    double linfx_exp,
    double linfy_exp,
    double linvar_exp
  );
  /* Calls {image_window_op_apply} to evaluate the linear fit operators 
    {linf,linfx,linfy,linvar,lindev} on a window
    whose samples are produced by the function {W}. Checks whether
    their values match the given expected values {linf_exp}, {linfx_exp}, etc.
    
    The {squared} flags are passed to {image_window_op_apply}.  Assumes that the 
    operator {linf} is {average}, and {linfx,linfy} are the smoothed {fx,fy}.
    If {squared} is true, compares the results with {linf_exp^2}, {linfx_exp^2}, etc. */

void tiwo_check_op(char *name, bool_t smoothed, bool_t squared, double v_cmp, double v_exp);
  /* Checks whether the computed value {v_cmp} of an operator
    matches its expected value {v_exp} (or {v_exp^2} if {squared} is true.
    The parameters {name} and {smoothed} are used only in the error message
    if the check fails. */

void tiwo_dump_window_samples(int32_t nwx, int32_t nwy, double wsmp[]);
  /* Prints the window samples {wsmp[0..nwx*nwy-1]} to {stderr}. */
  
/* IMPLEMENTATIONS */

int32_t main (int32_t argn, char **argv)
  {
    /* Size of individual test images: */
    int32_t NX = 11, NY = 11;
    
    for (int32_t ik = 0; ik < image_window_op_NUM; ik++)
      { image_window_op_t op = (image_window_op_t)ik;
        tiwo_test_operator(op, NX, NY); 
      }
      
    tiwo_test_unsmoothed_diff_ops();
    tiwo_test_smoothed_diff_ops();
    tiwo_test_linfit_ops();

    return 0;
  }

void tiwo_test_operator(image_window_op_t op, int32_t NX, int32_t NY)
  {
    const char *oname = image_window_op_to_string(op);
    fprintf(stderr, "=== testing operator %s ============\n", oname);

    tiwo_check_op_range(op);
    
    int32_t nt =  NUM_TEST_IMAGES; /* Number of test images. */

    /* Output images of one operator: */
    int32_t NCO = 3;        /* Color-coded output. */
    int32_t NXO = nt*(NX+1);    /* Test images are placed side by side. */
    int32_t NYO = 5*NY;     /* From top: test image, op, op squared, op smoothed, op smoothed squared. */
    float_image_t *oimg = float_image_new(NCO, NXO, NYO);

    /* Test and result images: */
    int32_t NC = 1;
    float_image_t *timg = float_image_new(NC, NX, NY);
    float_image_t *rimg = float_image_new(NC, NX, NY);
    
    for (int32_t it = 0; it < nt; it++)
      { fprintf(stderr, "--- testing with test image t%02d ---\n", it);
        int32_t xs = it*(NX+1);  /* Column in {oimg} where sub-images should go: */
        int32_t ys = NYO;  /* Row in {oimg} where sub-images should go: */
        double lo, hi;

        tiwo_fill_test_image(it, timg);
        lo = 0; hi = 1;
        ys = ys - NY;
        tiwo_paste_image(timg, lo, hi, xs, ys, FALSE, oimg);
        for (int32_t im = 0; im < 2; im++)
          { bool_t smoothed = (bool_t)im;
            for (int32_t iq = 0; iq < 2; iq++)
              { bool_t squared = (bool_t)iq;
                tiwo_do_single_test(timg, op, smoothed, squared, rimg);
                image_window_op_get_range(op, smoothed, squared, &lo, &hi);
                ys = ys - NY;
                tiwo_paste_image(rimg, lo, hi, xs, ys, (lo < 0), oimg);
              }
          }
      }
      
    /* Write the test and result images: */
    tiwo_write_image(oname, oimg);
    
    float_image_free(timg);
    float_image_free(rimg);
    float_image_free(oimg);

    fprintf(stderr, "=============================================================\n");
  }

void tiwo_check_op_range(image_window_op_t op)
  { 
    for (int32_t im = 0; im < 2; im++)
      { bool_t smoothed = (bool_t)im;
        for (int32_t iq = 0; iq < 2; iq++)
          { bool_t squared = (bool_t)iq;
          
            fprintf(stderr, "checking range for smoothed = %c squared = %c\n", "FT"[smoothed], "FT"[squared]);

            double vlo_exp, vhi_exp;
            image_window_op_get_range(op, smoothed, squared, &vlo_exp, &vhi_exp);
            fprintf(stderr, "expected range = [%.12f _ %.12f]\n", vlo_exp, vhi_exp);
    
            /* Allocate the test window: */
            int32_t nwx = 3, nwy = 3;
            int32_t nsmp = nwx*nwy;
            double wsmp[nsmp]; /* The window samples, linearized by rows. */
            int32_t iwctr = nsmp/2; /* Index of center sample in window. */
            
            double vlo_cmp = +INF, vhi_cmp = -INF; /* Actual observed operator range. */
            int32_t q = 3; /* Number of levels per sample. */

            /* Generate all possible 3x3 windows with {q} values per sample. */
            int64_t nw = ipow(q, nsmp); /* Number of test windows to generate. */
            for (int64_t iw = 0; iw < nw; iw++)
              { /* Split {iw} into {nsmp} base-{q} digits, map each to {0_1]: */
                int64_t b = iw;
                for (int32_t y = -1; y <= +1; y++)
                  { for (int32_t x = -1; x <= +1; x++)
                      { double val = ((double)(b % q))/((double)(q-1));
                        wsmp[iwctr + x + nwx*y] = val;
                        b = b/q;
                      }
                  }
                /* Evaluate the operator on that window: */
                double val = image_window_op_apply(op, smoothed, squared, iwctr, nwx, wsmp);
                
                if ((val < vlo_exp) || (val > vhi_exp))
                  { fprintf(stderr, "** op result %.12f outside stated range\n", val);
                    tiwo_dump_window_samples(nwx, nwy, wsmp);
                    assert(FALSE);
                  }
                if (val < vlo_cmp) { vlo_cmp = val; }
                if (val > vhi_cmp) { vhi_cmp = val; }
              }
            if (fabs(vlo_cmp - vlo_exp) > 1.0e-12*(fabs(vlo_cmp) + fabs(vlo_exp)))
              { fprintf(stderr, "** lower bound error exp = %.12e cmp = %.12e\n", vlo_exp, vlo_cmp);
                assert(FALSE);
              }
            if (fabs(vhi_cmp - vhi_exp) > 1.0e-12*(fabs(vhi_cmp) + fabs(vhi_exp)))
              { fprintf(stderr, "** upper bound error exp = %.12e cmp = %.12e\n", vhi_exp, vhi_cmp);
                assert(FALSE);
              }
          }
      }
  }

void tiwo_do_single_test
  ( float_image_t *timg,
    image_window_op_t op,
    bool_t smoothed,
    bool_t squared,
    float_image_t *rimg
  )
  {
    fprintf(stderr, "applying operator %s", image_window_op_to_string(op));
    fprintf(stderr, " smoothed = %c", "FT"[smoothed]);
    fprintf(stderr, " squared = %c\n", "FT"[squared]);
   
    int32_t NC = (int32_t)timg->sz[0]; assert(NC == rimg->sz[0]);
    int32_t NX = (int32_t)timg->sz[1]; assert(NX == rimg->sz[1]);
    int32_t NY = (int32_t)timg->sz[2]; assert(NY == rimg->sz[2]);
    assert(NC == 1);
    
    int32_t nwx, nwy, iwctr;
    image_window_op_get_window_size(op, smoothed, &nwx, &nwy, &iwctr);
    int32_t nsmp = nwx*nwy;
    float fsmp[nsmp];
    double wsmp[nsmp];
        
    for (int32_t ix = 0; ix < NX; ix++)
      { for (int32_t iy = 0; iy < NY; iy++) 
          { float_image_get_window_samples(timg, 0,ix,iy, nwx, nwy, FALSE, fsmp);
            for (int32_t k = 0; k < nsmp; k++) { wsmp[k] = fsmp[k]; }
            double res = image_window_op_apply(op, smoothed, squared, iwctr, nwx, wsmp);
            float_image_set_sample(rimg, 0,ix,iy, (float)res);
          }
      }
  }  
  

void tiwo_fill_test_image(int32_t it, float_image_t *timg)
  {
    auto void gen_proc(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
      /* Computes the value of test image at point {p} of domain. */
      
    /* Get the test function and mapping: */
    int32_t ip, ig;
    tiwo_get_function_and_mapping(it, &ip, &ig);

    fprintf(stderr, "--- testing with function f%02d mapping g%02d ---\n", ip, ig);
    fprintf(stderr, "filling test image...\n");

    /* Create the isolated test image: */
    int32_t NC = (int32_t)timg->sz[0];
    assert(NC == 1);
    float_image_test_paint(timg, gen_proc, 4);
    
    /* Rescale values to span the range {[0_1]}: */
    float vMin = +INF, vMax = -INF;
    float_image_update_sample_range(timg, 0, &vMin, &vMax);
    vMin += -1.0e-38f; /* To void {NAN} values if {img} is all zeros. */
    vMax += +1.0e-38f; /* To void {NAN} values if {img} is all zeros. */
    float_image_rescale_samples(timg, 0, vMin,vMax, 0.0,1.0);

    void gen_proc(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[])
      { 
        assert(NC == NC);
        assert(NX == timg->sz[1]);
        assert(NY == timg->sz[2]);
        
        double X = p->c[0]/NX - 0.5; /* Pixel center X, in [-1 _ +1]. */
        double Y = p->c[1]/NY - 0.5; /* Pixel center Y, in [-1 _ +1]. */

        /* Perform coordinate conversion from {X,Y} to {x,y}: */
        double x, y;
        tiwo_apply_domain_mapping(X, Y, ig, &x, &y);

        /* Compute the image as function {ip} of {x,y}: */
        fs[0] = (float)tiwo_compute_test_function(x, y, ip);

      }
  }
  
void tiwo_get_function_and_mapping(int32_t it, int32_t *ipP, int32_t *igP)  
  {
    switch(it)
      {
        case  0: (*ipP) = 0; (*igP) = 0; break; /* Constant. */
        
        case  1: (*ipP) = 1; (*igP) = 0; break; /* Linear ramp. */
        case  2: (*ipP) = 1; (*igP) = 1; break; /* Linear ramp. */
        case  3: (*ipP) = 1; (*igP) = 2; break; /* Linear ramp. */
        case  4: (*ipP) = 1; (*igP) = 3; break; /* Linear ramp. */
        case  5: (*ipP) = 1; (*igP) = 4; break; /* Linear ramp. */
        
        case  6: (*ipP) = 2; (*igP) = 0; break; /* Concave hump. */
        case  7: (*ipP) = 2; (*igP) = 4; break; /* Concave hump. */
        
        case  8: (*ipP) = 3; (*igP) = 0; break; /* Convex hump. */
        case  9: (*ipP) = 3; (*igP) = 4; break; /* Convex hump. */
        
        case 10: (*ipP) = 4; (*igP) = 0; break; /* Concave cylinder. */
        case 11: (*ipP) = 4; (*igP) = 1; break; /* Concave cylinder. */
        case 12: (*ipP) = 4; (*igP) = 3; break; /* Concave cylinder. */
        case 13: (*ipP) = 4; (*igP) = 4; break; /* Concave cylinder. */

        case 14: (*ipP) = 5; (*igP) = 0; break; /* Convex cylinder. */
        case 15: (*ipP) = 5; (*igP) = 1; break; /* Convex cylinder. */
        case 16: (*ipP) = 5; (*igP) = 3; break; /* Convex cylinder. */
        case 17: (*ipP) = 5; (*igP) = 4; break; /* Convex cylinder. */

        case 18: (*ipP) = 6; (*igP) = 0; break; /* Saddle. */ 
        case 19: (*ipP) = 6; (*igP) = 1; break; /* Saddle. */ 
        case 20: (*ipP) = 6; (*igP) = 4; break; /* Saddle. */ 

        case 21: (*ipP) = 7; (*igP) = 0; break; /* Skew saddle. */ 
        case 22: (*ipP) = 7; (*igP) = 1; break; /* Skew saddle. */  
        case 23: (*ipP) = 7; (*igP) = 4; break; /* Skew saddle. */ 

        case 24: (*ipP) = 8; (*igP) = 0; break; /* Pseudorandom. */ 

        case 25: (*ipP) = 9; (*igP) = 0; break; /* Disk. */ 

        case 26: (*ipP) = 10; (*igP) = 0; break; /* Butterfly. */ 
        case 27: (*ipP) = 10; (*igP) = 3; break; /* Skew butterfly, 45 degrees. */ 
        case 28: (*ipP) = 10; (*igP) = 5; break; /* Skew butterfly, 60 degrees. */ 
        
        default: demand(FALSE, "invalid {it}");
      }
  }

void tiwo_apply_domain_mapping(double X, double Y, int32_t ig, double *xP, double *yP)
  {
    double x, y;
    switch(ig)
      { 
        case  0: /* No change: */
          x = X; y = Y; 
          break;
        case  1:  /* Rotate 90 degrees: */
          x = +Y; y = -X; 
          break;
        case  2: /* Rotate 180 degrees: */ 
          x = -X; y = -Y; 
          break;
        case  3: /* Rotate 45 degrees: */ 
          x = X+Y; y = Y-X; 
          break;
        case  4: /* Squeeze and rotate by 60 degrees: */ 
          x = 2.0*(3*X + 2*Y); y = 1.0*(2*X - 3*Y); 
          break;
        case  5: /* Rotate by 60 degrees: */ 
          x = 0.5000*X - 0.8660*Y; y = 0.8660*X + 0.5000*Y; 
          break;
        default: demand(FALSE, "invalid {ig}");
      }
    (*xP) = x;
    (*yP) = y;
  }
    
double tiwo_compute_test_function(double x, double y, int32_t ip)
  {

    switch(ip)
      { 
        case  0: /* Constant: */
          return 0.75;
        case  1:  /* Linear ramp: */
          return x; 
        case  2: /* Concave hump: */ 
          return x*x + y*y; 
        case  3: /* Convex hump: */ 
          return -(x*x + y*y);
        case  4: /* Concave cylinder: */ 
          return x*x; 
        case  5: /* Convex cylinder: */ 
          return -x*x; 
        case  6: /* Saddle: */ 
          return x*y; 
        case  7: /* Skew saddle: */ 
          return x*x-y*y; 
        case  8: /* Pseudorandom: */ 
          { x = x*1000; y = y*1000;
            return sin(32715*(x*x+1) + 44319*(y*y-2)) + cos(16777*(x*x+2) - 32111*(y*y-7));
          }
        case  9: /* Disk: */ 
          { double r2 = x*x + y*y;
            return (r2 < 0.15 ? 1.0 : 0.0);
          }
        case  10: /* Butterfly: */ 
          { return (x*y >= 0 ? 1.0 : 0.0);
          }
        default: demand(FALSE, "invalid {ip}");
      }
  }

void tiwo_paste_image
  ( float_image_t *iimg, 
    double lo, 
    double hi, 
    int32_t xs, 
    int32_t ys, 
    bool_t sgn, 
    float_image_t *oimg
  )
  {
    fprintf(stderr, "pasting image with signed = %c range = [ %+9.6f _ %+9.6f ]\n", "FT"[sgn], lo, hi);
    
    int32_t ichns = (int32_t)iimg->sz[0]; assert(ichns == 1);
    int32_t icols = (int32_t)iimg->sz[1];
    int32_t irows = (int32_t)iimg->sz[2];

    int32_t NCO = (int32_t)oimg->sz[0]; assert(NCO == 3);
    int32_t NXO = (int32_t)oimg->sz[1];
    int32_t NYO = (int32_t)oimg->sz[2];
    
    assert(ichns == 1);
    assert(NCO == 3);
    assert((xs >= 0) && (xs + icols <= NXO));
    assert((ys >= 0) && (ys + irows <= NYO));
    
    /* Check value range: */
    if (sgn) 
      { assert(lo <= 0.0); assert(hi >= 0); } 
    else
      { assert(lo == 0.0); }
    assert(lo < hi);
    double eps = 1.0e-8 * (hi - lo); 
    double tad = 0.001 * (hi - lo); 
    
    int32_t x, y;
    for (y = 0; y < irows; y++) 
      { for (x = 0; x < icols; x++) 
          { float ival = float_image_get_sample(iimg, 0,x,y);
            frgb_t opix;
            if (isnan(ival))
              { /* Substitute gray: */
                opix.c[0] = opix.c[1] = opix.c[2] = 0.5;
              }
            else if (ival < lo-eps)
              { /* Flag too low values with green: */
                opix.c[0] = 0.0000; opix.c[1] = 1.0000; opix.c[2] = 0.0000;
              }
            else if (ival > hi+eps)
              { /* Flag too high values with magenta: */
                opix.c[0] = 1.0000; opix.c[1] = 0.0000; opix.c[2] = 1.0000; }
            else
              { /* Convert from {lo_hi} to luminosity in {[0_1]}: */
                if (sgn)
                  { frgb_t bck = (frgb_t){{ 0.0000f, 0.0000f, 0.0000f }};
                    frgb_t blu = (frgb_t){{ 0.0000f, 0.6572f, 1.0000f }};
                    frgb_t red = (frgb_t){{ 1.0000f, 0.3428f, 0.0000f }};
                    frgb_t wht = (frgb_t){{ 1.0000f, 1.0000f, 1.0000f }};
                    /* Convert {[lo_hi]} to colored RGB: */
                    if (ival < -tad)
                      { double s = (ival - lo)/(0.0 - lo);
                        opix = frgb_mix(1-s, &bck, s, &blu);
                      }
                    else if (ival > +tad)
                      { double s = (ival - 0.0)/(hi - 0.0);
                        opix = frgb_mix(1-s, &red, s, &wht);
                      }
                    else
                      { double s = (ival + tad)/(2*tad);
                        opix = frgb_mix(1-s, &blu, s, &red);
                      }
                  }
                else
                  { /* Convert from {lo_hi} to luminosity in {[0_1]}: */
                    double lum = (ival - lo)/(hi - lo); 
                    if (lum < 0) { lum = 0; }
                    if (lum > 1) { lum = 1; }
                    opix.c[0] = opix.c[1] = opix.c[2] = (float)lum;
                  }
              }
            float_image_set_pixel(oimg, xs+x, ys+y, opix.c);
          }
      }
  }

void tiwo_write_image(const char *name, float_image_t *img)
  {
    int32_t NC = (int32_t)img->sz[0];
    assert(NC == 3);
    char *fname = NULL;
    asprintf(&fname, "out/t-%s.ppm", name);
    FILE *wr = open_write(fname, TRUE);
    bool_t yup = TRUE;
    bool_t verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(img, isMask, NC, NULL, NULL, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }

void tiwo_test_unsmoothed_diff_ops(void)
  { 
    double A = dabrandom(-1.0, +1.0);
    double B = dabrandom(-1.0, +1.0);
    double C = dabrandom(-1.0, +1.0);
    double D = dabrandom(-1.0, +1.0);
    double E = dabrandom(-1.0, +1.0);
    double F = dabrandom(-1.0, +1.0);
    double G = dabrandom(-1.0, +1.0);
    double H = dabrandom(-1.0, +1.0);
    double I = dabrandom(-1.0, +1.0);
    
    auto double P(int32_t x, int32_t y);
      /* The quadratic polynomial assumed by the UNSMOOTHED differential operators. */
    
    tiwo_test_diff_ops_on_quadratic(P, "P", FALSE, FALSE, A, B, C, 2*D, E, 2*F, 2*G, 2*H, 4*I);
    tiwo_test_diff_ops_on_quadratic(P, "P", FALSE, TRUE,  A, B, C, 2*D, E, 2*F, 2*G, 2*H, 4*I);
    
    return;
    
    double P(int32_t x, int32_t y)
      { return 
          A + 
          B*x + C*y + 
          D*x*x + E*x*y + F*y*y + 
          G*x*x*y + H*x*y*y + 
          I*x*x*y*y;
      }
  }     
  
void tiwo_test_smoothed_diff_ops(void)
  {
    double A = dabrandom(-1.0, +1.0);
    double B = dabrandom(-1.0, +1.0);
    double C = dabrandom(-1.0, +1.0);
    double D = dabrandom(-1.0, +1.0);
    double E = dabrandom(-1.0, +1.0);
    double F = dabrandom(-1.0, +1.0);
    double G = dabrandom(-1.0, +1.0);
    double H = dabrandom(-1.0, +1.0);
    double I = dabrandom(-1.0, +1.0);
    
    auto double Q(int32_t x, int32_t y);
      /* The quadratic polynomial assumed by the SMOOTHED differential operators. */
    
    tiwo_test_diff_ops_on_quadratic(Q, "Q", TRUE, FALSE, A-D-F+I, B, C, 4*D, E, 4*F, 4*G, 4*H, 16*I);
    tiwo_test_diff_ops_on_quadratic(Q, "Q", TRUE, TRUE,  A-D-F+I, B, C, 4*D, E, 4*F, 4*G, 4*H, 16*I);
    
    return;
    
    double Q(int32_t x, int32_t y)
      { return 
          A + 
          B*x + C*y + 
          D*(2*x*x-1) + E*x*y + F*(2*y*y-1) + 
          G*(2*x*x-y*y)*y + H*(2*y*y-x*x)*x + 
          I*(1 - 2*(x*x-y*y)*(x*x-y*y));
      }
  }  

void tiwo_test_diff_ops_on_quadratic
  ( tiwo_func_t W,
    char *Wname,
    bool_t smoothed,
    bool_t squared,
    double f_exp,
    double fx_exp,
    double fy_exp,
    double fxx_exp,
    double fxy_exp,
    double fyy_exp,
    double fxxy_exp,
    double fxyy_exp,
    double fxxyy_exp
  )
  {
    fprintf(stderr, "=== quadratic %s smoothed = %c squared = %c ===\n", Wname, "FT"[smoothed], "FT"[squared]);
    
    /* Fill the window: */
    int nwx=3, nwy = 3;
    double wsmp[nwx*nwy];
    int32_t iwctr = nwx*nwy/2; /* Index of center sample. */
    for (int32_t y = -1; y <= +1; y++) 
      { for (int32_t x = -1; x <= +1; x++) 
          { wsmp[iwctr + y*nwx + x] = W(x,y); }
      }
      
    /* Compute the differential operators: */
    double f_cmp =     image_window_op_apply(image_window_op_F,     smoothed, squared, iwctr, nwx, wsmp);
    double fx_cmp =    image_window_op_apply(image_window_op_FX,    smoothed, squared, iwctr, nwx, wsmp);
    double fy_cmp =    image_window_op_apply(image_window_op_FY,    smoothed, squared, iwctr, nwx, wsmp);
    double fxx_cmp =   image_window_op_apply(image_window_op_FXX,   smoothed, squared, iwctr, nwx, wsmp);
    double fxy_cmp =   image_window_op_apply(image_window_op_FXY,   smoothed, squared, iwctr, nwx, wsmp);
    double fyy_cmp =   image_window_op_apply(image_window_op_FYY,   smoothed, squared, iwctr, nwx, wsmp);
    double fxxy_cmp =  image_window_op_apply(image_window_op_FXXY,  smoothed, squared, iwctr, nwx, wsmp);
    double fxyy_cmp =  image_window_op_apply(image_window_op_FXYY,  smoothed, squared, iwctr, nwx, wsmp);
    double fxxyy_cmp = image_window_op_apply(image_window_op_FXXYY, smoothed, squared, iwctr, nwx, wsmp);
    
    double laplacian_exp =  fxx_exp + fyy_exp;                  
    double orthicity_exp =  fxx_exp - fyy_exp;                  
    double elongation_exp = hypot(2*fxy_exp, fxx_exp - fyy_exp);
    
    double laplacian_cmp =  image_window_op_apply(image_window_op_LAPLACIAN,  smoothed, squared, iwctr, nwx, wsmp);
    double orthicity_cmp =  image_window_op_apply(image_window_op_ORTHICITY,  smoothed, squared, iwctr, nwx, wsmp);
    double elongation_cmp = image_window_op_apply(image_window_op_ELONGATION, smoothed, squared, iwctr, nwx, wsmp);
    
    tiwo_check_op("F",     smoothed, squared, f_cmp,     f_exp    );
    tiwo_check_op("FX",    smoothed, squared, fx_cmp,    fx_exp   );
    tiwo_check_op("FY",    smoothed, squared, fy_cmp,    fy_exp   );
    tiwo_check_op("FXX",   smoothed, squared, fxx_cmp,   fxx_exp  );
    tiwo_check_op("FXY",   smoothed, squared, fxy_cmp,   fxy_exp  );
    tiwo_check_op("FYY",   smoothed, squared, fyy_cmp,   fyy_exp  );
    tiwo_check_op("FXXY",  smoothed, squared, fxxy_cmp,  fxxy_exp );
    tiwo_check_op("FXYY",  smoothed, squared, fxyy_cmp,  fxyy_exp );
    tiwo_check_op("FXXYY", smoothed, squared, fxxyy_cmp, fxxyy_exp);
    
    tiwo_check_op("LAPLACIAN",  smoothed, squared, laplacian_cmp,  laplacian_exp);  
    tiwo_check_op("ORTHICITY",  smoothed, squared, orthicity_cmp,  orthicity_exp);  
    tiwo_check_op("ELONGATION", smoothed, squared, elongation_cmp, elongation_exp); 
  }
 
void tiwo_test_linfit_ops(void)
  {
    /* Assumes that the {linf,linfx,linfy} operators turn out to be obtained 
      by fitting the same quadratic {Q} as for the smoothed derivative
      operators, except that {linf} is {A} rather than {A-D-F+I}.
      The residual is then 
       
        { R(x,y) = D*BD(x,y) + E*BE(x,y) + F*BF(x,y) +
          G*BG(x,y) + H*BH(x,y) + I*BI(x,y) }
          
      where {BD(x,y)=2*x^2-1}, {BE(x,y)=x*y}, {BF(x,y)=2*y^2-1},
      {BG(x,y)=(2*x^2-y^2)*y}, {BH(x,y)=(2*y^2-x^2)*x},
      and {BI(x,y)=1-2*(x^2-y^2)^2}.
      
      The {linvar} operator is {<R|R>} where {<|>} denotes the inner
      product with 1/16-2/16-4/16 weights. The six basis functions
      {BD,BE,BF,BG,BH,BI} above turn out to be orthogonal under that
      inner product, with {<BD|BD>=<BF|BF>=<BI|BI>=1}, {<BE|BE>=1/4},
      {<BG|BG>=<BH|BH>=1/2}. Therefore the {linvar} operator is
      
        {linvar = D^2 + E^2/4 + F^2 + G^2/2 +H^2/2 + I^2}
        
      */
  
    double A = dabrandom(-1.0, +1.0);
    double B = dabrandom(-1.0, +1.0);
    double C = dabrandom(-1.0, +1.0);
    double D = dabrandom(-1.0, +1.0);
    double E = dabrandom(-1.0, +1.0);
    double F = dabrandom(-1.0, +1.0);
    double G = dabrandom(-1.0, +1.0);
    double H = dabrandom(-1.0, +1.0);
    double I = dabrandom(-1.0, +1.0);
    
    auto double Q(int32_t x, int32_t y);
      /* The same quadratic polynomial of the smoothed differential operators. */
    
    double linvar_exp =  D*D + E*E/4 + F*F + G*G/2 +H*H/2 + I*I;
    
    tiwo_test_linfit_ops_on_quadratic(Q, "Q", FALSE, A, B, C, linvar_exp);
    tiwo_test_linfit_ops_on_quadratic(Q, "Q", TRUE,  A, B, C, linvar_exp);
    
    return;
    
    double Q(int32_t x, int32_t y)
      { return 
          A + 
          B*x + C*y + 
          D*(2*x*x-1) + E*x*y + F*(2*y*y-1) + 
          G*(2*x*x-y*y)*y + H*(2*y*y-x*x)*x + 
          I*(1 - 2*(x*x-y*y)*(x*x-y*y));
      }
  }

void tiwo_test_linfit_ops_on_quadratic
  ( tiwo_func_t W,
    char *Wname,
    bool_t squared,
    double linf_exp,
    double linfx_exp,
    double linfy_exp,
    double linvar_exp
  )
  {
    fprintf(stderr, "=== affine %s squared = %c ===\n", Wname, "FT"[squared]);
    
    /* Fill the window: */
    int nwx=3, nwy = 3;
    double wsmp[nwx*nwy];
    int32_t iwctr = nwx*nwy/2; /* Index of center sample. */
    for (int32_t y = -1; y <= +1; y++) 
      { for (int32_t x = -1; x <= +1; x++) 
          { wsmp[iwctr + y*nwx + x] = W(x,y); }
      }
      
    /* Compute the linfit operators: */
    double linf_cmp =     image_window_op_apply(image_window_op_AVERAGE, FALSE, squared, iwctr, nwx, wsmp);
    double linfx_cmp =    image_window_op_apply(image_window_op_FX,      TRUE,  squared, iwctr, nwx, wsmp);
    double linfy_cmp =    image_window_op_apply(image_window_op_FY,      TRUE,  squared, iwctr, nwx, wsmp);
    double linvar_cmp =   image_window_op_apply(image_window_op_LINVAR,  FALSE, squared, iwctr, nwx, wsmp);
    
    double lindev_exp =  sqrt(linvar_exp);                  
    double lindev_cmp =   image_window_op_apply(image_window_op_LINDEV,  FALSE, squared, iwctr, nwx, wsmp);

    /* The {average,linvar,lindev} should not depend on {smoothed}: */
    double linf_cmp_s =   image_window_op_apply(image_window_op_AVERAGE, TRUE, squared, iwctr, nwx, wsmp);
    double linvar_cmp_s = image_window_op_apply(image_window_op_LINVAR,  TRUE, squared, iwctr, nwx, wsmp);
    double lindev_cmp_s = image_window_op_apply(image_window_op_LINDEV,  TRUE, squared, iwctr, nwx, wsmp);

    assert(linf_cmp_s == linf_cmp);
    assert(linvar_cmp_s == linvar_cmp);
    assert(lindev_cmp_s == lindev_cmp);
        
    tiwo_check_op("LINF",   FALSE, squared, linf_cmp,   linf_exp  );
    tiwo_check_op("LINFX",  FALSE, squared, linfx_cmp,  linfx_exp );
    tiwo_check_op("LINFY",  FALSE, squared, linfy_cmp,  linfy_exp );
    tiwo_check_op("LINVAR", FALSE, squared, linvar_cmp, linvar_exp);
    tiwo_check_op("LINDEV", FALSE, squared, lindev_cmp, lindev_exp);
  }

void tiwo_check_op(char *name, bool_t smoothed, bool_t squared, double v_cmp, double v_exp)
  { if (squared) { v_exp = v_exp*v_exp; }
    if (fabs(v_cmp - v_exp) > 1.0e-12*(fabs(v_cmp) + fabs(v_exp)))
      { fprintf(stderr, "** operator %s", name);
        if (smoothed || squared)
          { fprintf(stderr, "(");
            if (smoothed) { fprintf(stderr, "smoothed"); }
            if (smoothed && squared) { fprintf(stderr, ","); }
            if (squared) { fprintf(stderr, "squared"); }
            fprintf(stderr, ")");
          }
        fprintf(stderr, ": expected %.10f computed %.10f\n", v_exp, v_cmp);
        assert(FALSE);
      }
  }

void tiwo_dump_window_samples(int32_t nwx, int32_t nwy, double wsmp[])
  { int32_t iwctr = (nwx*nwy)/2;
    for (int32_t y = -1; y <= +1; y++)
      { fprintf(stderr, "  [");
        for (int32_t x = -1; x <= +1; x++)
          { fprintf(stderr, " %+16.12f", wsmp[iwctr + y*nwx + x]); }
        fprintf(stderr, " ]\n");
      }
  }
            
