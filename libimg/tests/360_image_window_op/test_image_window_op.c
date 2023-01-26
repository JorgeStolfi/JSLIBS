#define PROG_NAME "test_image_window_op"
#define PROG_DESC "test of {image_window_op.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-14 01:07:58 by stolfi */ 
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
    at position {(x,y)}. Assumes that the samples of {iimg} range in {[llo_hi]}.
    
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

int32_t main (int32_t argn, char **argv)
  {
    /* Size of individual test images: */
    int32_t NX = 11, NY = 11;
    
    for (int32_t ik = 0; ik < image_window_op_NUM_VALUES; ik++)
      { image_window_op_t op = (image_window_op_t)ik;
        tiwo_test_operator(op, NX, NY); 
      }

    return 0;
  }

void tiwo_test_operator(image_window_op_t op, int32_t NX, int32_t NY)
  {
    const char *oname = image_window_op_to_string(op);
    fprintf(stderr, "=== testing operator %s ============\n", oname);

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
                image_window_op_get_range(op, squared, &lo, &hi);
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
    
    int32_t nwx, nwy, wctr;
    image_window_op_get_window_size(op, smoothed, &nwx, &nwy, &wctr);
    int32_t nsmp = nwx*nwy;
    float fsmp[nsmp];
    double dsmp[nsmp];
        
    for (int32_t ix = 0; ix < NX; ix++)
      { for (int32_t iy = 0; iy < NY; iy++) 
          { float_image_get_window_samples(timg, 0,ix,iy, nwx, nwy, FALSE, fsmp);
            for (int32_t k = 0; k < nsmp; k++) { dsmp[k] = fsmp[k]; }
            double res = image_window_op_apply(op, smoothed, squared, wctr, nwx, dsmp);
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
          return x*x - y*y; 
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
