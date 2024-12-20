#define PROG_NAME "test_transform"
#define PROG_DESC "test of {float_image_transform.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-20 18:19:32 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_transform_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <affirm.h>
#include <ix_reduce.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <r2x2.h>
#include <r2.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image_test.h>
#include <float_image_transform.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

int main(int argn, char **argv);

void do_test
  ( float_image_t *iimg,
    ix_reduce_mode_t red, 
    float_image_test_generator_t *gen_proc,
    char *gen_name,
    float_image_t *oimg,
    r2_map_jacobian_t *map_proc,
    char *map_name,
    bool_t avg,
    int order,
    r2_pred_t *debugp
  );

void write_image(char *gen_name, char *map_name, bool_t avg, int order, char *io, float_image_t *img);

int main (int argn, char **argv)
  {
    /* Choose a min wavelength {PMin} for test images (pixels, integer): */
    int PMin = 4;
    
    /* Input image for testing (size must be multiple of {3*Pmin}): */
    int ichns = 3, icols = 9*(3*PMin), irows = 9*(3*PMin);
    float_image_t *iimg = float_image_new(ichns, icols, irows);
    
    /* Output (transformed) image: */
    int ochns = ichns, ocols = 200, orows = 150;
    float_image_t *oimg = float_image_new(ochns, ocols, orows);

    /* Get input and output image dimensions, floated: */
    double iwx = (double)icols;
    double iwy = (double)irows;
    double owx = (double)ocols;
    double owy = (double)orows;

    /* Compute circumradii {iR,oR} of input and output images (in pixels): */
    double iR = 0.5*hypot(icols, irows);
    double oR = 0.5*hypot(ocols, orows);
    fprintf(stderr, "iR = %7.2f  oR = %7.2f\n", iR, oR);
    
    /* Compute inradii {iH,oH} of input and output images (in pixels): */
    double iH = 0.5*(icols < irows ? icols : irows);
    double oH = 0.5*(ocols < orows ? ocols : orows);
    fprintf(stderr, "iH = %7.2f  oH = %7.2f\n", iH, oH);
    
    auto bool_t debp(r2_t *p);
      /* Returns true if the output pixel with center {p} should be debugged. */
    
    auto bool_t is_close(r2_t *p, double rx, double ry);
      /* Returns true if {p} is the output pixel at about {rx,ry} of domain. */

    /* Geometry transformations: */
    auto void map_ident(r2_t *p, r2x2_t *J);
    auto void map_shift(r2_t *p, r2x2_t *J);
    auto void map_magnf(r2_t *p, r2x2_t *J);
    auto void map_reduc(r2_t *p, r2x2_t *J);
    auto void map_shear(r2_t *p, r2x2_t *J);
    auto void map_quadr(r2_t *p, r2x2_t *J);
    auto void map_polar(r2_t *p, r2x2_t *J);
    auto void map_compo(r2_t *p, r2x2_t *J);

    bool_t debp(r2_t *p)
      { /* Debug some output pixels: */
        if (is_close(p, 0.1, 0.1)) return TRUE;
        if (is_close(p, 0.1, 0.9)) return TRUE;
        if (is_close(p, 0.9, 0.1)) return TRUE;
        if (is_close(p, 0.9, 0.9)) return TRUE;
        if (is_close(p, 0.4, 0.4)) return TRUE;
        return FALSE;
      }

    bool_t is_close(r2_t *p, double rx, double ry)
      { /* Debug some output pixels: */
        r2_t cp = (r2_t){{ floor(rx*owx)+0.5, floor(ry*owy)+0.5 }}; 
        return (r2_dist(p, &cp) <= 0.5000001);
      }

    ix_reduce_mode_t red1 = ix_reduce_mode_SINGLE; /* !!! Should test other options. !!! */
    ix_reduce_mode_t redM = ix_reduce_mode_REPEAT; /* !!! Should test other options. !!! */

    int order;
    for (order = 0; order <= 1; order++)
      { /* Run some tests: */
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_ident, "1-ident", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_magnf, "2-magnf", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_shift, "3-shift", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_reduc, "4-reduc", TRUE,  order, &debp);
        do_test(iimg, redM, &float_image_test_gen_stripes, "1-stripes", oimg, &map_shear, "5-shear", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_quadr, "6-quadr", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_stripes, "1-stripes", oimg, &map_quadr, "6-quadr", FALSE, order, &debp);
        do_test(iimg, redM, &float_image_test_gen_stripes, "1-stripes", oimg, &map_polar, "7-polar", TRUE,  order, &debp);
        do_test(iimg, redM, &float_image_test_gen_stripes, "1-stripes", oimg, &map_compo, "8-compo", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_ripples, "2-ripples", oimg, &map_quadr, "6-quadr", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_checker, "3-checker", oimg, &map_quadr, "6-quadr", TRUE,  order, &debp);
        do_test(iimg, red1, &float_image_test_gen_chopsea, "4-chopsea", oimg, &map_quadr, "6-quadr", TRUE,  order, &debp);
      }
    
    return 0;
    
    /* DEFINITION OF THE GEOMETRIC MAPS */

    void map_ident(r2_t *p, r2x2_t *J)
      {
        /* The identity map (copies the image aligned at bottom left corner). */
        /* Do nothing to {p}. */
        /* The Jacobian is the identity. */
      }
  
    void map_shift(r2_t *p, r2x2_t *J)
      {
        /* Shifts the input image by fractions of a pixel. */
        
        /* Get the coordinates {(ox,oy)} of {op}, unscaled: */
        double ox = p->c[0];
        double oy = p->c[1];
        
        /* Define the input-to-output image shifts: */
        double sx = 10 + 0.50000000;
        double sy = 10 + 0.33333333;
        
        /* Apply the shifts: */
        double ix = ox - sx;
        double iy = oy - sy;
        
        /* Set the input coords: */
        p->c[0] = ix;
        p->c[1] = iy;
        
        /* The Jacobian is the identity. */
      }
  
    void map_shear(r2_t *p, r2x2_t *J)
      {
        /* Shears the image in the X direction. */
        
        /* Get the coordinates {(ox,oy)} of {op}, unscaled: */
        double ox = p->c[0];
        double oy = p->c[1];
        
        /* Define the shear coefficient: */
        double sx = -0.4;
        
        /* Magnify by inverse scale factors: */
        double ix = ox + sx*oy;
        double iy = oy;
        
        /* Set the input coords: */
        p->c[0] = ix;
        p->c[1] = iy;
        
        if (J != NULL)
          { /* Update the Jacobian: */
            r2x2_t K;
            K.c[0][0] = 1;   /* dix_dox */ 
            K.c[0][1] = 0.0; /* diy_dox */ 
            K.c[1][0] = sx;  /* dix_doy */ 
            K.c[1][1] = 1;   /* diy_doy */ 
            r2x2_mul(J, &K, J);
          }
      }

    void map_reduc(r2_t *p, r2x2_t *J)
      {
        /* Reduces the input image by fixed factor along each axis. */
        
        /* Get the coordinates {(ox,oy)} of {op}, unscaled: */
        double ox = p->c[0];
        double oy = p->c[1];
        
        /* Define the output-to-input scale factors: */
        double sx = 3.0;
        double sy = 1.1;
        
        /* Magnify by scale factors: */
        double ix = sx*ox;
        double iy = sy*oy;
        
        /* Set the input coords: */
        p->c[0] = ix;
        p->c[1] = iy;
        
        if (J != NULL)
          { /* Update the Jacobian: */
            r2x2_t K;
            K.c[0][0] = sx;  /* dix_dox */ 
            K.c[0][1] = 0.0; /* diy_dox */ 
            K.c[1][0] = 0.0; /* dix_doy */ 
            K.c[1][1] = sy;  /* diy_doy */ 
            r2x2_mul(J, &K, J);
          }
      }
  
    void map_magnf(r2_t *p, r2x2_t *J)
      {
        /* Expands the input image by fixed factor along each axis. */
        
        /* Get the coordinates {(ox,oy)} of {op}, unscaled: */
        double ox = p->c[0];
        double oy = p->c[1];
        
        /* Define the output-to-input scale factors: */
        double sx = 5.0;
        double sy = 1.1;
        
        /* Magnify by the scale factors: */
        double ix = ox/sx;
        double iy = oy/sy;
        
        /* Set the input coords: */
        p->c[0] = ix;
        p->c[1] = iy;
        
        if (J != NULL)
          { /* Update the Jacobian: */
            r2x2_t K;
            K.c[0][0] = 1/sx;  /* dix_dox */ 
            K.c[0][1] = 0.0;   /* diy_dox */ 
            K.c[1][0] = 0.0;   /* dix_doy */ 
            K.c[1][1] = 1/sy;  /* diy_doy */ 
            r2x2_mul(J, &K, J);
          }
      }
  
    void map_quadr(r2_t *p, r2x2_t *J)
      {
        /* A quadratc map with a "hole". */
        
        /* Get the coordinates {(ox,oy)} of {op}, in {[0_1]} coordinates: */
        double ox = p->c[0]/owx;
        double oy = p->c[1]/owy;
        
        /* Compute the coordinates {(ix,iy)} of {ip}, in {[0_1]} coordinates: */
        double Ax_xx = +1.00;
        double Ax_xy = 00.00;
        double Ax_yy = 00.00;
        double Bx_x =  00.00;
        double Bx_y =  00.00;
        double Cx =    00.00;
        
        double Ay_xx = 00.00;
        double Ay_xy = 00.00;
        double Ay_yy = +1.00;
        double By_x =  00.00;
        double By_y =  00.00;
        double Cy =    00.00;
        
        double ix = Ax_xx*ox*ox + 2*Ax_xy*ox*oy + Ax_yy*oy*oy + Bx_x*ox + Bx_y*oy + Cx;
        double iy = Ay_xx*ox*ox + 2*Ay_xy*ox*oy + Ay_yy*oy*oy + By_x*ox + By_y*oy + Cy;

        /* Convert input coords to pixels, with overflow: */
        double mrg = 0.1;         /* Relative overflow along each edge. */
        double mag = 1.0 + 2*mrg; /* Scaling factor. */
        
        /* Set the input coords: */
        p->c[0] = (mag*ix - mrg)*iwx ;
        p->c[1] = (mag*iy - mrg)*iwy;
        
        /* Arbitrarily define some input pixels as invalid: */
        double bx = 0.6;
        double by = 0.6;
        bool_t invalid = (hypot((ix-bx)*iwx, (iy-by)*iwy) <= 6.0);
        if (invalid)
          { (*p) = (r2_t){{ NAN, NAN }}; }
        else
          { if (J != NULL)
              { /* Update the Jacobian: */
                r2x2_t K;
                K.c[0][0] = 2*(Ax_xx*ox + Ax_xy*oy) + Bx_x; /* dix_dox */ 
                K.c[0][1] = 2*(Ay_xx*ox + Ay_xy*oy) + By_x; /* diy_dox */  
                K.c[1][0] = 2*(Ax_xy*ox + Ax_yy*oy) + Bx_y; /* dix_doy */
                K.c[1][1] = 2*(Ay_xy*ox + Ay_yy*oy) + By_y; /* diy_doy */ 

                /* Adjust {J} to account for the input and output scalings: */
                K.c[0][0] *= mag*iwx/owx;  /* diX_doX */ 
                K.c[0][1] *= mag*iwy/owx;  /* diY_doX */ 
                K.c[1][0] *= mag*iwx/owy;  /* diX_doY */ 
                K.c[1][1] *= mag*iwy/owy;  /* diY_doY */ 
                r2x2_mul(J, &K, J);
              }
          }
      }
  
    void map_polar(r2_t *p, r2x2_t *J)
      {
        /* Inverse of the quadratic polar map with pole at center of output image. */
        
        /* Get the coordinates {(ox,oy)} of {op}, in {[-1:+1]} coordinates: */
        double ox = (p->c[0] - owx/2)/oH;
        double oy = (p->c[1] - owy/2)/oH;
        
        double tiny = 0.5/oH;
        
        /* Convert to polar coordinates: */
        double or2 = ox*ox + oy*oy;
        if (or2 < tiny*tiny)
          { /* Undefined at pole: */
            (*p) = (r2_t){{ NAN, NAN }};
          }
        else
          { 
            double ot = atan2(oy, ox);

            /* Compute the coordinates {(ix,iy)} of {ip}, in {[0_1]} coordinates: */
            double ix = 0.5*(ot/M_PI + 1);
            ix = ix - floor(ix);
            double iy = or2;

            /* Convert input coords to pixels: */
            p->c[0] = ix*iwx;
            p->c[1] = iy*iwy;
            
            if (J != NULL)
              { /* Update the Jacobian: */
                r2x2_t K;
                K.c[0][0] = -oy/or2*(0.5/M_PI); /* dix_dox */
                K.c[0][1] = 2*ox;               /* diy_dox */  
                K.c[1][0] = +ox/or2*(0.5/M_PI); /* dix_doy */ 
                K.c[1][1] = 2*oy;               /* diy_doy */ 

                /* Adjust {J} to account for the input and output scalings: */
                K.c[0][0] *= iwx/oH;  /* diX_doX */  
                K.c[0][1] *= iwy/oH;  /* diY_doX */
                K.c[1][0] *= iwx/oH;  /* diX_doY */ 
                K.c[1][1] *= iwy/oH;  /* diY_doY */ 
                r2x2_mul(J, &K, J);
              }
          }
      }
      
    void map_compo(r2_t *p, r2x2_t *J)
      {
        /* Composite of a polar map and a shear map: */
        map_polar(p, J);
        map_shear(p, J);
      }
      
  }

void do_test
  ( float_image_t *iimg,
    ix_reduce_mode_t red, 
    float_image_test_generator_t *gen_proc,
    char *gen_name,
    float_image_t *oimg,
    r2_map_jacobian_t *map_proc,
    char *map_name,
    bool_t avg,
    int order,
    r2_pred_t *debugp
  )
  {
    fprintf(stderr, "=== test gen=%s map=%s avg=%c order=%d ============\n", gen_name, map_name, "FT"[avg], order);
    fprintf(stderr, "creating input image \"%s\" ...\n", gen_name);
    float_image_test_paint(iimg, gen_proc, 3);
    float def = (float)(1/M_PI); /* An arbitrary value. */
    fprintf(stderr, "filling output image with %8.6f ...\n", def);
    float_image_fill(oimg, def);
    float undef = (float)(1/M_SQRT2); /* Another arbitrary value. */
    fprintf(stderr, "applying map \"%s\", undef = %8.6f avg = %c...\n", map_name, undef, "FT"[avg]);
    float_image_transform_all(iimg, red, map_proc, undef, avg, order, debugp, oimg);
    /* Write the input and output images: */
    write_image(gen_name, map_name, avg, order, "in", iimg);
    write_image(gen_name, map_name, avg, order, "ot", oimg);
    fprintf(stderr, "=============================================================\n");
  }  

void write_image(char *gen_name, char *map_name, bool_t avg, int order, char *io, float_image_t *img)
  {
    if (map_name == NULL) { map_name = "0-nomap"; }
    char *fname = jsprintf("out/test-%s-%s-%s-%c-%d.ppm", map_name, gen_name, io, "FT"[avg], order);
    FILE *wr = open_write(fname, TRUE);
    int chns = (int)img->sz[0];
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pimg = float_image_to_uint16_image(img, isMask, chns, NULL, NULL, NULL, 255, yup, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    uint16_image_free(pimg);
    fclose(wr);
    free(fname);
  }
