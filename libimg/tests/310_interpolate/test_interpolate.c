#define PROG_NAME "test_interpolate"
#define PROG_DESC "test of {float_image_interpolate.h}"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-26 16:59:22 by stolfilocal */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define test_interpolate_COPYRIGHT \
  "Copyright © 2009  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image_interpolate.h>
#include <r2.h>
#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <indexing.h>

int main(int argn, char **argv);

void do_incremental_shift_tests(float_image_t *img0, char *name);
void do_boundary_condition_tests(float_image_t *img0, char *name);
void do_plot_tests(char *name);

void do_shift_test
  ( char *name,
    char *ttag,          /* Test type tag. */
    float_image_t *iimg, /* Input image */
    r2_t ish,            /* Amount already shifted. */
    r2_t dsh,            /* Amount to shift */
    int order,           /* Continuity order. */
    ix_reduction_t red   /* Index reduction policy. */
  );
  /* Shifts image {iimg} by {dsh}, writes it
   to disk and puts it back into {iimg}. 
   Assumes that {iimg} has already been shifted by {ish}. */
  
void do_plot_test
  ( char *name,
    float_image_t *img,  /* Test image. */
    int kx,              /* Column of pulse. */
    int ky,              /* Row of pulse. */
    ix_reduction_t red   /* Index reduction policy. */
  );
  /* Fills {img} with a blip at column {ix}, row {iy}; then writes a 
    {gnuplot/splot} file with the result of interpolating  that image 
    over its domain. */

float_image_t *get_test_image(char *name, int NC);
  /* Reads a test image from "in/{name}-orig.{ext}" where {ext} is "pgm" for {NC==1},
    "ppm" for {NC==3}.  Then paints over it four angle brackets at the corners
    and a cross at the center. */ 

void write_image(char *name, char *ttag, r2_t *ish, r2_t *dsh, r2_t *osh, int order, ix_reduction_t red, float_image_t *img);
  /* Writes an image {img} that was shifted from {ish} by {dsh} to {osh}. 
    If {ish} is NULL assumes that it is the original image. */

int main (int argc, char **argv)
  {
    demand(argc == 3, "wrong num of parameters");
    char *name = argv[1];
    int NC = atoi(argv[2]);
    
    float_image_t *img0 = get_test_image(name, NC);
    write_image(name, "in", NULL, NULL, NULL, 0, FALSE, img0);

    do_incremental_shift_tests(img0, name);
    
    do_boundary_condition_tests(img0, name);
    
    float_image_free(img0);

    do_plot_tests(name);

    return 0;
  }
 
#define NSHIFTS 6
  /* Number of shift steps in {do_shift_tests}. */

void do_incremental_shift_tests(float_image_t *img0, char *name)
  {
    /* Input image for testing: */
    int NC = (int)img0->sz[0];
    int NX = (int)img0->sz[1];
    int NY = (int)img0->sz[2];

    r2_t dp[NSHIFTS] = 
      { (r2_t){{ 0.00, 0.00 }},
        (r2_t){{ 5.00, 9.00 }},
        (r2_t){{ 3.30, 2.25 }},
        (r2_t){{ 2.30, 3.50 }},
        (r2_t){{ 3.40, 3.25 }}
      };
    
    float_image_t *imgA = float_image_new(NC, NX, NY);
            
    ix_reduction_t red = ix_reduction_REPEAT;
    int order;
    for (order = -1; order <= 1; order++)
      { /* Shift little by little: */
        float_image_assign(imgA, img0);
        r2_t tsh = (r2_t){{ 0.00, 0.00 }};
        int i;
        for (i = 1; i < NSHIFTS; i++)
          { r2_t dsh = dp[i]; 
            do_shift_test(name, "is", imgA, tsh, dsh, order, red); 
            r2_add(&dsh, &tsh, &tsh);
          }

        /* Shift all at once: */
        float_image_assign(imgA, img0);
        r2_t zsh = (r2_t){{ 0.00, 0.00 }};
        do_shift_test(name, "is", imgA, zsh, tsh, order, red);
      }
    float_image_free(imgA);
  }   
 
void do_boundary_condition_tests(float_image_t *img0, char *name)
  {

    /* Input image for testing: */
    int NC = (int)img0->sz[0];
    int NX = (int)img0->sz[1];
    int NY = (int)img0->sz[2];

    r2_t dsh = (r2_t){{ 7.40, 9.25 }};
    r2_t zsh = (r2_t){{ 0.00, 0.00 }};
            
    float_image_t *imgA = float_image_new(NC, NX, NY);
            
    ix_reduction_t red;
    for (red = ix_reduction_FIRST; red <= ix_reduction_LAST; red++)
      { int order;
        for (order = -1; order <= 1; order++)
          { float_image_assign(imgA, img0);
            do_shift_test(name, "bd", imgA, zsh, dsh, order, red); 
          }
      }
    float_image_free(imgA);
  }   
 
void do_shift_test
  ( char *name,          /* Test name, for output filenames. */
    char *ttag,          /* Tag of type of test. */
    float_image_t *iimg, /* Input image. */
    r2_t ish,            /* Amount already shifted. */
    r2_t dsh,            /* Amount to shift */
    int order,           /* Continuity order. */
    ix_reduction_t red   /* Index reduction policy. */
  )
  {
    r2_t osh;
    r2_add(&ish, &dsh, &osh);
    fprintf(stderr, "shifting image from ");
    r2_print(stderr, &ish);
    fprintf(stderr, " by ");
    r2_print(stderr, &dsh);
    fprintf(stderr, " to ");
    r2_print(stderr, &osh);
    
    /* Get image dimensions: */
    int NC = (int)iimg->sz[0];
    int NX = (int)iimg->sz[1];
    int NY = (int)iimg->sz[2];
    
    float_image_t *oimg = float_image_new(NC, NX, NY);
    
    /* Allocate a pixel's worth of samples: */
    double v[NC];
    float f[NC];
    
    /* Scan pixels: */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Get coordinates of pixel center, relative to image center: */
            r2_t p = (r2_t){{ ix + 0.5, iy + 0.5 }};
            /* Get coordinates of source point for this pixel: */
            r2_t q; r2_sub(&p, &dsh, &q);
            /* Interpolate input image: */
            float_image_interpolate_pixel(iimg, q.c[0], q.c[1], order, red, v);
            /* Store into output pixel: */
            int i; for (i = 0; i < NC; i++) { f[i] = (float)v[i]; }
            float_image_set_pixel(oimg, ix, iy, f);
          }
      }
    /* Write the output image: */
    write_image(name, ttag, &ish, &dsh, &osh, order, red, oimg);
    /* Return into original image: */
    float_image_assign(iimg, oimg);
    float_image_free(oimg);
  }  

float_image_t *get_test_image(char *name, int NC)
  {
    demand((NC == 1) || (NC == 3), "bad num of channels");
    
    char *fname = NULL;
    asprintf(&fname, "in/%s-orig.%s", name, (NC == 3 ? "ppm" : "pgm"));
    bool_t isMask = FALSE; /* Assume pixels have a smooth distribution. */
    float_image_t *base = float_image_read_pnm_named(fname, isMask, 1/0.4500, 0.0327, TRUE, TRUE, FALSE);

    int BNX = (int)base->sz[1];
    int BNY = (int)base->sz[2];
    
    /* Extract a sub-image with even width and odd height: */
    int NX = BNX - (BNX % 2);
    int NY = BNY - (1 - (BNY % 2));
    float_image_t *img = float_image_crop(base, 0,NC, 0,NX, 0,NY, 0.5000);

    /* Paint corner brackets and a cross at center: */
    int c;
    for (c = 0; c < NC; c++)
      { int ip;
        int CX = NX/2, CY = NY/2;
        float val = (float)(c == NC/2 ? 1.00 : 0.33);
        for (ip = 0; ip < 10; ip++)
          { 
            float_image_set_sample(img, c, 0,     ip,      val);
            float_image_set_sample(img, c, NX-1,  ip,      val);
            float_image_set_sample(img, c, 0,     NY-1-ip, val);
            float_image_set_sample(img, c, NX-1,  NY-1-ip, val);
            
            float_image_set_sample(img, c, CX,       CY-ip,  val);
            float_image_set_sample(img, c, CX-ip,    CY,       val);
            float_image_set_sample(img, c, NX-CX,    NY-CY+ip, val);
            float_image_set_sample(img, c, NX-CX+ip, NY-CY,    val);
            
            float_image_set_sample(img, c, ip,      0,     val);
            float_image_set_sample(img, c, ip,      NY-1,  val);
            float_image_set_sample(img, c, NX-1-ip, 0,     val);
            float_image_set_sample(img, c, NX-1-ip, NY-1,  val);
          }
      }
    free(fname);
    float_image_free(base);
    return img;
  }
  
void write_image(char *name, char *ttag, r2_t *ish, r2_t *dsh, r2_t *osh, int order, ix_reduction_t red, float_image_t *img)
  {
    int NC = (int)img->sz[0];
    
    auto int shtoi(r2_t *sh, int j);
      /* Converts coordinate {j} of {sh} from {[0 _ 1]} to an integer in {0..1000}: */
    
    auto char *makename(char *ext);
    
    { char *fname = makename("fni");
      FILE *wr = open_write(fname, TRUE);
      float_image_write(wr,img);
      fclose(wr);
      free(fname);
    }
    
    { char *fname = makename((NC == 3 ? "ppm" : "pgm"));
      bool_t isMask = FALSE; /* Assume pixels have a smooth distribution. */
      float_image_write_pnm_named(fname, img, isMask, 0.4500, 0.0327, TRUE, TRUE, TRUE);
      free(fname);
    }

    int shtoi(r2_t *sh, int j)
      { return (int)floor(sh->c[j]*1000 + 0.5); }
    
    char *makename(char *ext)
      { 
        char *fname = NULL;
        if (ish == NULL)
          { asprintf(&fname, "out/%s-orig.%s", name, ext); }
        else
          { asprintf
              ( &fname,
                "out/%s-%s-C%c-R%c--%05d-%05d--%05d-%05d--%05d-%05d.%s", 
                name,
                ttag,
                "n01"[order+1],
                "SERMP"[red],
                shtoi(osh,0), shtoi(osh,1),
                shtoi(ish,0), shtoi(ish,1),
                shtoi(dsh,0), shtoi(dsh,1), ext
              );
          }
        return fname;
      }
  }

void do_plot_tests(char *name)
  {
    int NC = 1, NX = 6, NY = 7;
    int CX = NX/2, CY = NY/2;
    float_image_t *img = float_image_new(NC, NX, NY);
    ix_reduction_t red;
    for (red = ix_reduction_FIRST; red <= ix_reduction_LAST; red++)
      { do_plot_test(name, img, CX, CY, red);

        int k;
        for (k = 0; k < 2; k++)
          {
            do_plot_test(name, img, CX,     k,      red);
            do_plot_test(name, img, CX,     NY-1-k, red);
            do_plot_test(name, img, k,      CY,     red);
            do_plot_test(name, img, NX-1-k, CY,     red);

            do_plot_test(name, img, 0,      k,      red);
            do_plot_test(name, img, 0,      NY-1-k, red);
            do_plot_test(name, img, NX-1,   k,      red);
            do_plot_test(name, img, NX-1,   NY-1-k, red);

            if (k != 0)
              { do_plot_test(name, img, k,      0,      red);
                do_plot_test(name, img, NX-1-k, 0,      red);

                do_plot_test(name, img, k,      NY-1,   red);
                do_plot_test(name, img, NX-1-k, NY-1,   red);
              }
          }
      }
    float_image_free(img);
  }

void do_plot_test
  ( char *name,
    float_image_t *img,  /* Test image. */
    int kx,              /* Column of pulse. */
    int ky,              /* Row of pulse. */
    ix_reduction_t red   /* inedx reduction policy. */
  )
  {
    fprintf(stderr, "plotting blip at element [%3d,%3d]\n", kx, ky);
    
    /* Get image dimensions: */
    int NC = (int)img->sz[0]; assert(NC == 1);
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    /* Fill image with blip: */
    float_image_fill_channel(img, 0, 0.0);
    float_image_set_sample(img, 0, kx, ky, 1.0);
    
    /* Generate 2D plot: */
    char *fname = NULL;
    asprintf(&fname, "out/%s-R%c-%04d-%04d.txt", name, "SERMP"[red], kx, ky);
    FILE *plot = open_write(fname, TRUE);
    int m = 8; /* Subsampling rate */
    int G = 3; /* Extra plot margin in pixels. */
    double eps = 1.0/1024; /* Fudge to detect discontinuities. */
    int ix, iy, dx, dy;
    for (iy = -G; iy <= NY+G; iy++)
      { for (dy = 0; dy <= m; dy++)
          { for (ix = -G; ix <= NX+G; ix++)
              { for (dx = 0; dx <= m; dx++)
                  { /* Get coordinates of interpolation point: */
                    double x = ix + ((double)dx + eps)/((double)m + 2*eps);
                    double y = iy + ((double)dy + eps)/((double)m + 2*eps);
                    /* Interpolate image: */
                    float vn = (float)float_image_interpolate_sample(img, 0, x, y, -1, red);
                    float v0 = (float)float_image_interpolate_sample(img, 0, x, y,  0, red);
                    float v1 = (float)float_image_interpolate_sample(img, 0, x, y,  1, red);
                    /* Characteristic function of domain: */
                    int s = ((x >= 0) && (x <= NX) && (y >= 0) && (y <= NY));
                    /* Wrte to plot file: */
                    fprintf(plot, "%8.5f %8.5f %10.7f %10.7f %10.7f %d\n", x, y, vn, v0, v1, s);
                  }
              }
            /* Blank line between scanlines: */
            fprintf(plot, "\n");
          }
      }
    /* Cleanup: */
    fclose(plot);
    free(fname);
  }  

