/* See float_image_average.h */
/* Last edited on 2024-11-23 05:55:58 by stolfi */ 

#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <r3.h>
#include <interval.h>
#include <ix_reduce.h>
#include <float_image.h>
#include <float_image_interpolate.h>

#include <float_image_average.h>

/* INTERNAL PROTOTYPES */

void float_image_average_adjust_sampling(int nsNew, int *nsP, double *stepP, bool_t debug);
  /* Sets the number of samples {*nsP} to {nsNew}, then adjusts
    {*stepP} so that the sampling span {(*nsP)*(*stepP)}: */
    
void float_image_average_accum_pixel(int chns, double vs[], double wt, double v[], double *wtotP);
  /* Adds {wt*vs[c]} to {v[c]} for {c} in {0..chns-1}. Also adds {wt} to {*wtotP}. */

void float_image_average_fill_pixel(int chns, float v[], float und);
  /* Sets the {float} samples {v[0..chns-1]} to the value {und}. */

void float_image_average_scale_pixel(int chns, double s, double v[], float f[]);
  /* Scales {double} samples {v[0..chns-1]} by {s}, stored in {float} samples {f[0..chns-1]}. */

double float_image_average_get_sample_weight(double d, double sigma);
  /* Returns the relative weight of a sample that is displaced by 
    a distance {d} from the center sample, assuming that 
    the samples are to be smoothed by a Gaussian with 
    deviation {sigma}. */

void fiav_debug_eigen(char *label, r2_t *e, char *tail);
  /* Prints the eigenvalues {e.c[0..1]} to {stderr}. */

void fiav_debug_rows(char *label, r2x2_t *R, char *tail);
void fiav_debug_cols(char *label, r2x2_t *R, char *tail);
  /* Prints the rows (resp columns) of {R} as vectors to {stderr}. Each 
    row (resp.column) is tagged with the string {label} and its index,
    and terminated by the string {tail}. */

/* EXPORTED FUNCTIONS */

#define float_image_average_MAXSMP_PER_AXIS (1000)
  /* Max number of samples to take along each sampling axis when
    computing one output pixel. */

#define float_image_average_MAXSMP_PER_PIXEL (2000)
  /* Max total number of samples to take to compute one
    output pixel. */

void float_image_average_parallelogram
  ( float_image_t *img, 
    ix_reduce_mode_t red,
    r2_t *p,
    r2x2_t *J,
    bool_t avg,
    int order,
    float f[],
    bool_t debug
  )
  {
    /* Get image dimensions: */
    int chns = (int)img->sz[0];
    
    /* Get coordinates {(x,y)} of center point: */
    double x = p->c[0], y = p->c[1];
    
    if (debug) { fprintf(stderr, "    evaluating img(%9.4f,%9.4f)\n", x, y); }

    /* If {J} expands the space along some direction, we may need to
      interpolate at various points near {p}, spaced at most 1 unit
      apart, and combine those samples with the proper weights. To
      ensure the proper spacing, we place the samples at a unit-step
      rectangular grid aligned with the vectors {u*J} and {v*J} where
      {u} and {v} are the eigenvectors of {J*J^t}. These are the
      directions in the input domain along which {J} has the maximum
      and minimum stretch. */

    /* Display the vectors of the image that correcpond to axes of argument: */
    if (debug) { fiav_debug_cols("      J", J, "\n"); }

    /* Compute the eigenvalues {e.c[0..1]}} and eigenvectors {R.c[0..1]} of {J*J^t}: */
    r2_t e;   /* Eigenvalues of {J*J^t}. */
    r2x2_t R; /* Rows are the corresponding eigenvectors (eventually). */
    r2x2_mul_tr(J,J,&R);         /* Now {R = J*J^t}. */
    r2x2_sym_eigen(&R, &e, &R);  /* Now the rows of {R} are the eigenvectors of {J*J^t}. */

    if (debug) { fiav_debug_eigen("      e", &e, "\n"); }
    if (debug) { fiav_debug_rows("      R", &R, "\n"); }
    
    /* Compute the directions of maximum stretch: */
    r2x2_mul(&R, J, &R);        /* Now rows of {R} are the resampling axes. */
    if (debug) { fiav_debug_rows("      U", &R, "\n"); }

    /* Compute the number of steps {ns[i]} to take along each vector,
      and their weight factors. If there is no stretching along 
      the sampling axis {i}, then one sample suffices (with unit weight).
      If there is stretching by a factor {f > 1}, we take samples
      at unit distance apart, and average them, to remove the high
      frequencies. */

    double step[2]; /* Step size for smoothing of interpolated values. */
    int ns[2];          /* Will step from {-ns[i]} to {+ns[i]} along sampling axis {i}. */
    double sigma[2];    /* Amount of extra smoothing needed along sampling axis {i}. */
    int i;
    for (i = 0; i < 2; i++)
      { step[i] = 1.0; /* May be increased to avoid excessive sampling. */
        double str = hypot(R.c[i][0], R.c[i][1]); /* Stretch along sampling axis {i}. */
        if (debug) 
          { fprintf(stderr, "      str[%d] = %24.15e\n", i, str);
            double sqe = sqrt(e.c[i]);
            fprintf(stderr, "      sqe[%d] = %24.15e\n", i, sqe);
            fprintf(stderr, "      dif = %24.15e\n", str-sqe);
            affirm(fabs(str-sqe) <= 1.0e-100 + 1.0e-3*(str+sqe), "bug in eigenvalue");
          }
        
        if (str <= 1.0)
          { /* Input image is not being squeezed; just interpolate it. */
            sigma[i] = 0.0; ns[i] = 0;
            /* We don't need the sampling axis {i}: */
            R.c[i][0] = R.c[i][1] = 0.0;
          }
        else
          { /* Input image is being squeezed; must interpolate and smooth. */
            /* Normalize sampling axis vector {R[i]} to unit length: */
            R.c[i][0] /= str; R.c[i][1] /= str; 
            /* Compute width {sigma[i]} of Gaussian smoothing kernel: */
            sigma[i] = sqrt(str*str - 1.0);
            /* Compute number of steps needed for {sigma[i]}-smoothing: */
            ns[i] = (int)ceil(3*sigma[i]/step[i]) - 1;
            if (ns[i] < 0) { ns[i] = 0; }
            /* Guard against excessive samples along sampling axis {i}: */
            int nsMax = float_image_average_MAXSMP_PER_AXIS;
            if (ns[i] > nsMax)
              { float_image_average_adjust_sampling(nsMax, &(ns[i]), &(step[i]), debug); }
          }
      }
      
    /* Guard against excessive total samples: */
    int ts = (2*ns[0] + 1)*(2*ns[1] + 1);
    int tsMax = float_image_average_MAXSMP_PER_PIXEL;
    if (ts > tsMax)
      { double sMean = sqrt(((double)tsMax)/((double)ts));  /* {ns} shrink factor. */
        for (i = 0; i < 2; i++)
          { int nsNew = (int)ceil(sMean*ns[i]);
            float_image_average_adjust_sampling(nsNew, &(ns[i]), &(step[i]), debug);
          }
      }
      
    if (debug) 
      { /* Show sampling parameters: */
        for (i = 0; i < 2; i++)
          { fprintf(stderr, "      sigma[%d] = %8.2f", i, sigma[i]);
            fprintf(stderr, " ns[%d] = %d", i, ns[i]);
            fprintf(stderr, " step[%d] = %10.6f\n", i, step[i]);
          }
      }

    /* Temporary result pixel: */
    double v[chns];

    /* Clear the result pixel: */
    int ich;
    for (ich = 0; ich < chns; ich++) { v[ich] = 0.0; }

    /* Accumulate samples: */
    double vs[chns];  /* Pixel at one sample point. */
    double wtot = 1.0e-200;
    int k0, k1;
    for (k0 = -ns[0]; k0 <= +ns[0]; k0++)
      { double f0 = k0*step[0]; /* Sample displacement along evec {R[0]}. */
        double w0 = float_image_average_get_sample_weight(f0, sigma[0]); /* Sample weight for {f0}. */
        for (k1 = -ns[1]; k1 <= +ns[1]; k1++)
          { double f1 = k1*step[1]; /* Sample displacement along evec {R[1]}. */
            double w1 = float_image_average_get_sample_weight(f1, sigma[1]); /* Sample weight for {f1}. */
            /* Compute sample point {(xs,ys)}: */
            double xs = x + f0*R.c[0][0] + f1*R.c[1][0];
            double ys = y + f0*R.c[0][1] + f1*R.c[1][1];
            /* Evaluate source image at fractional point {(xs,ys)}: */
            float_image_interpolate_pixel(img, xs, ys, order, red, vs);
            if (debug) { float_image_debug_double_pixel("        is", xs, ys, chns, vs, ""); }
            double ws = w0*w1; /* Weight of sample. */
            if (debug) { fprintf(stderr, "        weight = %8.4f\n", ws); }
            /* Accumulate: */
            float_image_average_accum_pixel(chns, vs, ws, v, &wtot);
          }
      }
    if (debug) { float_image_debug_double_pixel("      ac", x, y, chns, v, "\n"); }
    if (debug) { fprintf(stderr, "      wtot = %8.4f\n", wtot); }
    
    /* Compute scale factor to account for non-unit-sum weights: */
    assert(wtot != 0);
    double scale = 1.0/wtot;
    
    /* If {avg} is FALSE, correct for the Jacobian's stretch/shrink: */
    if (! avg) { scale *= fabs(r2x2_det(J)); }

    /* Scale pixel value: */
    float_image_average_scale_pixel(chns, scale, v, f);
    if (debug) { float_image_debug_pixel("    ev", x, y, chns, f, "\n"); }
  }

void float_image_average_persp_rectangle
  ( float_image_t *img, /* Image to sample. */
    ix_reduce_mode_t red, /* Index reduction method. */ 
    interval_t tbox[],  /* Rectangle in true coordinates. */
    r3x3_t *T2I,        /* Projective map from true coords to image coords. */
    double mrg,         /* Safety border width (pixels). */
    double avg[]        /* (OUT) Pixel average. */
  )
  { /* Get image channel count: */
    int NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    /* Get the inverse (image-to-true) map: */
    r3x3_t I2T; r3x3_inv(T2I, &I2T);
    /* Get a bounding rectangle in image coordinates: */
    interval_t ibox[2];
    r2_get_persp_rectangle_bbox(tbox, T2I, ibox);
    /* Round image rectangle to integers: */
    int xlo = (int)ceil(ibox[0].end[0]);
    int xhi = (int)floor(ibox[0].end[1]);
    int ylo = (int)ceil(ibox[1].end[0]);
    int yhi = (int)floor(ibox[1].end[1]);
    /* Scan pixels fully contained in bounding rectangle: */
    int NP = 0;       /* Count of pixels inside disk. */
    double davg[NC];  /* Sum of pixel colors inside disk. */
    int c;
    for (c = 0; c < NC; c++) { davg[c] = 0.0; }
    int x, y;
    for (y = ylo; y < yhi; y++)
      { for (x = xlo; x < xhi; x++)
          { /* Reduce indices to image domain: */
            int xr = (int)ix_reduce(x, NX, red); if (xr < 0) { continue; }
            int yr = (int)ix_reduce(y, NY, red); if (yr < 0) { continue; }
            /* Is pixel (fattened by {mrg}) inside the disk image?: */
            if (r2_pixel_is_inside_persp_rectangle(xr, yr, mrg, &I2T, tbox))
              { /* Accumulate its color in {davg,NP}: */
                float pix[NC];
                float_image_get_pixel(img, xr, yr, pix);
                for (c = 0; c < NC; c++) { davg[c] += pix[c]; }
                NP++;
              }
          }
      }
    /* Compute average color: */
    for (c = 0; c < NC; c++) { avg[c] = davg[c]/((double)NP); }
  }

void float_image_average_persp_disk
  ( float_image_t *img, /* Input image. */
    ix_reduce_mode_t red, /* Index reduction method. */ 
    r2_t *ctr,          /* Disk center in true coordinates. */
    double rad,         /* Disk radius in true coordinates. */
    r3x3_t *T2I,        /* True-to-image projective map matrix. */
    double mrg,         /* Safety border width (pixels). */
    float avg[]         /* (OUT) average disk color. */
  )
  { 
    /* Get image channel count: */
    int NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    /* Get the inverse (image-to-true) map: */
    r3x3_t I2T; r3x3_inv(T2I, &I2T);
    /* Get a bounding rectangle in image coordinates: */
    interval_t ibox[2];
    r2_get_persp_disk_bbox(ctr, rad, T2I, ibox);
    /* Round rectangle to integers: */
    int xlo = (int)ceil(ibox[0].end[0]);
    int xhi = (int)floor(ibox[0].end[1]);
    int ylo = (int)ceil(ibox[1].end[0]);
    int yhi = (int)floor(ibox[1].end[1]);
    /* Scan pixels fully contained in bounding rectangle: */
    int NP = 0;       /* Count of pixels inside disk. */
    double davg[NC];  /* Sum of pixel colors inside disk. */
    int c;
    for (c = 0; c < NC; c++) { davg[c] = 0.0; }
    int x, y;
    for (y = ylo; y < yhi; y++)
      { for (x = xlo; x < xhi; x++)
          { /* Reduce indices to image domain: */
            int xr = (int)ix_reduce(x, NX, red); if (xr < 0) { continue; }
            int yr = (int)ix_reduce(y, NY, red); if (yr < 0) { continue; }
            /* Is pixel (fattened by {mrg}) inside the disk image?: */
            if (r2_pixel_is_inside_persp_disk(xr, yr, mrg, &I2T, ctr, rad))
              { /* Accumulate its color in {davg,NP}: */
                float pix[NC];
                float_image_get_pixel(img, xr, yr, pix);
                for (c = 0; c < NC; c++) { davg[c] += pix[c]; }
                NP++;
              }
          }
      }
    /* Compute average color: */
    for (c = 0; c < NC; c++) { avg[c] = (float)(davg[c]/((double)NP)); }
  }
 
/* INTERNAL IMPLEMENTATIONS */

void float_image_average_adjust_sampling(int nsNew, int *nsP, double *stepP, bool_t debug)
  { 
    double s = ((double)nsNew)/((double)(*nsP)); /* {ns[i]} shrink factor. */
    double stepNew = (*stepP) / s;
    if (debug) 
      { fprintf(stderr, "        adjusting {ns} from %d to %d\n", (*nsP), nsNew);
        fprintf(stderr, "        adjusting {step} from %.2f to %.2f\n", (*stepP), stepNew );
      }
    (*stepP) = stepNew;
    (*nsP) = nsNew;
  }

double float_image_average_get_sample_weight(double d, double sigma)
  { 
    if ((sigma == 0.0) || (d == 0.0))
      { return 1.0; }
    else
      { double z = d/sigma;
        return exp(-z*z/2);
      }
  }

void float_image_average_accum_pixel(int chns, double vs[], double wt, double v[], double *wtotP)
  {
    if (wt != 0.0)
      { bool_t undef = FALSE;
        for (int ich = 0; ich < chns; ich++)
          { if (isnan(vs[ich])) { undef = TRUE; } }
        if (! undef)
          { for (int ich = 0; ich < chns; ich++)
              { v[ich] += wt*vs[ich]; }
            (*wtotP) += wt;
          }
      }
  }

void float_image_average_fill_pixel(int chns, float *v, float val)
  {
    int ich;
    for (ich = 0; ich < chns; ich++) 
      { v[ich] = val; }
  }

void float_image_average_scale_pixel(int chns, double s, double v[], float f[])
  {
    int ich;
    for (ich = 0; ich < chns; ich++) { f[ich] = (float)(s*v[ich]); }
  }

void fiav_debug_rows(char *label, r2x2_t *R, char *tail)
  {
    int i;
    for (i = 0; i < 2; i++)
      { fprintf(stderr, "    %s", label);
        fprintf
          ( stderr, 
            "[%d] = [%8.4f %8.4f]", 
            i, R->c[i][0], R->c[i][1]
          );
        fprintf(stderr, "%s", tail);
      }
  }

void fiav_debug_cols(char *label, r2x2_t *R, char *tail)
  {
    int j;
    for (j = 0; j < 2; j++)
      { fprintf(stderr, "    %s", label);
        fprintf
          ( stderr, 
            "[*,%d] = [%8.4f %8.4f]", 
            j, R->c[0][j], R->c[1][j]
          );
        fprintf(stderr, "%s", tail);
      }
  }

void fiav_debug_eigen(char *label, r2_t *e, char *tail)
  {
    fprintf(stderr, "    %s", label);
    fprintf(stderr, " = [%24.15e %24.15e]", e->c[0], e->c[1]);
    fprintf(stderr, "%s", tail);
  }
    
