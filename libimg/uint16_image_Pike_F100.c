/* See uint16_image_Pike_F100.h */
/* Last edited on 2024-10-25 22:51:05 by stolfi */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_Pike_F100.h>
#include <float_image.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <sample_conv.h>

uint16_image_t *uint16_image_Pike_F100_read(char *name, bool_t verbose)
  { FILE *rd = open_read(name, verbose);
    uint16_image_t *img = uint16_image_Pike_F100_fread(rd);
    /* Read and print the trailer for diagnostics: */
    uint16_image_Pike_F100_read_trailer(rd, verbose);
    if (rd != stdin) { fclose(rd); }
    return img;
  }
  
void uint16_image_Pike_F100_read_trailer(FILE *rd, bool_t verbose)
  {
    int i;
    for (i = 0; i < 2944; i += 2)
      { int c1 = getc(rd); if (c1 == EOF) { pnm_error("unexpected EOF in sample"); }
        int c2 = getc(rd); if (c2 == EOF) { pnm_error("unexpected EOF in sample"); }
        if (verbose)
          { fprintf(stderr, "%4d ( %02x = %3d  %02x = %3d ) = %5d\n", i, c1, c1, c2, c2, (c1 << 8) | c2); }
      }
  }

uint16_image_t *uint16_image_Pike_F100_fread(FILE *rd)
  { /* The AVC Pike F-100 Raw15 format is 1000x1000, 16 bits per pixel. */
    /* Create the file header: */
    int cols = 1000;
    int rows = 1000;
    int chns = 1;
    /* Allocate image and set maxval: */
    uint16_image_t *img = uint16_image_new(cols, rows, chns);
    int32_t maxval = 65535;
    img->maxval = (uint16_t)maxval;
    /* Read pixels: */
    bool_t raw = TRUE;   /* Samples are binary, not ASCII. */
    bool_t bits = FALSE; /* Samples are byte-aligned, not bit-aligned. */
    int row;
    for (row = 0; row < rows; row++)
      { pnm_read_pixels(rd, img->smp[row], cols, chns, (uint16_t)maxval, raw, bits); }
    return(img);
  }

int uint16_image_Pike_F100_sample_channel(int col, int row)
  { int ct = (col % 2);
    int rt = (row % 2);
    int chn = (1 - ct) + rt;
    return chn;
  }
               
uint16_image_t *uint16_image_Pike_F100_debayer(uint16_image_t *img, bool_t squeeze, bool_t verbose)  
  {
    demand(img->chns == 1, "input must be single-channel");
    uint16_image_t *omg;
    if (! squeeze)
      { /* Retain the original size, just separate the channels: */
        int rows = img->rows;
        int cols = img->cols;
        omg = uint16_image_new(cols, rows, 3);
        omg->maxval = img->maxval;
        int col, row;
        for (row = 0; row < rows; row++)
          { for (col = 0; col < cols; col++)
              { /* Decide which channel {chn} this sample belongs to: */
                int chn = uint16_image_Pike_F100_sample_channel(col, row);
                /* Copy sample, set other channels to 0: */
                uint16_t smp = img->smp[row][col];
                uint16_t *po = &(omg->smp[row][3*col]);
                int k;
                for (k = 0; k < 3; k++) 
                  { (*po) = (uint16_t)(k == chn ? smp : 0); po++; }
              }
          }
      }
    else
      { /* Allocate a half-size 3-channel image: */
        omg = uint16_image_new(img->cols/2, img->rows/2, 3);
        omg->maxval = img->maxval;
        /* Fill it with the merged channels: */
        int col,row;
        for (row = 0; row < omg->rows; row++)
          { uint16_t *G0p = &(img->smp[2*row][0]);
            uint16_t *R0p = &(img->smp[2*row][1]);
            uint16_t *B0p = &(img->smp[2*row+1][0]);
            uint16_t *G1p = &(img->smp[2*row+1][1]);
            uint16_t *op = omg->smp[row];
            for (col = 0; col < omg->cols; col++)
              { (*op) = (*R0p); op++; R0p++; R0p++;
                int G0v = (int)(*G0p); G0p++; G0p++;
                int G1v = (int)(*G1p); G1p++; G1p++;
                /* Apply a simple dither matrix to avoid downwards bias in rounding: */
                int dither = (col + row) % 2;
                int Gv = (G0v + G1v + dither)/2;
                assert((Gv >= 0) && (Gv <= omg->maxval));
                (*op) = (uint16_t)Gv; op++;
                (*op) = (*B0p); op++; B0p++; B0p++;
              }
          }
      }
    return omg;
  }

uint16_image_t *uint16_image_Pike_F100_extract_bayer_channel(uint16_image_t *img, int bayer_col, int bayer_row, bool_t verbose)
  {
    int rows = img->rows;
    int cols = img->cols;
    demand(img->chns == 1, "input must be single-channel");
    uint16_image_t *omg = uint16_image_new(img->cols/2, img->rows/2, 1);
    omg->maxval = img->maxval;
    int col, row;
    for (row = bayer_row; row < rows; row += 2)
      { for (col = bayer_col; col < cols; col += 2)
          { uint16_t *pi = &(img->smp[row][col]);
            uint16_t *po = &(omg->smp[row >> 1][col >> 1]);
            (*po) = (*pi);
          }
      }
    return omg;
  }

void uint16_image_Pike_F100_color_balance(uint16_image_t *img, float bal[], bool_t verbose) 
  {
    int rows = img->rows;
    int cols = img->cols;
    int chns = img->chns;
    
    auto void print_bugs(char *title, int count[]);
    
    int nov[3] = {0,0,0}; /* Number of overflowed pixels per channel. */
    int nun[3] = {0,0,0}; /* Number of zeroed pixels per channel. */
    int col, row;
    for (row = 0; row < rows; row++)
      { for (col = 0; col < cols; col++)
          { int k;
            for (k = 0; k < chns; k++)
              { /* Choose the balance factor {fac}: */
                double fac;
                if (chns == 1)
                  { /* Image has not been debayered. */
                    /* Decide which channel {chn} this sample belongs to: */
                    int chn = uint16_image_Pike_F100_sample_channel(col, row);
                    fac = (double)bal[chn];
                  }
                else if (chns == 3)
                  { /* Image has been debayered. */
                    fac = (double)bal[k];
                  }
                /* Get hold of sample in channel {k}: */
                uint16_t *pi = &(img->smp[row][chns*col + k]);
                /* Apply factor: */
                double old = (double)(*pi);
                double new = floor(fac*old + 0.5);
                /* Count overflows/underflows: */
                if (new > img->maxval) { new = img->maxval; nov[k]++; }
                if ((new == 0) & (old > 0)) { nun[k]++; }
                (*pi) = (uint16_t)floor(new);
              }
          }
      }
    if (verbose || (nov > 0)) { print_bugs("oversaturated pixels", nov); }
    if (verbose || (nun > 0)) { print_bugs("zeroed pixels       ", nun); }
    return;
    
    void print_bugs(char *title, int count[])
      { 
        fprintf(stderr, "oversaturated pixels:");
        int k;
        for (k = 0; k < chns; k++)
          { fprintf(stderr, " %c = %8d", "RGB"[k], count[k]); }
        fprintf(stderr, "\n");
      }
  }
    
void uint16_image_Pike_F100_output_amp_balance(uint16_image_t *img, int split_col, double gain0, double gain1, bool_t verbose)
  {
    demand(gain0 >= 0, "invalid left-half gain");
    demand(gain1 >= 0, "invalid right-half gain");
    int col, row, chn;
    for (chn = 0; chn < img->chns; chn++)
      { int nov = 0; /* Number of overflowed samples. */
        for (row = 0; row < img->rows; row++)
          { uint16_t *p = img->smp[row] + chn;
            for (col = 0; col < img->cols; col++)
              { int val = (*p);
                double gain = (col < split_col ? gain0 : gain1);
                val = (int)floor(gain*((double)val) + 0.5);
                if (val > img->maxval)
                  { if (verbose) { fprintf(stderr, "  chn %d col %3d row %3d overflow val = %6d\n", chn, col, row, val); }
                    nov++;
                    val = img->maxval;
                  }
                (*p) = (uint16_t)val;
              }
          }
        if (nov > 0) { fprintf(stderr, "  chn = %d has %6d overflowed samples\n", chn, nov); }
      }
  }

void uint16_image_Pike_F100_interpolate_bayer(uint16_image_t *img, bool_t verbose)    
  { 
    int col, row;
    for (row = 0; row < img->rows; row++)
      { for (col = 0; col < img->cols; col++)
          { uint16_image_Pike_F100_interpolate_bayer_pixel(img, col, row); }
      }
  }
      
void uint16_image_Pike_F100_interpolate_bayer_pixel(uint16_image_t *img, int col, int row)   
  { 
    int rows = img->rows;
    int cols = img->cols;
    int chns = img->chns;
    demand(chns == 3, "cannot do interpolation, de-Bayer first");
    
    /* Get pointers to the previous, current, and next row: */
    uint16_t *prm = (row > 0 ? img->smp[row-1] : NULL);
    uint16_t *pro = img->smp[row];
    uint16_t *prp = (row < rows-1 ? img->smp[row+1] : NULL);
    
    /* Determine the Bayer pattern indices {ct,rt}: */
    int ct = (col % 2);
    int rt = (row % 2);
    
    /* Determine the channel {chn} that has actual data: */
    int chn = uint16_image_Pike_F100_sample_channel(col, row);
    /* Interpolate the remaining channels: */
    int k;
    for (k = 0; k < chns; k++)
      { 
        if (k != chn) 
          { /* Sample is missing. */
            /* Select the samples to average and their weights: */
            double v1 = 0, v2 = 0, v3 = 0, v4 = 0;
            double w1 = 0, w2 = 0, w3 = 0, w4 = 0;

            if ((k == 0) || (k == 2))
              { int par = (k == 0 ? 0 : 1);  /* Parity bit. */
                int cs = ct ^ par;
                int rs = rt ^ par;
                if ((cs == 0) && (rs == 0))
                  { /* Average of two row neighbors: */
                    if (col > 0)      { v1 = pro[3*(col-1)+k]; w1 = 1; } 
                    if (col < cols-1) { v2 = pro[3*(col+1)+k]; w2 = 1; } 
                  }
                else if ((cs == 0) && (rs == 1))
                  { /* Average of four diagonal neighbors: */
                    if ((col > 0)      && (row > 0))      { v1 = prm[3*(col-1)+k]; w1 = 1; } 
                    if ((col < cols-1) && (row > 0))      { v2 = prm[3*(col+1)+k]; w2 = 1; } 
                    if ((col > 0)      && (row < rows-1)) { v3 = prp[3*(col-1)+k]; w3 = 1; } 
                    if ((col < cols-1) && (row < rows-1)) { v4 = prp[3*(col+1)+k]; w4 = 1; } 
                  }
                else if ((cs == 1) && (rs == 1))
                  { /* Average of two column neighbors. */
                    if (row > 0)      { v1 = prm[3*col+k]; w1 = 1; } 
                    if (row < rows-1) { v2 = prp[3*col+k]; w2 = 1; }
                  }
                else
                  { assert(FALSE); }
              }
            else if (k == 1)
              { if (ct != rt)
                  { /* Average of four rook neighbors. */
                    if (row > 0)      { v1 = prm[3*col+k]; w1 = 1; } 
                    if (row < rows-1) { v2 = prp[3*col+k]; w2 = 1; } 
                    if (col > 0)      { v3 = pro[3*(col-1)+k]; w3 = 1; } 
                    if (col < cols-1) { v4 = pro[3*(col+1)+k]; w4 = 1; } 
                  }
              }

            /* Compute average: */
            double wtot = w1 + w2 + w3 + w4;
            assert (wtot > 0);
            double vmed = (w1*v1 + w2*v2 + w3*v3 + w4*v4)/wtot;
            /* Store it: */
            int iv = (int)floor(vmed + 0.5);
            assert(iv >= 0);
            assert(iv <= img->maxval);
            pro[3*col+k] = (uint16_t)iv;
          }
      }
  }

uint16_image_t *pnm_Pike_F100_bayer_channel_white_mask(uint16_image_t *img, bool_t verbose)
  {
    /* Convert the image to floating point: */
    bool_t isMask = FALSE; /* Assume smooth distr. of pixel values in encoding/decoding. */
    float_image_t *fim = float_image_from_uint16_image(img, isMask, NULL, NULL, TRUE, FALSE);
    /* Allocate an image for the high-frequency signal: */
    float_image_t *hif = float_image_new(img->chns, img->cols, img->rows);
    /* Allocate a monochrome working image for the low freq signal: */
    float_image_t *tmp = float_image_new(1, img->cols, img->rows);
    /* Create a very broad Gaussian kernel: */
    int hw = 64;     /* Window half-width */
    int nw = 2*hw-1; /* Window width */
    double wt[nw];   /* One-dimensional window weights. */
    wt_table_binomial_fill(nw, wt, NULL);
    wt_table_normalize_sum(nw, wt);
    double vRef = 4.0/256.0; /* A reasonable black level. */
    double base = 2.0;
    int chn;
    for (chn = 0; chn < img->chns; chn++)
      { 
        /* Convert channe {chn} of {fim} to log scale with a fixed vRef: */
        float_image_rescale_samples(fim, chn, 0.0f, 1.0f, (float)vRef, (float)(1.0+vRef));
        float_image_log_scale(fim, chn, 0.0, vRef, base);
        /* Extract the high-freq part, ignoring outliers: */
        int maxiter = 3;
        int iter = 0;       /* Iterations. */
        int nout;           /* Number of outliers in last iteration. */
        double avg, dev;    /* Average and deviation of channel {chn} of {hif}. */
        double maxd = 4;    /* Maximum deviations from mean allowed. */
        do 
          { /* Copy the relevant channel to {hif}: */
            float_image_set_channel(hif, chn, fim, chn);
            /* Blur {hif} by a very broad Gaussian kernel, put in {tmp}: */
            float_image_local_avg_var(hif, chn, hw, wt, tmp, 0, NULL, 0);
            /* Subtract the blurred image from {hif} giving the log of the mask: */
            float_image_mix_channels(-1.0, tmp, 0, 1.0, hif, chn, hif, chn);
            /* Find the average and variance of {hif}: */
            float_image_compute_sample_avg_dev(hif, chn, &avg, &dev);
            if (verbose) { fprintf(stderr, "  iteration %d avg = %+12.6f  dev = %12.6f\n", iter, avg, dev); }
            /* Check {hif} for outliers, fix them in {fim}: */
            int row, col;
            nout = 0;
            for (row = 0; row < img->rows; row++)
              { for (col = 0; col < img->cols; col++)
                  { float v = float_image_get_sample(hif, chn, col, row); 
                    if (fabs(v - avg) > maxd*dev)
                      { if (verbose)
                          { fprintf(stderr, "  %3d %3d", col, row);
                            fprintf(stderr, "  v = %12.6f = avg%+12.6f*dev", v, (v-avg)/dev);
                          }
                        /* Pull the sample in {fim} towards the local mean: */
                        float m = float_image_get_sample(tmp, 0, col, row); 
                        float f = (float)(m + avg + 0.75*maxd*dev);
                        float_image_set_sample(fim, chn, col, row, f);
                        nout++;
                      }
                  }
              }
            if (verbose || (nout > 0))
              { fprintf(stderr, "  found %d outliers\n", nout); }
            iter++;
          }
        while ((nout > 0) && (iter < maxiter));
        /* Convert the high-frequency signal from log to lin scale, so that {avg + maxd*dev --> 1}: */ 
        double uRef = 1.0/sample_conv_undo_log((float)(avg + maxd*dev), 0.0, 1.0, log(base));
        float_image_undo_log_scale(hif, chn, 0.0, uRef, base);
      }
    /* Convert back to PNM image: */
    uint16_image_t *pot = float_image_to_uint16_image(hif, isMask, img->chns, NULL, NULL, NULL, uint16_image_MAX_SAMPLE, TRUE, verbose);
    float_image_free(fim);
    float_image_free(hif);
    float_image_free(tmp);
    return pot;
  }
