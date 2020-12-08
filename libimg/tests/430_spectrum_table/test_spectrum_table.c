#define PROG_NAME "test_spectrum_table"
#define PROG_DESC "test of {spectrum_table_exact.h}, {spectrum_table_binned.h}, {spectrum_table_convert.h}"
#define PROG_VERS "1.0"

/* Last edited on 2020-10-11 01:52:04 by jstolfi */ 
/* Created on 2008-10-05 by J. Stolfi, UNICAMP */

#define test_spectrum_table_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <bool.h>
#include <gauss_table.h>
#include <spectrum_table_convert.h>
#include <spectrum_table_binned.h>
#include <spectrum_table_exact.h>
#include <sample_conv.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <float_image_paint.h>
#include <float_image_hartley.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

#define BT_GAMMA (0.450) 
#define BT_BIAS (sample_conv_BT709_BIAS) 
/* Values of {gamma} and {bias} parameter that approximate
  the BT.709 encoding.  */

typedef enum { kind_PEAK, kind_WAVE, kind_BUMP, kind_REAL } kind_t;
#define kind_LAST (kind_REAL)

static char *kind_name[kind_LAST+1] =
  { [kind_PEAK] = "peak",
    [kind_WAVE] = "wave", 
    [kind_BUMP] = "bump",
    [kind_REAL] = "real"
  };

int main(int argn, char **argv);

void do_test
  ( int chns,
    int cols,
    int rows,
    kind_t kind,
    int kx, 
    int ky,
    int bins
  );
  /* Tests the spectrum table functions on a test image of type {kind}
    and parameters {kx,ky}. The parameters {kx,ky} will be
    automatically reduced to the image's domain. Uses {bins} bins in
    the binned spectrum table.

    If {kind} is {kind_REAL} the image is read from
    "in/real-{kx}-{ky}-0-orig.ppm" or "in/real-{kx}-{ky}-0-orig.pgm";
    the given {cols,rows} is ignored. Otherwise the function creates a
    synthetic image with size {cols,rows}.
    
    Writes all output files with names
    "out/{kind_name[kind]}-{nx}x{ny}-{kx}-{ky}-{S}-{TAG}-{FIG}.{EXT}" where
    "{S}-{TAG}" is "0-orig", "1-pwrs", "2-tbex" or "3-tbbn", {FIG} is
    a figure label ("b", "cr", etc.), and {EXT} is "txt", "ppm", or "pgm". */
  
float_image_t *read_image(char *dir, int kx, int ky, int chns, char *suffix);
  /* Reads an image file called "{dir}/real-{kx}-{ky}-{suffix}.{EXT}",
    where {kx,ky} are formatted as "%04d", and {EXT} is "ppm" or "pgm". */
  
void write_image(char *dir, kind_t kind, int kx, int ky, char *suffix, float_image_t *img, bool_t pwr);
  /* Writes the image {img} as file "{dir}/{kind_name[kind]}-{kx}-{ky}-{suffix}.{ext}" 
     where {kx,ky} are formatted as "%04d" and {ext} is "ppm" or "pgm".
     If {pwr} is TRUE, uses logarithmic scale.  Adjust the max intensity to 1.0
     and applies teh standard gamma encoding. */ 

void write_spectrum_exact
  ( char *dir, 
    kind_t kind, 
    int cols, 
    int rows, 
    int kx, 
    int ky, 
    char *suffix, 
    spectrum_table_exact_t *tx
  );
  /* Writes {tx} to "{dir}/{kind_name[kind]}-{cols}x{rows}-{kx}-{ky}-{suffix}.txt". */

void write_spectrum_binned
  ( char *dir, 
    kind_t kind, 
    int cols, 
    int rows, 
    int kx, 
    int ky, 
    char *suffix, 
    spectrum_table_binned_t *tb
  );
  /* Writes {tb} to "{dir}/{kind_name[kind]}--{cols}x{rows}{kx}-{ky}-{suffix}.txt". */

int main (int argn, char **argv)
  {
    /* Get the test image type {kind}: */
    kind_t kind = 0; 
    while ((kind < kind_LAST) && (strcmp(kind_name[kind],argv[1])!= 0)) { kind++; }
    assert(strcmp(kind_name[kind],argv[1])== 0);
    /* Get the image dimensions {chns,cols,rows}: */
    int chns = atoi(argv[2]);
    int cols = atoi(argv[3]);
    int rows = atoi(argv[4]);
    /* Get the image indices {kx, ky}: */
    int kx = atoi(argv[5]);
    int ky = atoi(argv[6]);
    /* Get the bin count {bins}: */
    int bins = atoi(argv[7]);

    /* Run the test: */
    do_test(chns, cols, rows, kind, kx, ky, bins);

    return 0;
  }

void do_test
  ( int chns,
    int cols,
    int rows,
    kind_t kind,
    int kx, 
    int ky,
    int bins
  )
  {
    int c, k;
    
    bool_t verbose = (cols <= 5) & (rows <= 5);
    
    auto void reduce_indices(int *kxp, int *kyp);
      /* Reduces {*kxp} to {0..cols-1}, {*kyp} to {0..rows-1}. */
       
    void reduce_indices(int *kxp, int *kyp)
      { (*kxp) = (*kxp) % cols; if ((*kxp) < 0) { (*kxp) += cols; }
        (*kyp) = (*kyp) % rows; if ((*kyp) < 0) { (*kyp) += rows; }
      }
    
    auto void fit_bell(double *cp, double dp, int32_t n);
    /* Adjusts {*cp} so that a Gaussian with center {*cp} and deviation {dp} fits inside {[0_n]}. */
    
    void fit_bell(double *cp, double dp, int32_t n)
      { double max_devs = gauss_table_BIG_ARG; /* Max devs that is worth considering. */
        double mrg = max_devs * dp + 2.0;
        if ((*cp) < mrg) { (*cp) = mrg; }
        if ((*cp) >= n - mrg) { (*cp) = n - mrg; }
      }

    fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "testing with image of type %s\n", kind_name[kind]); 
    
    fprintf(stderr, "creating the test image...\n");
    float_image_t *iimg;
    switch(kind)
      {
        case kind_PEAK:
          { reduce_indices(&kx, &ky);
            iimg = float_image_new(chns, cols, rows);
            float_image_fill(iimg, 0.0);
            double amp = 1.0; /* so that the total channel energy is 1. */
            float_image_fill_pixel(iimg, kx, ky, (float)amp);
          }
          break;
          
        case kind_WAVE:
          { reduce_indices(&kx, &ky);
            iimg = float_image_new(chns, cols, rows);
            double amp = sqrt(2.0); /* So that the channel energy is {cols*rows}. */
            float_image_hartley_wave(iimg, kx, ky, (float)amp);
          }
          break;
          
        case kind_BUMP:
          { reduce_indices(&kx, &ky);
            iimg = float_image_new(chns, cols, rows);
            float_image_fill(iimg, 0.0);
            double dx = 0.040*((double)cols);
            double dy = 0.025*((double)cols);
            /* Make sure that the whole smudge is inside the image: */
            double cx = kx + 0.5; fit_bell(&cx, dx, cols);
            double cy = ky + 0.5; fit_bell(&cy, dy, rows);
            for (c = 0; c < chns; c++)
              { float val = (float)(0.1 + 0.4*c); /* Intensity in channel {c}. */
                (void)float_image_paint_smudge(iimg, c, cx, cy, dx, dy, val, 3);
              }
          }
          break;
          
        case kind_REAL:
          { demand((kx >= 0) && (ky >= 0), "bad image number");
            iimg = read_image("in", kx, ky, chns, "0-orig");
            assert(chns == iimg->sz[0]);
            /* Ignore the given {cols,rows} parameters: */
            cols = (int)iimg->sz[1];
            rows = (int)iimg->sz[2];
          }
          break;
      }
    fprintf(stderr, "image size chns = %d  cols = %d  rows = %d\n", chns, cols, rows); 
    write_image("out", kind, kx, ky, "0-orig", iimg, FALSE);

    /* Compute the total energy of the original image: */
    int ient = cols*rows;
    double itrm = 0.0;
    double ierg = 0.0;
    for (c = 0; c < chns; c++) 
      { itrm += cols*rows;
        ierg += float_image_compute_total_energy(iimg,c,0.0);
      }
    fprintf(stderr, "input image:\n");
    fprintf(stderr, "  pixels = %9d\n", ient);
    fprintf(stderr, "  samples = %24.16e\n", itrm);
    fprintf(stderr, "  total energy = %24.16e\n", ierg);
    
    /* Allocate the transform/spectrum image: */
    float_image_t *timg = float_image_new(chns, cols, rows);
    
    auto void check_table_totals(char *tname, int zent, double ztrm, double zerg);
    /* Prints the numbr of entries {zent}, the number of terms {ztrm}
       and total energy {zerg} of table {tbname}, and checks whether
       they match those of the input image {itrm,ierg}. */

    void check_table_totals(char *tname, int zent, double ztrm, double zerg)
      {
        fprintf(stderr, "%s:\n", tname);
        fprintf(stderr, "  entries = %9d\n", zent);
        fprintf(stderr, "  Hartley terms = %12.3f\n", ztrm);
        fprintf(stderr, "  total energy = %24.16e\n", zerg);
        double rtrm = sqrt(ztrm/itrm);
        double rerg = sqrt(zerg/ierg);
        if ((fabs(rtrm - 1.0) > 1.0e-6) || (fabs(rerg - 1.0) > 1.0e-6))
          { fprintf(stderr, "  relative errors:");
            fprintf(stderr, "  terms = %18.10f", rtrm);
            fprintf(stderr, "  energy = %18.10f\n", rerg);
            demand(FALSE, "terms and/or energy mismatch");
          }
      }

    fprintf(stderr, "computing its energy spectrum...\n");
    float_image_hartley_transform(iimg, timg);
    int tent = cols*rows;
    double ttrm = 0.0;
    double terg = 0.0;
    for (c = 0; c < chns; c++) 
      { ttrm += cols*rows;
        terg += float_image_compute_total_energy(timg,c,0.0);
        float_image_square_samples(timg, c);
      }
    write_image("out", kind, kx, ky, "1-pwrs", timg, TRUE);
    check_table_totals("power spectrum", tent, ttrm, terg);

    fprintf(stderr, "gathering the exact spectrum table...\n");
    spectrum_table_exact_t tx = spectrum_table_exact_new(0);
    for (c = 0; c < chns; c++) 
      { spectrum_table_exact_append_all(timg, c, &tx, verbose); }
    spectrum_table_exact_sort(&tx, verbose);
    write_spectrum_exact("out", kind, cols, rows, kx, ky, "2-tbex", &tx);

    fprintf(stderr, "checking exact spectrum table totals...\n");
    int xent = tx.ne;
    double xtrm = 0.0;
    double xerg = 0.0;
    for (k = 0; k < tx.ne; k++) 
      { xtrm += tx.e[k].nTerms; xerg += tx.e[k].power; }
    check_table_totals("exact spectrum table", xent, xtrm, xerg);
    
    fprintf(stderr, "gathering the binned spectrum table...\n");
    spectrum_table_binned_t tb = spectrum_table_binned_make(bins);
    for (c = 0; c < chns; c++) 
      { spectrum_table_binned_add_all(timg, c, &tb, verbose); }
    write_spectrum_binned("out", kind, cols, rows, kx, ky, "3-tbbn", &tb);

    fprintf(stderr, "checking binned spectrum table totals...\n");
    int bent = tb.ne;
    double btrm = 0.0;
    double berg = 0.0;
    for (k = 0; k < tb.ne; k++) 
      { btrm += tb.e[k].nTerms; berg += tb.e[k].power; }
    check_table_totals("binned spectrum table", bent, btrm, berg);
    
    fprintf(stderr, "checking conversion of exact spectrum to binned spectrum...\n");
    spectrum_table_binned_t tc = spectrum_table_convert_exact_to_binned(&tx, cols, rows);
    write_spectrum_binned("out", kind, cols, rows, kx, ky, "4-tbcv", &tc);
 
    fprintf(stderr, "checking converted spectrum table totals...\n");
    int cent = tc.ne;
    double ctrm = 0.0;
    double cerg = 0.0;
    for (k = 0; k < tc.ne; k++) 
      { ctrm += tc.e[k].nTerms; cerg += tc.e[k].power; }
    check_table_totals("converted spectrum table", cent, ctrm, cerg);
    
    free(tc.e);
    free(tb.e);
    free(tx.e);
   
    float_image_free(timg);
    float_image_free(iimg);

    fprintf(stderr, "------------------------------------------------------------\n");
  }  

float_image_t *read_image(char *dir, int kx, int ky, int chns, char *suffix)
  { char *fname = NULL;
    char *ext = (chns == 1 ? "pgm" : "ppm");
    asprintf(&fname, "%s/real-%04d-%04d-%s.%s", dir, kx, ky, suffix, ext);
    FILE *rd = open_read(fname, TRUE);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    fclose(rd);
    free(fname);
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yup, verbose);
    uint16_image_free(pim);
    int c;
    for (c = 0; c < fim->sz[0]; c++) 
      { float_image_apply_gamma(fim, c, 1/BT_GAMMA, BT_BIAS); }
    return fim;
  }

void write_image(char *dir, kind_t kind, int kx, int ky, char *suffix, float_image_t *img, bool_t pwr)
  {
    int chns = (int)img->sz[0];
    int cols = (int)img->sz[1];
    int rows = (int)img->sz[2];
    
    /* Make a copy of the image, to avoid spoiling it: */
    float_image_t *fim = float_image_copy(img);
    
    /* Find true sample range {vMin,vMax}: */
    int c;
    float vMin = 0.0, vMax = -INF;
    for (c = 0; c < chns; c++)  { float_image_update_sample_range(fim, c, &vMin, &vMax); }
    
    if (pwr)
      { /* Nicefy and map to log scale: */
        double vRef = 1.0e-20*fmax(1.0e-100,vMax); /* Reference value for log scale. */
        double vBase = 10.0;                       /* Base for log scale. */
        fprintf(stderr, "vRef for log scale = %24.16e\n", vRef);
        for (c = 0; c < chns; c++) 
          { /* Remove mean value: */
            float_image_set_sample(fim, c, 0, 0, 0.0);
            /* Shift so that the mean value is at the center: */
            float_image_shift(fim, c, cols/2, rows/2);
            /* convert to log scale: */
            float_image_log_scale(fim, c, vRef, vBase);
          }
        /* Recompute {vMin,vMax}: */
        vMin = 0.0; vMax = -INF;
        for (c = 0; c < chns; c++) { float_image_update_sample_range(fim, c, &vMin, &vMax); }
      }
    
    /* Map 0 to 0, {vMax} to 1, and apply view gamma: */
    vMax = (float)fmax(1.0e-38, vMax); 
    fprintf(stderr, "vMax before rescaling = %24.16e\n", (double)vMax);
    for (c = 0; c < chns; c++) 
      { float_image_rescale_samples(fim, c, 0.0, vMax, 0.0, 1.0);
        float_image_apply_gamma(fim, c, BT_GAMMA, BT_BIAS);
      }
    
    /* Quantize: */
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, NULL, NULL, NULL, 255, yup, verbose);
    
    /* Write to PPM file: */
    char *fname = NULL;
    char *ext = (chns == 1 ? "pgm" : "ppm");
    asprintf(&fname, "%s/%s-%04dx%04d-%04d-%04d-%s.%s", dir, kind_name[kind], cols, rows, kx, ky, suffix, ext);
    FILE *wr = open_write(fname, TRUE);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    fclose(wr);
    
    /* Cleanup: */
    float_image_free(fim);
    uint16_image_free(pim);
    free(fname);
  }

void write_spectrum_exact
  ( char *dir, 
    kind_t kind, 
    int cols, 
    int rows, 
    int kx, 
    int ky, 
    char *suffix, 
    spectrum_table_exact_t *tx
  )
  {
    char *fname = NULL;
    asprintf(&fname, "%s/%s-%04dx%04d-%04d-%04d-%s.txt", dir, kind_name[kind], cols, rows, kx, ky, suffix);
    FILE *wr = open_write(fname, TRUE);
    int k;
    for (k = 0; k < tx->ne; k++)
      { spectrum_table_exact_entry_t *txk = &(tx->e[k]);
        fprintf(wr, " %20lu %20lu", txk->freq2.num, txk->freq2.den); 
        double freq_app = sqrt(((double)txk->freq2.num)/((double)txk->freq2.den));
        fprintf(wr, " %24.16e", freq_app); 
        fprintf(wr, " %24.16e %24.16e", txk->nTerms, txk->power);
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }

void write_spectrum_binned
  ( char *dir, 
    kind_t kind, 
    int cols, 
    int rows, 
    int kx, 
    int ky, 
    char *suffix, 
    spectrum_table_binned_t *tb
  )
  {
    char *fname = NULL;
    asprintf(&fname, "%s/%s-%04dx%04d-%04d-%04d-%s.txt", dir, kind_name[kind], cols, rows, kx, ky, suffix);
    FILE *wr = open_write(fname, TRUE);
    int k;
    for (k = 0; k < tb->ne; k++)
      { spectrum_table_binned_entry_t *tbk = &(tb->e[k]);
        fprintf(wr, " %24.16e %24.16e", tbk->fmin, tbk->fmax); 
        fprintf(wr, "  %24.16e", tbk->fmid); 
        fprintf(wr, " %24.16e %24.16e", tbk->nTerms, tbk->power);
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }
