#define PROG_NAME "test_spectrum_table"
#define PROG_DESC "test of {spectrum_table_exact.h}, {spectrum_table_binned.h}, {spectrum_table_convert.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-30 08:04:30 by stolfi */ 
/* Created on 2008-10-05 by J. Stolfi, UNICAMP */

#define test_spectrum_table_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <gauss_table.h>
#include <spectrum_table_convert.h>
#include <spectrum_table_binned.h>
#include <spectrum_table_exact.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <float_image_paint.h>
#include <float_image_hartley.h>
#include <float_image_hartley_spectrum.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <float_image.h>

#define BT_ENC_EXPO sample_conv_gamma_BT709_ENC_EXPO
#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS 
  /* Values of {expo} and {bias} parameters for {sample_conv_gamma}
    that approximate the BT.709 encoding.  */

typedef enum { kind_PEAK, kind_WAVE, kind_BUMP, kind_REAL } kind_t;
#define kind_LAST (kind_REAL)

static char *kind_name[kind_LAST+1] =
  { [kind_PEAK] = "peak",
    [kind_WAVE] = "wave", 
    [kind_BUMP] = "bump",
    [kind_REAL] = "real"
  };

int32_t main(int32_t argn, char **argv);

void do_test
  ( int32_t chns,
    int32_t cols,
    int32_t rows,
    kind_t kind,
    int32_t kx, 
    int32_t ky,
    int32_t bins
  );
  /* Tests the spectrum table functions on a test image of type {kind}
    and parameters {kx,ky}. See {get_input_image} below.
    
    Uses {bins} bins in the binned spectrum table.
    
    Writes all output files with names
    "out/{kind_name[kind]}-{nx}x{ny}-{kx}-{ky}-{S}-{TAG}-{FIG}.{EXT}" where
    "{S}-{TAG}" is "0-orig", "1-pwrs", "2-tbex" or "3-tbbn", {FIG} is
    a figure label ("b", "cr", etc.), and {EXT} is "txt", "ppm", or "pgm". */
  
float_image_t *get_test_image(int32_t chns, int32_t cols, int32_t rows, kind_t kind, int32_t kx, int32_t ky);
  /* Generates an image of type {kind} and parameters {kx,ky} for tests.
    The parameters {kx,ky} will be automatically reduced to the image's
    domain.

    If {kind} is {kind_REAL} the image is read from
    "in/real-{kx}-{ky}-0-orig.ppm" or "in/real-{kx}-{ky}-0-orig.pgm";
    the given {cols,rows} is ignored. Otherwise the function creates a
    synthetic image with size {cols,rows}. 
    
    Also writes the image to "out/{kind_name[kind]}-{nx}x{ny}-{kx}-{ky}-0-orig.{EXT}"
    where {EXT} is "ppm" or "pgm" depending on {chns}. */

void compute_total_image_energy(float_image_t *img, bool_t squared, double *nsmp_P, double *terg_P);
  /* Computes the total number of samples {nsmp} in all channels of {img} and a sum of 
    {terg} of all those samples, if {squared} is false, or of their squares, if {squared} is true.
    Returns them in {*nsmp_P,*terg_P}. */

void check_power_totals(char *tname, int32_t zent, double zntrm, double zterg, double insmp, double iterg);
  /* Prints the number of entries {zent}, the number of terms {zntrm}
     and total energy {zterg} of table {tbname}, and checks whether
     they match those of the input image {insmp,iterg}. */
    
float_image_t *read_image(char *dir, int32_t kx, int32_t ky, int32_t chns, char *suffix);
  /* Reads an image file called "{dir}/real-{kx}-{ky}-{suffix}.{EXT}",
    where {kx,ky} are formatted as "%04d", and {EXT} is "ppm" or "pgm". */
  
void write_image(char *dir, kind_t kind, int32_t kx, int32_t ky, char *suffix, float_image_t *img, bool_t pwr, bool_t centered);
  /* Writes the image {img} as file "{dir}/{kind_name[kind]}-{kx}-{ky}-{suffix}.{ext}" 
     where {kx,ky} are formatted as "%04d" and {ext} is "ppm" or "pgm".
     
     If {pwr} is false, assumes that {img} is an ordinary image or Fourier transform. If all samples
     are non-negative, maps 0.0 to 0.0 and the max sample value to 1.0. If there are negative samples,
     maps {-vMax} to 0 and {+vMax} to 1.0, where {vMax} is the max absolute sample value.
     The {centered} parameter is ignored.
     
     If {pwr} is true, assumes that it is a power spectrum image, and uses logarithmic scale.  
     Adjust the scale so that the max entry, excluding the constant term, is mapped to a 
     bit over 1.0. Assumes that the constant term is at indices {0,0} if {centered}
     is false, and {cols/2,rows/2} if {centered} is true.
     
     In either case, the image is then adjusted by the standard PNM gamma encoding
     and discretized with 0.0 mapped to 0 and 1.0 mapped to {maxval}. */ 

void write_spectrum_exact
  ( char *dir, 
    kind_t kind, 
    int32_t cols, 
    int32_t rows, 
    int32_t kx, 
    int32_t ky, 
    char *suffix, 
    spectrum_table_exact_t *tx
  );
  /* Writes {tx} to "{dir}/{kind_name[kind]}-{cols}x{rows}-{kx}-{ky}-{suffix}.txt". */

void write_spectrum_binned
  ( char *dir, 
    kind_t kind, 
    int32_t cols, 
    int32_t rows, 
    int32_t kx, 
    int32_t ky, 
    char *suffix, 
    spectrum_table_binned_t *tb
  );
  /* Writes {tb} to "{dir}/{kind_name[kind]}--{cols}x{rows}{kx}-{ky}-{suffix}.txt". */

int32_t main (int32_t argn, char **argv)
  {
    /* Get the test image type {kind}: */
    kind_t kind = 0; 
    while ((kind < kind_LAST) && (strcmp(kind_name[kind],argv[1])!= 0)) { kind++; }
    assert(strcmp(kind_name[kind],argv[1])== 0);
    /* Get the image dimensions {chns,cols,rows}: */
    int32_t chns = atoi(argv[2]);
    int32_t cols = atoi(argv[3]);
    int32_t rows = atoi(argv[4]);
    /* Get the image indices {kx, ky}: */
    int32_t kx = atoi(argv[5]);
    int32_t ky = atoi(argv[6]);
    /* Get the bin count {bins}: */
    int32_t bins = atoi(argv[7]);

    /* Run the test: */
    do_test(chns, cols, rows, kind, kx, ky, bins);

    return 0;
  }

void do_test
  ( int32_t chns,
    int32_t cols,
    int32_t rows,
    kind_t kind,
    int32_t kx, 
    int32_t ky,
    int32_t bins
  )
  { bool_t verbose = ((cols <= 9) && (rows <= 9));
    
    fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "testing with image of type %s\n", kind_name[kind]); 
    
    float_image_t *iimg = get_test_image(chns, cols, rows, kind, kx, ky);

    double insmp;
    double iterg;
    compute_total_image_energy(iimg, TRUE, &insmp, &iterg);

    fprintf(stderr, "computing the Hartley transform ...\n");
    float_image_t *fimg = float_image_new(chns, cols, rows); 
    float_image_hartley_transform(iimg, fimg);
    /* Check of total energy is preserved: */
    double fnsmp; /* Total samples in Hartley transform. */
    double fterg; /* Sum of squares of samples in Hartley transform. */
    compute_total_image_energy(fimg, TRUE, &fnsmp, &fterg);
    int32_t fent = chns*cols*rows;
    check_power_totals("Hartley transform image", fent, fnsmp, fterg, insmp, iterg);
    
    fprintf(stderr, "computing the power spectrum ...\n");
    float_image_t *pimg = float_image_new(chns, cols, rows);
    bool_t center = (cols + rows + kx + ky) % 2 == 0;
    fprintf(stderr, "spectrum centering = %c\n", "FT"[center]);
    float_image_hartley_spectrum(fimg, pimg, center);
    write_image("out", kind, kx, ky, "1-pwrs", pimg, TRUE, center);
    double pnsmp; /* Total samples in pwoer spectrum image. */
    double pterg; /* Sum of entries in power spectrum image. */
    compute_total_image_energy(pimg, FALSE, &pnsmp, &pterg);
    int32_t pent = chns*cols*rows;
    check_power_totals("power spectrum image", pent, pnsmp, pterg, insmp, iterg);

    fprintf(stderr, "gathering the exact spectrum table...\n");
    spectrum_table_exact_t tx = spectrum_table_exact_new(0);
    for (int32_t c = 0;  c < chns; c++) 
      { spectrum_table_exact_append_all(pimg, c, center, &tx, verbose); }
    spectrum_table_exact_sort(&tx, verbose);
    write_spectrum_exact("out", kind, cols, rows, kx, ky, "2-tbex", &tx);

    fprintf(stderr, "checking exact spectrum table totals...\n");
    int32_t xent = (int32_t)tx.ne;
    double xnsmp = 0.0;
    double xterg = 0.0;
    for (int32_t k = 0;  k < tx.ne; k++) 
      { xnsmp += tx.e[k].nTerms; 
        xterg += tx.e[k].power;
      }
    check_power_totals("exact spectrum table", xent, xnsmp, xterg, insmp, iterg);
    
    fprintf(stderr, "gathering the binned spectrum table...\n");
    spectrum_table_binned_t tb = spectrum_table_binned_make((uint32_t)bins);
    for (int32_t c = 0;  c < chns; c++) 
      { spectrum_table_binned_add_all(pimg, c, center, &tb, verbose); }
    write_spectrum_binned("out", kind, cols, rows, kx, ky, "3-tbbn", &tb);

    fprintf(stderr, "checking binned spectrum table totals...\n");
    int32_t bent = (int32_t)tb.ne;
    double bnsmp = 0.0;
    double bterg = 0.0;
    for (int32_t k = 0;  k < tb.ne; k++) 
      { bnsmp += tb.e[k].nTerms; 
        bterg += tb.e[k].power; 
      }
    check_power_totals("binned spectrum table", bent, bnsmp, bterg, insmp, iterg);
    
    fprintf(stderr, "checking conversion of exact spectrum to binned spectrum...\n");
    spectrum_table_binned_t tc = spectrum_table_convert_exact_to_binned(&tx, (uint32_t)cols, (uint32_t)rows);
    write_spectrum_binned("out", kind, cols, rows, kx, ky, "4-tbcv", &tc);
 
    fprintf(stderr, "checking converted spectrum table totals...\n");
    int32_t cent = (int32_t)tc.ne;
    double cnsmp = 0.0;
    double cterg = 0.0;
    for (int32_t k = 0;  k < tc.ne; k++) 
      { cnsmp += tc.e[k].nTerms; 
        cterg += tc.e[k].power;
      }
    check_power_totals("converted spectrum table", cent, cnsmp, cterg, insmp, iterg);
    
    free(tc.e);
    free(tb.e);
    free(tx.e);
   
    float_image_free(pimg);
    float_image_free(fimg);
    float_image_free(iimg);

    fprintf(stderr, "------------------------------------------------------------\n");
  }  

float_image_t *get_test_image(int32_t chns, int32_t cols, int32_t rows, kind_t kind, int32_t kx, int32_t ky)
  {
    auto void reduce_indices(int32_t *kxp, int32_t *kyp);
      /* Reduces {*kxp} to {0..cols-1}, {*kyp} to {0..rows-1}. */
 
    auto void fit_bell(double *cp, double dp, int32_t n);
    /* Adjusts {*cp} so that a Gaussian with center {*cp} and deviation {dp} fits inside {[0_n]}. */
    
    fprintf(stderr, "creating the test image...\n");
    float_image_t *img;
    switch(kind)
      {
        case kind_PEAK:
          { reduce_indices(&kx, &ky);
            img = float_image_new(chns, cols, rows);
            float_image_fill(img, 0.0);
            double amp = 1.0; /* so that the total channel energy is 1. */
            float_image_fill_pixel(img, kx, ky, (float)amp);
          }
          break;
          
        case kind_WAVE:
          { reduce_indices(&kx, &ky);
            img = float_image_new(chns, cols, rows);
            double amp = sqrt(2.0); /* So that the channel energy is {cols*rows}. */
            float_image_hartley_wave(img, kx, ky, (float)amp);
          }
          break;
          
        case kind_BUMP:
          { reduce_indices(&kx, &ky);
            img = float_image_new(chns, cols, rows);
            float_image_fill(img, 0.0);
            double dx = 0.040*((double)cols);
            double dy = 0.025*((double)cols);
            /* Make sure that the whole smudge is inside the image: */
            double cx = kx + 0.5; fit_bell(&cx, dx, cols);
            double cy = ky + 0.5; fit_bell(&cy, dy, rows);
            for (int32_t c = 0;  c < chns; c++)
              { float val = (float)(0.1 + 0.4*c); /* Intensity in channel {c}. */
                (void)float_image_paint_smudge(img, c, cx, cy, dx, dy, val, 3);
              }
          }
          break;
          
        case kind_REAL:
          { demand((kx >= 0) && (ky >= 0), "bad image number");
            img = read_image("in", kx, ky, chns, "0-orig");
            assert(chns == img->sz[0]);
            /* Ignore the given {cols,rows} parameters: */
            cols = (int32_t)img->sz[1];
            rows = (int32_t)img->sz[2];
          }
          break;
      }
    fprintf(stderr, "image size chns = %d  cols = %d  rows = %d\n", chns, cols, rows); 
    write_image("out", kind, kx, ky, "0-orig", img, FALSE, FALSE);
       
    return img;
    
    /* INTERNAL IMPLEMENTATIONS */

    void reduce_indices(int32_t *kxp, int32_t *kyp)
      { (*kxp) = (*kxp) % cols; if ((*kxp) < 0) { (*kxp) += cols; }
        (*kyp) = (*kyp) % rows; if ((*kyp) < 0) { (*kyp) += rows; }
      }
    
    void fit_bell(double *cp, double dp, int32_t n)
      { double max_devs = gauss_table_BIG_ARG; /* Max devs that is worth considering. */
        double mrg = max_devs * dp + 2.0;
        if ((*cp) < mrg) { (*cp) = mrg; }
        if ((*cp) >= n - mrg) { (*cp) = n - mrg; }
      }
  }
    
void compute_total_image_energy(float_image_t *img, bool_t squared, double *nsmp_P, double *terg_P)
  {
    int32_t chns = (int32_t)img->sz[0];
    int32_t cols = (int32_t)img->sz[1];
    int32_t rows = (int32_t)img->sz[2];

    /* Compute the total energy of the original image: */
    int32_t npix = cols*rows;
    double nsmp = 0.0;
    double terg = 0.0;
    for (int32_t c = 0;  c < chns; c++) 
      { nsmp += cols*rows;
        if (squared) 
          { terg += float_image_compute_squared_sample_sum(img, c, 0.0, NULL); }
        else
          { terg += float_image_compute_sample_sum(img, c, NULL); }
      }
    fprintf(stderr, "input image:\n");
    fprintf(stderr, "  pixels = %9d\n", npix);
    fprintf(stderr, "  samples = %24.16e\n", nsmp);
    fprintf(stderr, "  total energy = %24.16e\n", terg);
    
    (*nsmp_P) = nsmp;
    (*terg_P) = terg;
  }

void check_power_totals(char *tname, int32_t zent, double zntrm, double zterg, double insmp, double iterg)
  {
    fprintf(stderr, "%s:\n", tname);
    fprintf(stderr, "  entries = %9d\n", zent);
    fprintf(stderr, "  power terms = %12.3f\n", zntrm);
    fprintf(stderr, "  total energy = %24.16e\n", zterg);
    double rnsmp = sqrt(zntrm/insmp);
    double rterg = sqrt(zterg/iterg);
    if ((fabs(rnsmp - 1.0) > 1.0e-6) || (fabs(rterg - 1.0) > 1.0e-6))
      { fprintf(stderr, "  relative errors:");
        fprintf(stderr, "  terms = %18.10f", rnsmp);
        fprintf(stderr, "  energy = %18.10f\n", rterg);
        demand(FALSE, "terms and/or energy mismatch");
      }
  }


float_image_t *read_image(char *dir, int32_t kx, int32_t ky, int32_t chns, char *suffix)
  { char *ext = (chns == 1 ? "pgm" : "ppm");
    char *fname = jsprintf("%s/real-%04d-%04d-%s.%s", dir, kx, ky, suffix, ext);
    FILE *rd = open_read(fname, TRUE);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    fclose(rd);
    free(fname);
    bool_t yUp = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yUp, verbose);
    uint16_image_free(pim);
    for (int32_t c = 0;  c < fim->sz[0]; c++) 
      { float_image_apply_gamma(fim, c, 1/BT_ENC_EXPO, BT_ENC_BIAS); }
    return fim;
  }

void write_image(char *dir, kind_t kind, int32_t kx, int32_t ky, char *suffix, float_image_t *img, bool_t pwr, bool_t centered)
  {
    int32_t chns = (int32_t)img->sz[0];
    int32_t cols = (int32_t)img->sz[1];
    int32_t rows = (int32_t)img->sz[2];
    
    /* Make a copy of the image, to avoid spoiling it: */
    float_image_t *fim = float_image_copy(img);
    
    float vMin, vMax;
    if (pwr)
      { vMin = (float)1.0e-38; vMax = vMin;
        for (int32_t c = 0;  c < chns; c++) 
          { if (! centered) 
              { /* Shift so that the mean value is at the center: */
                float_image_shift(fim, c, cols/2, rows/2);
              }
            vMax = fmaxf(vMax, float_image_spectrum_max_sample(fim, c, TRUE));
          }
        /* Nicefy and map to log scale: */
        double vRef = 1.0e-5*fmax(1.0e-38,vMax); /* Reference value for log scale. */
        double vBase = 10.0;                       /* Base for log scale. */
        fprintf(stderr, "vRef for log scale = %24.16e\n", vRef);
        for (int32_t c = 0;  c < chns; c++) 
          { float_image_log_scale(fim, c, 0.0, vRef, vBase); }
        vMin = 0.0;
        vMax = (float)(log(1.125*vMax/vRef)/log(vBase)); 
      }
    else
      { vMin = +INF; vMax = -INF;
        for (int32_t c = 0;  c < chns; c++)
          { float_image_update_sample_range(fim, c, &vMin, &vMax); }
        if ((vMin < 0) || (vMax < 0))
          { /* Range contains negative values, make it symmetric: */
            vMax = (float)fmax(1.0e-38, fmax(fabs(vMin), fabs(vMax)));
            vMin = -vMax;
          }
        else
          { /* Set the low end to 0.0: */
            vMin = 0.0;
          }
      }
    
    /* Map {vMin} to 0, {vMax} to 1, and apply view gamma: */
    fprintf(stderr, "nominal range before rescaling = [ %24.16e _ %24.16e]\n", (double)vMin, (double)vMax);
    for (int32_t c = 0;  c < chns; c++) 
      { float_image_rescale_samples(fim, c, vMin, vMax, 0.0, 1.0);
        float_image_apply_gamma(fim, c, BT_ENC_EXPO, BT_ENC_BIAS);
      }
    
    /* Quantize: */
    bool_t yUp = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, NULL, NULL, NULL, 255, yUp, verbose);
    
    /* Write to PPM file: */
    char *ext = (chns == 1 ? "pgm" : "ppm");
    char *fname = jsprintf("%s/%s-%04dx%04d-%04d-%04d-%s.%s", dir, kind_name[kind], cols, rows, kx, ky, suffix, ext);
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
    int32_t cols, 
    int32_t rows, 
    int32_t kx, 
    int32_t ky, 
    char *suffix, 
    spectrum_table_exact_t *tx
  )
  {
    char *fname = jsprintf("%s/%s-%04dx%04d-%04d-%04d-%s.txt", dir, kind_name[kind], cols, rows, kx, ky, suffix);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t k = 0;  k < tx->ne; k++)
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
    int32_t cols, 
    int32_t rows, 
    int32_t kx, 
    int32_t ky, 
    char *suffix, 
    spectrum_table_binned_t *tb
  )
  {
    char *fname = jsprintf("%s/%s-%04dx%04d-%04d-%04d-%s.txt", dir, kind_name[kind], cols, rows, kx, ky, suffix);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t k = 0;  k < tb->ne; k++)
      { spectrum_table_binned_entry_t *tbk = &(tb->e[k]);
        fprintf(wr, " %24.16e %24.16e", tbk->fmin, tbk->fmax); 
        fprintf(wr, "  %24.16e", tbk->fmid); 
        fprintf(wr, " %24.16e %24.16e", tbk->nTerms, tbk->power);
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
  }
