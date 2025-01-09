#define PROG_NAME "test_float_image_hist"
#define PROG_DESC "test of {float_image_hist.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-24 16:31:00 by stolfi */ 
/* Created on 2008-10-05 by J. Stolfi, UNICAMP */

#define test_float_image_hist_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <float_image.h>
#include <float_image_hist.h>
#include <uint16_image_read_pnm.h>
#include <sample_conv_gamma.h>
#include <float_image_from_uint16_image.h>

#define BT_ENC_EXPO sample_conv_gamma_BT709_ENC_EXPO
#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS 
  /* Values of {expo} and {bias} parameters for {sample_conv_gamma}
    that approximate the BT.709 encoding.  */

int32_t main(int32_t argn, char **argv);

void do_test(char *imgName, uint32_t N, bool_t doGamma, bool_t allBad);
  /* Tests the functions {float_image_hist_make} and 
    {float_image_hist_write_named} with the image file "in/imgName.ppm".
    Writes the histogram of each channel {c} to "in/imgName-ab{AB}-{NNNN}-{CC}.hist" where
    {AB} is {allBad} converted to 'T' or 'F', {CC} is the channel index {c} formated as "%02d",
    and {NNNN} is the histogram size formatted as "%04d". */
  
float_image_t *read_image(char *imgName, bool_t doGamma);
  /* Reads image 'in/{imgName}.ppm' and converts to a {float_image_t}.
    If {doGamma} is true assumes the BT.709 encoding,
    else assumes linear encoding (gamma = 1). */

int32_t main (int32_t argn, char **argv)
  {
    srandom(4615*417);
    
    do_test("frutasgd", 100, TRUE,  FALSE);
    do_test("iclogobk", 100, FALSE, FALSE);
    do_test("coloramp", 100, TRUE,  FALSE);
    do_test("langeva0", 100, FALSE, FALSE);
    do_test("noisegra", 100, FALSE, FALSE);
    do_test("noisegre", 100, FALSE, FALSE);
    do_test("gridsvga", 100, FALSE, FALSE);
    do_test("uniforan", 100, FALSE, FALSE);

    do_test("noisegra", 100, FALSE, TRUE);
    do_test("coloramp",   2, TRUE,  FALSE);

    return 0;
  }

void do_test(char *imgName, uint32_t N, bool_t doGamma, bool_t allBad)
  { fprintf(stderr, "------------------------------------------------------------\n");
    fprintf(stderr, "testing with imgName = '%s' doGamma = %d", imgName, "FT"[doGamma]); 
    fprintf(stderr, " N = %d allBad = %c\n", N, "FT"[allBad]); 
    
    float_image_t *img = read_image(imgName, doGamma);
    
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    
    double tol = 1.0e-12; /* Tolerance for {double} roundoff. */
    
    for (int32_t c = 0; c < NC; c++)
      { /* Set some values to {NAN} or {±INF}: */
        uint32_t nbad_exp = 0;
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { if (allBad || (drandom() < 0.05))
                  { float bad = (drandom() < 0.333 ? NAN : (drandom() < 0.5 ? -INF : +INF));
                    float_image_set_sample(img, c, x, y, bad);
                    nbad_exp++;
                  }
              }
          }
      
        /* Try sometimes with {cumh}, sometimes without: */
        double_vec_t hist = double_vec_new(0);
        double_vec_t cumh = double_vec_new(0);
        double_vec_t *hist_P = &(hist);
        double_vec_t *cumh_P = (drandom() < 0.5 ? &(cumh) : NULL);
        double hmin, hmax;
        uint32_t ngud, nbad;
        float_image_hist_build(img, c, N, &hmin, &hmax, hist_P, cumh_P, &ngud, &nbad);
        fprintf(stderr, "  channel %d: [%16.8e _ %16.8e] gud = %d bad = %d\n", c, hmin, hmax, ngud, nbad);
        demand((ngud + nbad) == NX*NY, "invalid {ngud,nbad}");
        demand(nbad == nbad_exp, "incorrect bad value count {nbad}");
        demand(hist.ne == N, "wrong {hist} size");
        if (cumh_P != NULL) { demand(cumh.ne == N, "wrong {cumh} size"); }
        demand(isfinite(hmin) && isfinite(hmax) && (hmin < hmax), "bad {hmin,hmax}");
        double sum = 0;
        for (uint32_t k = 0; k < N; k++)
          { double hk = hist.e[k];
            demand(isfinite(hk) && (hk >= 0), "bad {hist} value");
            sum += hk;
            if (cumh_P != NULL) { demand(fabs(cumh.e[k] - sum) <= tol, "bad {cumh} value"); }
          }
        demand(fabs(sum - (double)ngud) <= ngud*tol, "sum mismatch");
        
        /* Write the histograms: */
        char *fname_hist = jsprintf("out/%s-ab%c-%04d-%02d.hist", imgName, "FT"[allBad], N, c);
        fprintf(stderr, "  writing histograms to %s ...\n", fname_hist);
        float_image_hist_write_named(fname_hist, hmin, hmax, &hist);
        /* Cleanup: */
        free(hist.e);
        if (cumh.e != NULL) { free(cumh.e); }
      }
        
    float_image_free(img);

    fprintf(stderr, "------------------------------------------------------------\n");
  }  

float_image_t *read_image(char *imgName, bool_t doGamma)
  { char *fname = jsprintf("in/%s.ppm", imgName);
    FILE *rd = open_read(fname, TRUE);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    fclose(rd);
    free(fname);
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yup, verbose);
    if (doGamma)
      { for (int32_t c = 0;  c < fim->sz[0]; c++) 
          { float_image_apply_gamma(fim, c, 1/BT_ENC_EXPO, BT_ENC_BIAS); }
      }
    uint16_image_free(pim);
    return fim;
  }
