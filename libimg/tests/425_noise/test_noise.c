#define PROG_NAME "test_noise"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-26 05:32:01 by stolfi */ 
/* Created on 2023-01-10 by J. Stolfi, UNICAMP */

#define test_noise_COPYRIGHT \
  "Copyright � 2024  by the State University of Campinas (UNICAMP)"
  
#define PROG_INFO \
  "For various image sizes {NX,NY}, filter" \
  " frequencies {fxFilter,fyFilter}, and filter complement flag, the program creates" \
  " images \"out/{ext}-{sizeTag}/img-{filterTag}.{ext}\" showing" \
  " the images generated by {float_image_noise} with those filter" \
  " arguments.\n" \
  "\n" \
  "  The sampls of the noise image will be compressed to the range {[0 _ 1]} when" \
  " writing out, either by a sigmoid squashing" \
  " function or by  simple affine scaling and clipping.\n" \
  "\n" \
  "  The {sizeTag} will be \"{NNNX}x{NNNY}\" where" \
  " {NNNX} and {NNNY} are the image dimensions {NX} and {NY}, formatted" \
  " as \"%04d\"; and  {ext} is \"pgm\" or \"ppm\".\n" \
  "\n" \
  "  The {filterTag} will be \"{XF}-{YF}-cm{CM}-sq{SQ}\" where {XF,YF} are the relative filter" \
  " frequencies {fxFilter/NX,fyFilter/NY} formatted as \"%08.6f\", and {CM} and {SQ} are the" \
  " complement flag and the sigmoid squashing flag, both formatted as \"F\" or \"T\".\n" \
  "\n" \
  "  For each output image, the program will also write a" \
  " file \"out/{sizeTag}/pwr-{filterTag}\" shoing the power" \
  " spectrum of the image."

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_hartley.h>
#include <float_image_hartley_spectrum.h>
#include <float_image_write_gen.h>

#include <float_image_noise.h>

int main(int argn, char **argv);

void tnoi_create_images_from_filter
  ( double rfxF,
    double rfyF,
    bool_t verbose
  );
  /* Generates all images using the filter parameters {rfxF*NX,rfyF*NX}.
    Note that {rfxF,rfyF} are relative frequencies, not absolute ones. */

void tnoi_create_images_from_size_and_filter
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double rfxF,
    double rfyF,
    bool_t verbose
  );
  /* Generates a noise image using {tnoi_create_noise_image} with {NC}
    channels (1 or 3) and size {NX} by {NY} and filter parameters
    {rfxF*NX,rfxF*NY}.
    
    Creates up to four versions of the image with
    {tnoi_create_noise_image}: with and without complemented filter, and
    with linear and non-linear squashing. The complemented-filter images
    are suppressed if the filter is trivial ({rfxF=rfyf=0}). The
    procedure then computes the log-scaled power spectrum of each
    version with {tnoi_compute_spectrum_images}.
    
    For each version, bot the noise image and its spectrum are written
    to disk as ".png" files with the names specified above. Most samples
    will have been mapped almost linearly to {[0 _ 1]} before writing.
    
    Note that {rfxF,rfyF} are relative frequencies, not absolute ones.
    Typically they should be 0.25 or less. */

float_image_t *tnoi_create_noise_image
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double fxFilter,
    double fyFilter,
    bool_t complement,
    bool_t verbose
  );
  /* Generates a noise image with {float_image_noise} and parameters
    {NC,NX,NY,fxFilter,fyFilter,complement}.
    
    Note that here {fxFilter,fyFilter} are signed absolute (not
    relative) frequencies, as expected by {float_image_noise}. */
    
float_image_t *tnoi_remap_image(float_image_t *img, bool_t squash, bool_t verbose);
  /* Rescales the samples of {img}, assumed to span some range symmetric about 0,
    to the {[-1 _ +1]} range. 
    
    If {squash} is true, uses a non-linear sigmoid map that goes from
    {-1} to {+1} as the original sample goes from {-INF} to {+INF}.
    Otherwise uses an affine map followed by clipping to {[-1 _ +1]}. In
    both cases, zero samples are mapped to to 0.5, and the slope of the
    mapping at 0 is chosen so that most values are mapped mostly
    linearly. */
    
float_image_t *tnoi_compute_spectrum_image(float_image_t *img, bool_t verbose);
  /*  Computes the power spectrum image of {img}.  The zero-frequency
    power sample will be cleared out, the samples are converted to log scale
    and mapped linearly to {[0 _ 1]}, with clipping of extreme values. */ 
  
void tnoi_compute_avg_rms(float_image_t *img, double *avg_P, double *rms_P);
  /* Computes the average value {avg} and the RMS value {rms}
    of the samples in all channels of {img}.  Returns the 
    results in {*avg_P}, {*rms_P}. */
 
void tnoi_write_image(char *sizeTag, char *prefix, char *filterTag, float_image_t *img, float vMin, float vMax);
  /* Writes the image {img} to "out/png-{sizeTag}/{prefix}-{filterTag}.png".
    The samples of {img} are suposed to range in {[vMin _ vMax]}. This range
    is mapped to {[0 _ 1]} in the file, with clipping. */

int main (int argn, char **argv)
  {
    fprintf(stderr, "choosing wave parameters...\n");
    bool_t verbose = TRUE;
    
    tnoi_create_images_from_filter(0.000, 0.000, verbose);
    tnoi_create_images_from_filter(0.010, 0.010, verbose);
    tnoi_create_images_from_filter(0.100, 0.100, verbose);
    tnoi_create_images_from_filter(0.200, 0.200, verbose);
   
    return 0;
  }

void tnoi_create_images_from_filter
  ( double rfxF,
    double rfyF,
    bool_t verbose
  )
  {
    tnoi_create_images_from_size_and_filter(1,  640,  480, rfxF, rfyF, verbose);
                                                                       
    tnoi_create_images_from_size_and_filter(3,  400,  400, rfxF, rfyF, verbose);
                                                                       
    tnoi_create_images_from_size_and_filter(1,  128,  128, rfxF, rfyF, verbose);
    tnoi_create_images_from_size_and_filter(1,  256,  256, rfxF, rfyF, verbose);
    tnoi_create_images_from_size_and_filter(1,  512,  512, rfxF, rfyF, verbose);
    tnoi_create_images_from_size_and_filter(1, 1024, 1024, rfxF, rfyF, verbose);
    
  }

void tnoi_create_images_from_size_and_filter
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double rfxF,
    double rfyF,
    bool_t verbose
  )
  {
    double fxFilter = rfxF*NX;
    double fyFilter = rfyF*NY;

    char *sizeTag = jsprintf("%04dx%04d", NX, NY);
  
    for (uint32_t kcm = 0;  kcm <= 1; kcm++)
      { bool_t complement = (kcm == 1);
      
        if ((fxFilter != 0) || (fyFilter != 0) || (! complement))
          { 
            float_image_t *img = tnoi_create_noise_image(NC, NX, NY, fxFilter, fyFilter, complement, verbose);

            for (uint32_t ksq = 0;  ksq <= 1; ksq++)
              { bool_t squash = (ksq == 1);

                char *filterTag = NULL;
                char *filterTag = jsprintf("%08.6f-%08.6f-cm%c-sq%c", rfxF, rfyF, "FT"[complement], "FT"[squash]);

                float_image_t *sqz = tnoi_remap_image(img, squash, verbose);
                tnoi_write_image(sizeTag, "img", filterTag, sqz, -1.0, +1.0);

                float_image_t *pwr = tnoi_compute_spectrum_image(sqz, verbose);
                tnoi_write_image(sizeTag, "pwr", filterTag, pwr, 0.0, 1.0);

                free(filterTag);
                float_image_free(pwr);
                float_image_free(sqz);
              }

            float_image_free(img);
          }
      }
    free(sizeTag);
  }  

float_image_t *tnoi_create_noise_image
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    double fxFilter,
    double fyFilter,
    bool_t complement,
    bool_t verbose
  )
  {
    if (verbose) 
      { fprintf(stderr, "creating noise image NC = %d size = %d x %d", NC, NX, NY);
        fprintf(stderr, " filter = %+.4f x  %+.4f complement = %c\n", fxFilter, fyFilter, "FT"[complement]);
      }
      
    float_image_t *img = float_image_noise(NC, NX, NY, fxFilter, fyFilter, complement);
    
    return img;
  }
  
float_image_t *tnoi_remap_image(float_image_t *img, bool_t squash,bool_t verbose)
  {
    /* Get image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    
    /* Compute the RMS value {dev} of the image samples: */
    double rms;
    tnoi_compute_avg_rms(img, NULL, &rms);
    
    float_image_t *sqz = float_image_copy(img);

    if (squash)
      { /* Nonlinear squashing to {[0_1]}: */
        double vSquash = 3.0*rms;
        fprintf(stderr, "squashing samples to [-1 _ +1], vSquash = %.6f\n", vSquash);
        for (uint32_t ic = 0;  ic < NC; ic++)
          { for (uint32_t ix = 0;  ix < NX; ix++)
              { for (uint32_t iy = 0;  iy < NY; iy++)
                  { double smp = float_image_get_sample(img, ic, ix, iy);
                    smp = smp/vSquash;
                    smp = smp/hypot(1, smp);
                    float_image_set_sample(sqz, ic, ix, iy, (float)smp);
                  }
              }
          }
      }
    else
      { /* Affine scaling to {[0_1]}: */
        float vMax = (float)(3.0*rms);
        float vMin = -vMax;
        fprintf(stderr, "rescaling samples from [%8.6f _ %8.6f] to [-1 _ +1]\n", vMin, vMax);
        for (uint32_t kc = 0;  kc < NC; kc++)
          { float_image_rescale_samples(sqz, kc, vMin, vMax, -1.0f, +1.0f); }
      }

    return sqz;
  }

float_image_t *tnoi_compute_spectrum_image(float_image_t *img, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "computing power spectrum image\n"); }
    
    /* Get image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);

    /* Conpute the Hartley transform image: */
    float_image_t *har = float_image_new(NC, NX, NY);
    float_image_hartley_transform(img, har);

    /* Conpute the Hartley power spectrum image: */
    bool_t centered = TRUE;
    float_image_t *pwr = float_image_new(NC, NX, NY);
    float_image_hartley_spectrum(har, pwr, centered);
    
    /* Clear out the zero frequency term: */
    for (uint32_t ic = 0;  ic < NC; ic++)
      { float_image_set_sample(pwr, ic, NX/2, NY/2, 0.0); }
  
    /* Compute the RMS value {dev} of the image samples: */
    double avg;
    tnoi_compute_avg_rms(pwr, &avg, NULL);
    if (verbose) { fprintf(stderr, "average power spectrum sample = %13.5e\n", avg); }
    
    /* Get maximum power in original spectrum: */
    float vMax = 1.0e-38f;
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { float vmc = float_image_spectrum_max_sample(pwr, ic, centered);
        vMax = fmaxf(vMax, vmc);
      }

    /* Convert spectrum to logscale: */
    double bias = 1.0e-8*avg;
    double vRef = 1.0e-8*avg;
    if (verbose) { fprintf(stderr, "converting to log scale, bias = %13.5e vRef = %13.5e\n", bias, vRef); }
    double base = M_E;
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { float_image_log_scale(pwr, ic, bias, vRef, base); }
    
    /* Map original values below {vRef} to zero: */
    float vMin = sample_conv_log((float)vRef, bias, vRef, log(base));
    vMax = sample_conv_log(vMax, bias, vRef, log(base));
    if (verbose) { fprintf(stderr, "scaling log spectrum from [%+10.6f _ %+10.6f] to [0 _ 1]\n", vMin, vMax); }
    for (uint32_t ic = 0;  ic < NC; ic++) 
      { float_image_rescale_samples(pwr, ic, vMin, vMax, 0.0, 1.0); }
      
    float_image_free(har);
    
    return pwr;
  }
   
void tnoi_compute_avg_rms(float_image_t *img, double *avg_P, double *rms_P)
  { 
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    
    double sum_s = 0; /* Sum of samples. */
    double sum_s2 = 0; /* Sum of samples squared. */
    int32_t NS = 0; /* Number of finite samples found. */
    for (uint32_t ic = 0;  ic < NC; ic++)
      { int32_t NS_ch;
        double sum_s_ch = float_image_compute_sample_sum(img, ic, &NS_ch);
        double sum_s2_ch = float_image_compute_squared_sample_sum(img, ic, 0.0, NULL);
        sum_s += sum_s_ch;
        sum_s2 += sum_s2_ch;
        NS += NS_ch;
      }
    assert(NS == NC*NX*NY);
    double avg = sum_s/NS;
    double rms = sqrt(sum_s2/NS);
    assert(isfinite(avg));
    assert(isfinite(rms));
    if (avg_P != NULL) { (*avg_P) = avg; }
    if (rms_P != NULL) { (*rms_P) = rms; }
  }

void tnoi_write_image(char *sizeTag, char *prefix, char *filterTag, float_image_t *img, float vMin, float vMax)
  {
    char *ext = "png";
    
    mkdir("out", 0755);
    char *sizeDir = NULL; char *sizeDir = jsprintf("out/%s-%s", ext, sizeTag);
    mkdir(sizeDir, 0755);
    char *imgName = NULL; char *imgName = jsprintf("%s/%s-%s.%s", sizeDir, prefix, filterTag, ext);
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaEnd = 1.000;
    double bias = 0.000;
    bool_t verbose = TRUE;
    float_image_write_gen_named(imgName, img, ffmt, vMin, vMax, gammaEnd, bias, verbose);

    free(sizeDir);
    free(imgName);
  }
 
