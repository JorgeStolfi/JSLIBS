/* Test of jspng.h, uint16_image_io_png.h */
/* Last edited on 2023-02-12 09:58:11 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <fget.h>
#include <jsmath.h>
#include <jspnm.h>
#include <uint16_image.h>

#include <uint16_image_read_png.h>
#include <uint16_image_write_png.h>

#define MAX_CHNS (uint16_image_read_png_MAX_CHNS)
      
#define MIN_GAMMA (1.0/1024.0)
#define MAX_GAMMA (1024.0)

int32_t main (int32_t argc, char **argv);

void do_uint16_image_io_png_own_tests(char *iDir, char *oDir);
  /* Performs various PNG I/O tests.
  
    Namely reads files with names "{iDir}/test-{CSP}-A{ALP}-{MXV}.png"
    and writing files with names "{oDir}/test-{CSP}-A{ALP}-{MXV}.png",
    where {CSP} is the colorspace ("GRY", "RGB", or "MAP"), {ALP} is the
    alpha channel flag ("0" or "1"), and {MXV} is the 5-digit pixel
    maxval ("00001", "00003", "00015", "00255", or "65535"). */

void do_uint16_image_io_png_own_test(char *iDir, char *oDir, int32_t csp, int32_t alp, int32_t mxv, FILE *wr);
  /* Reads an image from file "{iDir}/test-{CSP}-{ALP}-{MXV}.png" and writes the
    frobnicated image as "{oDir}/test-{CSP}-A{alp}-{MXV}.png".

    The {csp} parameter indicates the colorspace. Currently the only
    valid choices are 0 (grayscale, in which case {CSP} is "GRY"), 1
    (RGB color, in which case {CSP} is "RGB"), or 2 (indexed RGB color,
    in which case CSP is "MAP").
    
    The {alp} flag indicates the presence (1) or absence (0) of an alpha
    channel.
    
    The {mxv} parameter is the pixel's max value (in any case {MXV} is
    {mxv} with five digits). */
    
void do_uint16_image_io_png_official_tests(char *iDir, char *oDir);
  /* Reads a list "test-images.dir" of image names and attribues 
    from directory {iDir} and tests {uint16_image_io_png_read}
    on each. */

void do_uint16_image_io_png_official_test
  ( char *iDir,
    char *iName, 
    int32_t NX,  /* Num of columns. */                                                       
    int32_t NY,  /* Num of rows. */                                                          
    int32_t NC,  /* Number of channels of image after mapping (1 to 4). */                   
    int32_t BY,  /* True bit depth of luminance channel (for GRAY or GRAY+ALPHA images). */  
    int32_t BR,  /* True bit depth of red channel (for RGB or RGB+ALPHA images). */          
    int32_t BG,  /* True bit depth of green channel (for RGB or RGB+ALPHA images). */        
    int32_t BB,  /* True bit depth of blue channel (for RGB or RGB+ALPHA images). */         
    int32_t BA,  /* True bit depth of alpha channel (for images with ALPHA). */              
    double iGamma, /* Gamma in file. */
    FILE *wr  /* File for summary. */
  );
  /* Reads image {iName} from directory {iDir}
    and checks whether it is consistent with {NX,NY,
    NC,BY,BR,BG,BB,BA,iGamma}. */

bool_t is_valid_comb(int32_t csp, int32_t alp, int32_t mxv);
  /* Returns true iff {csp,alp,mxv} is a valid combination
    for a PNG image. */

void frobnicate_image(uint16_image_t *img, uint16_image_t *omg);  
  /* Copies {img} into {omg}, complementing some pixels
    relative to {img.maxval}, flipping left-right and
    top-bottom. Also sets {omg.maxval = img.maxval}. */

void compare_images(uint16_image_t *img, uint16_image_t *omg);  
  /* Compares {img} and {omg}. Bombs out if different. */

int32_t main (int32_t argc, char **argv)
  { do_uint16_image_io_png_official_tests("00-DATA/official", "out/official");
    do_uint16_image_io_png_own_tests("00-DATA/stolfi", "out/stolfi");
    return 0;
  }
  
#define N_MAXVALS 5
static int32_t maxval[N_MAXVALS] = { 1, 3, 15, 255, 65535};

void do_uint16_image_io_png_own_tests(char *iDir, char *oDir)
  { 
    char *oName = NULL;
    asprintf(&oName, "%s/test-results.dir", oDir);
    FILE *wr = open_write(oName, TRUE);
    
    fprintf(wr, "# generated by {test_jspng.c}\n");
    fprintf(wr, "#\n");
    fprintf(wr, "# Filename               NX NY  NC BY BR BG BB BA Gamma\n");
    fprintf(wr, "# ---------------------  -- --  -- -- -- -- -- -- -------\n");

    int32_t csp; /* Basic color space: 0 = GRAY, 1 = explicit RGB, 2 = mampped RGB. */
    int32_t alp; /* Alpha channel: 0 = no, 1 = yes. */
    int32_t imx; /* The primary bit size is {2^imx}. */
    for (csp = 0; csp < 3; csp++)
      { for (alp = 0; alp < 2; alp++)
          { for (imx = 0; imx < N_MAXVALS; imx++)
              { int32_t mxv = maxval[imx]; /* Maxval of image. */
                if (is_valid_comb(csp,alp,mxv))
                  { do_uint16_image_io_png_own_test(iDir, oDir, csp, alp, mxv, wr); }
              }
          }
      }
    fclose(wr);
    free(oName);
  }

bool_t is_valid_comb(int32_t csp, int32_t alp, int32_t mxv)
  { 
    if (csp == 0)
      { /* Explicit grayscale: */
        if (alp == 0)
          { return (mxv == 1) || (mxv == 3) || (mxv == 15) || (mxv == 255) || (mxv == 65535); }
        else
          { return (mxv == 255) || (mxv == 65535); }
      }
    else if (csp == 1)
      { /* Explicit RGB: */
        return (mxv == 255) || (mxv == 65535);
      }
    else if (csp == 2)
      { /* Mapped RGB: */
        if (alp == 1) { /* We are unable to create test images, so: */ return FALSE; }
        return (mxv == 1) || (mxv == 3) || (mxv == 15) || (mxv == 255);
      }
    else 
      { assert(FALSE); }
  }

void do_uint16_image_io_png_own_test(char *iDir, char *oDir, int32_t csp, int32_t alp, int32_t mxv, FILE *wr)
  { 
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "testing color = %d  alpha = %d  maxval = %d\n", csp, alp, mxv);

    /* File names: */
    char *csp_tbl[3] = {"GRY", "RGB", "MAP"};
    char *iName = NULL;
    asprintf(&iName, "%s/test-%s-A%d-%05d.png", iDir, csp_tbl[csp], alp, mxv);
    char *oName = NULL;
    asprintf(&oName, "%s/test-%s-A%d-%05d.png", oDir, csp_tbl[csp], alp, mxv);

    /* Read input image: */
    fprintf(stderr, "reading test image...\n");
    double iGamma;
    uint32_t fMaxval[MAX_CHNS];
    uint16_image_t *img = uint16_image_read_png_named(iName, &iGamma, fMaxval, TRUE);
    fprintf(stderr, "------------------------------------------------------------\n");
    uint16_image_describe(stderr, iName, img);
    fprintf(stderr, "------------------------------------------------------------\n");
    
    /* Check attributes: */
    if ((csp == 0) || (csp ==1))
      { /* Explicit GRAY or RGB, possibly with ALPHA:  */
        demand(img->chns == 1 + 2*csp + alp, "unexpected {chns}");
        demand(img->maxval == mxv, "unexpected {maxval}");
      }
    else if (csp == 2)
      { /* Mapped RGB, possibly with ALPHA: */
        demand(img->chns == 3 + alp, "unexpected {chns}");
        demand(img->maxval == 255, "unexpected {maxval}");
      }
      
    /* Check Gamma: */
    demand(isnan(iGamma) || ((iGamma >= MIN_GAMMA) && (iGamma <= MAX_GAMMA)), "bad {gamma}"); 

    /* Modify it: */
    fprintf(stderr, "modifying image...\n");
    uint16_image_t *omg = uint16_image_new(img->cols, img->rows, img->chns);
    frobnicate_image(img, omg);
    fprintf(stderr, "------------------------------------------------------------\n");
    uint16_image_describe(stderr, "modified", omg);
    fprintf(stderr, "------------------------------------------------------------\n");

    /* Write it out: */
    fprintf(stderr, "writing it out...\n");
    double oGamma = iGamma;
    uint16_image_write_png_named(oName, omg, oGamma, TRUE);
    
    /* Read it back: */
    fprintf(stderr, "reading it back in...\n");
    double bGamma;
    uint32_t bmaxval[MAX_CHNS];
    uint16_image_t *bmg = uint16_image_read_png_named(oName, &bGamma, bmaxval, TRUE);
    fprintf(stderr, "------------------------------------------------------------\n");
    uint16_image_describe(stderr, oName, bmg);
    fprintf(stderr, "------------------------------------------------------------\n");
    
    /* Check attributes: */
    fprintf(stderr, "checking...\n");
    compare_images(omg, bmg);
    if (isnan(oGamma))
      { demand(isnan(bGamma), "gamma should be {NAN}"); }
    else
      { demand(fabs(bGamma-oGamma)/sqrt(bGamma*oGamma) < 0.000001, "inconsistent {gamma}"); }
   
    uint32_t iMaxval[MAX_CHNS]; /* Maxvals expected from {BY,BR,BG,BB,BA}. */
    int32_t NC = img->chns;
    int32_t icY = ( (NC == 1) || (NC == 2) ? 0 : -1 ); /* Index of intensity channel, or -1 if none. */
    int32_t icR = ( NC < 3 ? -1 : 0 ); /* Index of red channel, or -1 if none. */
    int32_t icG = ( NC < 3 ? -1 : 1 ); /* Index of red channel, or -1 if none. */
    int32_t icB = ( NC < 3 ? -1 : 2 ); /* Index of red channel, or -1 if none. */
    int32_t icA = ( (NC == 2) || (NC == 4) ? NC-1 : -1 ); /* Index of red channel, or -1 if none. */
    
    if (icY >= 0) { iMaxval[icY] = mxv; }
    if (icR >= 0) { iMaxval[icR] = mxv; }
    if (icG >= 0) { iMaxval[icG] = mxv; }
    if (icB >= 0) { iMaxval[icB] = mxv; }
    if (icA >= 0) { iMaxval[icA] = mxv; }
    
    for (int32_t ic = 0; ic < NC; ic++)
      { if (fMaxval[ic] != iMaxval[ic]) 
          { fprintf(stderr, "** %s: channel %d - maxval %u should be %u\n", iName, ic, fMaxval[ic], iMaxval[ic]); }
      }

    /* Write summary line: */
    fprintf(wr, "| %s %5d %5d", iName, img->rows, img->cols);
    fprintf(wr, "  ? ??  ??");
    fprintf(wr, "  %d", img->chns);
    fprintf(wr, " %02d", (icY < 0 ? 0 : minbits(fMaxval[icY])));
    fprintf(wr, " %02d", (icR < 0 ? 0 : minbits(fMaxval[icR])));
    fprintf(wr, " %02d", (icG < 0 ? 0 : minbits(fMaxval[icG])));
    fprintf(wr, " %02d", (icB < 0 ? 0 : minbits(fMaxval[icB])));
    fprintf(wr, " %02d", (icA < 0 ? 0 : minbits(fMaxval[icA])));
    fprintf(wr, "  %7.5f", iGamma);
    fprintf(wr, "\n");

    uint16_image_free(img);
    uint16_image_free(omg);
    uint16_image_free(bmg);
    
    free(iName);
    free(oName);
    fprintf(stderr, "test successful!\n");
    fprintf(stderr, "============================================================\n");
  }

void do_uint16_image_io_png_official_tests(char *iDir, char *oDir)
  { 
    bool_t debug = FALSE;
    char *dName = NULL;
    asprintf(&dName, "%s/test-images.dir", iDir);
    FILE *rd = open_read(dName, TRUE);
    
    char *oName = NULL;
    asprintf(&oName, "%s/test-results.dir", oDir);
    FILE *wr = open_write(oName, TRUE);
    
    fprintf(wr, "# generated by {test_jspng.c}\n");
    fprintf(wr, "#\n");
    fprintf(wr, "# Filename      NX NY  NC BY BR BG BB BA Gamma\n");
    fprintf(wr, "# ------------  -- --  -- -- -- -- -- -- -------\n");
    
    int32_t nline = 0;
    while (TRUE)
      { nline++;
        if (debug) { fprintf(stderr, "parsing line %d...\n", nline); }
        bool_t ok = fget_test_comment_or_eol(rd, '#');
        if (ok) { continue; }
        if (fget_test_eof(rd)) { break; }
        /* Found something, not space, break, or EOF: */ 
        char ch = fget_char(rd);
        if (debug) { fprintf(stderr, "found characters '%c' = \\%03o\n", ch, ch); }
        if (ch == '|') 
          { /* Read data about one test image: */
            char *iName = fget_string(rd); /*Image file name, with extension. */
            if (debug) { fprintf(stderr, "image = %s\n", iName); }
        
            int32_t NX = fget_int32(rd);   /* Num of columns. */
            int32_t NY = fget_int32(rd);   /* Num of rows. */
            
            /* For palette-mapped images before expansion: */
            int32_t NI = fget_int32(rd);   /* 1 if colormapped, 0 if true color */
            int32_t BI = fget_int32(rd);   /* If colormapped, bits per palette index; 0 otherwise. */
            
            /* For all images, including palette-mapped after expansion: */
            int32_t BS = fget_int32(rd);   /* Nominal bits per sample of image after mapping. */
            
            int32_t NC = fget_int32(rd);   /* Number of channels of image after mapping (1 to 4). */
            
            int32_t BY = fget_int32(rd);   /* True bit depth of luminance channel (for GRAY or GRAY+ALPHA images). */
            int32_t BR = fget_int32(rd);   /* True bit depth of red channel (for RGB or RGB+ALPHA images). */
            int32_t BG = fget_int32(rd);   /* True bit depth of green channel (for RGB or RGB+ALPHA images). */
            int32_t BB = fget_int32(rd);   /* True bit depth of blue channel (for RGB or RGB+ALPHA images). */
            int32_t BA = fget_int32(rd);   /* True bit depth of alpha channel (for images with ALPHA). */                    
            
            double gamma = fget_double(rd); /* Gamma . */
            
            if (debug) { fprintf(stderr, "skipping blanks comments eol...\n"); }
            
            fget_comment_or_eol(rd, '#');
            
            /* Write expected data summary to {stderr}: */
            fprintf(stderr, "| %s %5d %5d", iName, NX, NY);
            fprintf(stderr, "  %d %02d  %02d", NI, BI, BS);
            fprintf(stderr, "  %d %02d %02d %02d %02d %02d", NC, BY, BR, BG, BB, BA);
            fprintf(stderr, "  %7.5f", gamma);
            fprintf(stderr, "\n");
            
            do_uint16_image_io_png_official_test(iDir, iName, NX, NY, NC, BY, BR, BG, BB, BA, gamma, wr);
            
            free(iName);
          }
        else
          { fprintf(stderr, "%s:%d ** invalic character '%c' (%02x)\n", dName, nline, ch, (uint32_t)ch);
            ungetc(ch, rd);
            fget_skip_to_eol(rd);
          }
      }
      
    fclose(rd);
    fclose(wr);
    free(dName);
    free(oName);
  }
  
void do_uint16_image_io_png_official_test
  ( char *iDir,
    char *iName, 
    int32_t NX,  
    int32_t NY,  
    int32_t NC, 
    int32_t BY,  
    int32_t BR,  
    int32_t BG,  
    int32_t BB,  
    int32_t BA,  
    double iGamma,
    FILE *wr
  ) 
  {
    fprintf(stderr, "=== %s  =======================================================\n", iName);
    
    /* Read the image: */
    char *fName = NULL;
    asprintf(&fName, "%s/%s", iDir, iName); 
    double fGamma;
    uint32_t fMaxval[MAX_CHNS]; /* Maxvals read from image, accounting for 'sBits' */
    bool_t verbose = TRUE;
    uint16_image_t *img = uint16_image_read_png_named(fName, &fGamma, fMaxval, verbose);
    
    /* Check data: */
    if ((img->cols != NX) || (img->rows != NY)) 
      { fprintf(stderr, "** %s: size %d x %d should be %d x %d\n", iName, img->cols, img->rows, NX, NY); }
      
    if (img->chns != NC) 
      { fprintf(stderr, "** %s: num channels %d should be %d\n", iName, img->chns, NC); }
      
    bool_t iGamma_def = (! isnan(iGamma)) && (iGamma > 0);
    bool_t fGamma_def = (! isnan(fGamma)) && (fGamma > 0);
    if ((iGamma_def != fGamma_def) || (iGamma_def && fGamma_def && (fabs(fGamma-iGamma) >= 0.00001)))
      { fprintf(stderr, "** %s: gamma %.6f should be %.6f\n", iName, fGamma, iGamma); }
      
    uint32_t iMaxval[MAX_CHNS]; /* Maxvals expected from {BY,BR,BG,BB,BA}. */
      
    int32_t icY = ( (NC == 1) || (NC == 2) ? 0 : -1 ); /* Index of intensity channel, or -1 if none. */
    int32_t icR = ( NC < 3 ? -1 : 0 ); /* Index of red channel, or -1 if none. */
    int32_t icG = ( NC < 3 ? -1 : 1 ); /* Index of red channel, or -1 if none. */
    int32_t icB = ( NC < 3 ? -1 : 2 ); /* Index of red channel, or -1 if none. */
    int32_t icA = ( (NC == 2) || (NC == 4) ? NC-1 : -1 ); /* Index of red channel, or -1 if none. */
    
    if (icY >= 0) { iMaxval[icY] = (1u << BY) - 1; }
    if (icR >= 0) { iMaxval[icR] = (1u << BR) - 1; }
    if (icG >= 0) { iMaxval[icG] = (1u << BG) - 1; }
    if (icB >= 0) { iMaxval[icB] = (1u << BB) - 1; }
    if (icA >= 0) { iMaxval[icA] = (1u << BA) - 1; }
    
    for (int32_t ic = 0; ic < NC; ic++)
      { if (fMaxval[ic] != iMaxval[ic]) 
          { fprintf(stderr, "** %s: channel %d - maxval %u should be %u\n", iName, ic, fMaxval[ic], iMaxval[ic]); }
      }

    /* Write summary line: */
    fprintf(wr, "| %s %5d %5d", iName, img->rows, img->cols);
    fprintf(wr, "  ? ??  ??");
    fprintf(wr, "  %d", img->chns);
    fprintf(wr, " %02d", (icY < 0 ? 0 : minbits(fMaxval[icY])));
    fprintf(wr, " %02d", (icR < 0 ? 0 : minbits(fMaxval[icR])));
    fprintf(wr, " %02d", (icG < 0 ? 0 : minbits(fMaxval[icG])));
    fprintf(wr, " %02d", (icB < 0 ? 0 : minbits(fMaxval[icB])));
    fprintf(wr, " %02d", (icA < 0 ? 0 : minbits(fMaxval[icA])));
    fprintf(wr, "  %7.5f", fGamma);
    fprintf(wr, "\n");

    uint16_image_free(img);
    free(fName);

    fprintf(stderr, "=====================================================================\n\n");
  }

void frobnicate_image(uint16_image_t *img, uint16_image_t *omg)
  {
    int32_t cols = img->cols; assert(img->cols == omg->cols);
    int32_t rows = img->rows; assert(img->rows == omg->rows);
    int32_t chns = img->chns; assert(img->chns == omg->chns);

    omg->maxval = img->maxval;
    
    int32_t skip = 20;
    
    int32_t invchns = ((chns == 2) || (chns == 4) ? chns - 1 : chns);

    int32_t x, y, c;
    for (y = 0; y < rows; y++)
      { uint16_t *ip = img->smp[y];
        uint16_t *op = omg->smp[rows - 1 - y];
        bool_t yinv = ((y >= skip) && (y < rows-skip));
        for (x = 0; x < cols; x++)
          { uint16_t *iv = ip + chns*x;
            uint16_t *ov = op + chns*(cols - 1 - x);
            bool_t xinv = ((x >= skip) && (x < cols-skip));
            for (c = 0; c < chns; c++)
              { bool_t cinv = (c < invchns);
                bool_t inv = (xinv && yinv && cinv);
                ov[c] = (uint16_t)(inv ? img->maxval - iv[c] : iv[c]);
              }
          }
      }
  }

void compare_images(uint16_image_t *img, uint16_image_t *omg)
  {
    int32_t cols = img->cols; demand(img->cols == omg->cols, "wrong cols");
    int32_t rows = img->rows; demand(img->rows == omg->rows, "wrong rows");
    int32_t chns = img->chns; demand(img->chns == omg->chns, "wrong chns");
    demand(omg->maxval == img->maxval, "wrong maxval");
        
    for (int32_t y = 0; y < rows; y++)
      { uint16_t *ip = img->smp[y];
        uint16_t *op = omg->smp[y];
        int32_t k = 0;
        for (int32_t x = 0; x < cols; x++)
          { for (int32_t c = 0; c < chns; c++)
              { if (ip[k] != op[k])
                  { fprintf(stderr, "img[%d][%d] = %d  omg[%d][%d] = %d\n", y, x, ip[k], y, x, op[k]); 
                    demand(FALSE, "samples differ!");
                  }
                k++;
              }
          }
      }
  }
