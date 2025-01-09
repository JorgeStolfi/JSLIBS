#define PROG_NAME "test_map_channels"
#define PROG_DESC "test of {float_image_map_channels.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-01 16:06:11 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_map_channels_COPYRIGHT \
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
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image.h>

#include <float_image_map_channels.h>

int32_t main(int32_t argn, char **argv);

void tm_do_test(float_image_t *imgA, int32_t NCB);
  /* Tests conversion of {imgA} to {NCB} channels. */

float_image_t *tm_read_pnm(char *fname);
  /* Reads an image from a PNM (PGM or PPM) file with name "{fname}". */

void tm_write_fni(char *outPrefix, float_image_t *img);
  /* Writes {img} to disk as "{outPrefix}.fni". */

void tm_write_pnm(char *outPrefix, float_image_t *img);
  /* Writes {img} to disk as "{outPrefix}.pgm" if it has 1 channel,
    or "{outPrefix}.ppm" if it has 3 channels. */

int32_t main (int32_t argn, char **argv)
  {
    /* Read input image: */
    float_image_t *imgA = tm_read_pnm("in/test-RGB.ppm");
    int32_t NCA;
    float_image_get_size(imgA, &NCA, NULL, NULL);
    assert(NCA == 3);
    
    tm_do_test(imgA, 1);
    tm_do_test(imgA, NCA);
    tm_do_test(imgA, NCA+2);
    
    return 0;
  }
    
void tm_do_test( float_image_t *imgA, int32_t NCB)
  {
    int32_t NCA, NX, NY;
    float_image_get_size(imgA, &NCA, &NX, &NY);

    /* Create output image: */
    float_image_t *imgB = float_image_new(NCB, NX, NY);

    /* Test RGB to YUV conversion: */
    float_image_map_channels_RGB_to_YUV(imgA, imgB);
    
    char *outPrefix = jsprintf("out/test-YUV-%02d", NCB);
    
    tm_write_fni(outPrefix, imgB);
    if ((NCB == 1) || (NCB == 3)) { tm_write_pnm(outPrefix, imgB); }
    
    free(outPrefix);
    float_image_free(imgB);
  }

float_image_t *tm_read_pnm(char *fname)
  {
    bool_t isMask = FALSE;
    double gamma = 1.000;
    double bias = 0.000;
    bool_t yup = FALSE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named("in/test-RGB.ppm", isMask, gamma, bias, yup, warn, verbose);
    return img;
  }
  
void tm_write_fni(char *outPrefix, float_image_t *img)
  { char *fname = jsprintf("%s.fni", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, img);
    fclose(wr);
    free(fname);
  }

void tm_write_pnm(char *outPrefix, float_image_t *img)
  { int32_t NCB;
    float_image_get_size(img, &NCB, NULL, NULL);
    demand((NCB == 1) || (NCB == 3), "invalid channel count");
    char *ext = (NCB == 1 ? "pgm" : "ppm");
    char *fname = jsprintf("%s.%s", outPrefix, ext);
    bool_t isMask = FALSE;
    double expoDec = 1.000;
    double bias = 0.000;
    bool_t yup = FALSE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_write_pnm_named(fname, img, isMask, expoDec, bias, yup, warn, verbose);
  }
