/* Test of jspnm.h, uint16_image.h */
/* Last edited on 2024-12-26 15:12:35 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>

int32_t main (int32_t argc, char **argv);

void do_uint16_image_io_tests(char *name);
  /* Performs various image I/O tests, reading files
    with names "{name}-rd-{VAR}-{MXV}.{EXT}" and writing files
    with names "{name}-wr-{VAR}-{MXV}.{EXT}", where 
    {VAR} is "raw" or "txt", {MXV} is "uni" or "sma" or "big", and
    {EXT} is "pbm" or "pgm" or "ppm". */

void do_uint16_image_io_test(char *name, bool_t raw, bool_t big, uint32_t kind);
  /* Reads an image from file "{name}-rd-{VAR}-{MXV}.{EXT}" and writes 
    its complementary image as "{name}-wr-{VAR}-{MXV}.{EXT}".

    The {kind} may be 0, 1, or 2, meaning PBM, PGM, or PPM file format,
    respectively.
    
    If {kind} is 0, then {EXT} will be "pbm"; {big} must be false,
    {MXV} will be "uni", and the files should have {maxval==1}.
    
    If {kind} is 1 or 2, then {EXT} will be "pgm" or "ppm",
    respectively; {MXV} will be "big" if {big=TRUE} (expects {256 <=
    maxval <= 65535}), or "sma" if {big=FALSE} (expects {1 <= maxval
    <= 255}).
    
    In any case, {VAR} will be "raw" if {raw=TRUE} (assumes that the
    files have `raw' (binary) samples), or "txt" if {raw=FALSE}
    (assumes `plain' (decimal) samples). */
    
void ressublime_image(uint16_image_t *img, uint16_image_t *omg);  
  /* Copies {img} ino {omg}, complementing some pixels
    relative to {img.maxval}, flipping left-right and
    top-bottom. Also sets {omg.maxval = img.maxval}. */  

int32_t main (int32_t argc, char **argv)
  { do_uint16_image_io_tests("out/test");
    return 0;
  }
  
void do_uint16_image_io_tests(char *name)
  { for (uint32_t raw = 0; raw <= 1; raw++)
      for (uint32_t big = 0; big <= 1; big++)
        for (uint32_t kind = 0; kind <= 2; kind++)
          { do_uint16_image_io_test(name, (raw!=0), (big!=0), kind); }
  }
      
void do_uint16_image_io_test(char *name, bool_t raw, bool_t big, uint32_t kind)
  { 
    /* Filename extension: */
    char *extbl[3] = {"pbm", "pgm", "ppm"};
    char *EXT = extbl[kind]; 
    /* Tag indicating format variant: */
    char *VAR = (raw ? "raw" : "txt"); 
    /* Tag indicating maxval range: */
    char *MXV = (kind == 0 ? "uni" : (big ? "big" : "sma")); 
    
    fprintf(stderr, "-----------------------------------------\n");
    fprintf(stderr, "testing name = %s parms = %s %s %s\n", name, VAR, MXV, EXT);
    fprintf(stderr, "\n");

    char *fname = jsprintf("%s-rd-%s-%s.%s", name, VAR, MXV, EXT);
    uint16_image_t *img = uint16_image_read_pnm_named(fname, TRUE);
    uint16_image_describe(stderr, fname, img);
    free(fname);

    uint16_image_t *omg = uint16_image_new(img->cols, img->rows, img->chns);
    ressublime_image(img, omg);
    fprintf(stderr, "\n");

    fname = jsprintf("%s-wr-%s-%s.%s", name, VAR, MXV, EXT);
    uint16_image_describe(stderr, fname, omg);
    uint16_image_write_pnm_named(fname, omg, (!raw), TRUE);
    free(fname);
    fprintf(stderr, "-----------------------------------------\n");
  }

void ressublime_image(uint16_image_t *img, uint16_image_t *omg)
  {
    uint32_t cols = img->cols; assert(img->cols == omg->cols);
    uint32_t rows = img->rows; assert(img->rows == omg->rows);
    uint32_t chns = img->chns; assert(img->chns == omg->chns);

    omg->maxval = img->maxval;
    
    double rad = 0.45*(cols < rows ? cols : rows);

    for (int32_t y = 0; y < rows; y++)
      { double dy = y + 0.5 - 0.5*rows;
        uint16_t *ip = img->smp[y];
        uint16_t *op = omg->smp[(int32_t)rows - 1 - y];
        int32_t ik = 0;
        int32_t ok = (int32_t)(cols*chns) - 1;
        for (int32_t x = 0; x < cols; x++)
          { double dx = x + 0.5 - 0.5*cols;
            bool_t rev = (dx*dx + dy*dy < rad*rad);
            for (int32_t c = 0; c < chns; c++)
              { op[ok] = (uint16_t)(rev ? img->maxval - ip[ik] : ip[ik]);
                ik++; ok--;
              }
          }
      }
  }
