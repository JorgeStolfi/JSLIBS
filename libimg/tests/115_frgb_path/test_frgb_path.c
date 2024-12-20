#define PROG_NAME "test_frgb_path"
#define PROG_DESC "test of various functions from {frgb_path.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-05 07:59:42 by stolfi */
/* Created on 2023-03-03 by J. Stolfi, UNICAMP */

#define test_frgb_ops_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  convert(1), gimp(1), display(1), ppm(1), pgm(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2023-03-03 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_frgb_ops_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  This program reads nothing and writes some messages."

#define PROG_INFO_OPTS \
  "  None."

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include <frgb_ops.h>
#include <frgb_path.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jspnm.h>
#include <bool.h>
#include <affirm.h>
#include <argparser.h>

typedef struct options_t
  { bool_t op; /* placeholder. */
  } options_t;
  
typedef void tfpt_tr_t(frgb_t *p);
typedef double tfpt_fn_t(frgb_t *p);

int32_t main(int32_t argc, char **argv);

options_t *tfpt_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

void tfpt_do_tests(options_t *o);
  /* Creates a test image using various
    {float_image_gradient.h} tools. */
  
typedef frgb_t path_proc_t(double s, int32_t cycles, int32_t style);
  /* Type of a path procedure expected by {test_frgb_path}. */

void test_frgb_path(char *tag, path_proc_t path, int32_t cycles, int32_t style);

int32_t main(int32_t argc, char **argv)
  {
    /* Parse the command line options: */
    options_t *o = tfpt_parse_options(argc, argv);

    tfpt_do_tests(o);
    free(o); o = NULL;

    return 0;
  }

void tfpt_do_tests(options_t *o)
  {
    auto frgb_t unsigned_path(double s, int32_t cycles, int32_t style);
    auto frgb_t signed_path(double s, int32_t cycles, int32_t style);
    
    for (int32_t cycles = -5; cycles <= +5; cycles++)
      { int32_t max_u_style = frgb_path_map_unsigned_max_style();
        for (int32_t style = 0; style <= max_u_style; style ++)
          { test_frgb_path("u", unsigned_path, cycles, style); }
        int32_t max_s_style = frgb_path_map_signed_max_style();
        for (int32_t style = 0; style <= max_s_style; style ++)
          { if ((style != 0) || (cycles == +1))
              { test_frgb_path("s", signed_path, cycles, style); }
          }
      }
    return;
    
    frgb_t unsigned_path(double s, int32_t cycles, int32_t style)
      { assert((s >= 0) && (s <= 1));
        return frgb_path_map_unsigned(s, cycles, style);
      }
      
    frgb_t signed_path(double s, int32_t cycles, int32_t style)
      { assert((s >= -1) && (s <= +1));
        return frgb_path_map_signed(s, cycles, style);
      }
      
  }

void test_frgb_path(char *tag, path_proc_t path, int32_t cycles, int32_t style)
  {
    int32_t N = 1000; /* Number of sample points along path. */
    
    char *prefix = jsprintf("out/path_%s%d_c%+03d", tag, style, cycles);
    
    char *txname = jsprintf("%s.txt", prefix);
    FILE *txwr = open_write(txname, TRUE);
    
    char *imname = jsprintf("%s.ppm", prefix);
    FILE *imwr = open_write(imname, TRUE);
    int32_t NX = 20;
    int32_t NY = N;
    uint16_t maxval = PNM_FILE_MAX_MAXVAL;
    fprintf(imwr, "P3\n");
    fprintf(imwr, "%d %d\n", NX, NY);
    fprintf(imwr, "%u\n", maxval);
    
    for (int32_t iy = 0;  iy <= N; iy++)
      { double s = ((double)iy)/((double)N);
        if (tag[0] == 's') { s = 2*s - 1; }
        frgb_t f = path(s, cycles, style);
        frgb_t g = f; frgb_to_YUV(&g); frgb_YUV_to_YHS(&g);
        fprintf(txwr, "%5d %+12.8f", iy, s);
        fprintf(txwr, " %7.4f %7.4f %7.4f", f.c[0], f.c[1], f.c[2]);
        fprintf(txwr, " %7.4f %7.4f %7.4f\n", g.c[0], g.c[1], g.c[2]);
        uint16_t dv[3];
        for (int32_t ic = 0;  ic < 3; ic++) 
          { if ((f.c[ic] < 0.0) || (f.c[ic] > 1.0)) 
              { fprintf(stderr, "!! f[%d] = %+14.8f for s = %+14.8f\n", ic, f.c[ic], s); }
            int32_t dval = 1 + (int32_t)floor(f.c[ic] * (maxval-2) + 0.5);
            if (dval < 0) { dval = 0; }
            if (dval > maxval) { dval = maxval; }
            dv[ic] = (uint16_t)dval;
          }
        
        for (int32_t ix = 0;  ix < NX; ix++) 
          { for (int32_t ic = 0;  ic < 3; ic++) 
              { fprintf(imwr, " %u", dv[ic]); }
          }
        fprintf(imwr, "\n");
      }
    fclose(imwr);
    fclose(txwr);
    
    free(imname);
    free(txname);
    free(prefix);
  }

options_t *tfpt_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
