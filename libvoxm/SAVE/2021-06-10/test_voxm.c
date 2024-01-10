#define PROG_NAME "test_voxm"
#define PROG_DESC "Test the voxel-based 3D modeling routines of {voxm.h}"
#define PROG_VERS "1.0"

#define test_voxm_C_COPYRIGHT \
  "Copyright © 2016 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-06-09 16:29:59 by jstolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -size {SIZE} \\\n" \
  "    -object { objs | tubes } \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes to standard output a tomogram (3D voxel array) containing an antialiased model of some test objects.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {SIZE}\n" \
  "    This mandatory argument defines the number of voxels" \
  " along the X, Y, and Z axes, respectively.\n" \
  "\n" \
  "  -object {OBJ}\n" \
  "    This mandatory argument defines the test object to render.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  salamic(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2016-03-16 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  " 2016-03-16 J. Stolfi: Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_voxm_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <r3.h>
#include <i3.h>
#include <r3x3.h>
#include <r3_motion.h>
#include <ppv_array.h>
#include <ppv_io.h>

#include <voxm_obj.h>
#include <voxm_splat.h>
#include <voxm_splat_tube.h>

#include <test_voxm_objs.h>
#include <test_voxm_tubes.h>
#include <test_voxm_mark.h>

/* COMMAND-LINE OPTIONS */

typedef struct test_voxm_options_t
  { int32_t size;  /* Size of voxel array. */
    char *object;  /* Name of test object. */
  } test_voxm_options_t;

#define test_voxm_size_MIN (100)
  /* Min tomogram size along any axis. */

#define test_voxm_size_MAX (3000)
  /* Max tomogram size along any axis. */

/* INTERNAL PROTOTYPES */

    /* The unit for all linear dimensions is the voxel side. */ 
    
test_voxm_options_t *test_voxm_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {test_voxm_options_t}. */

void test_voxm_tomogram_write(char *outFile, ppv_array_t *a); 

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    test_voxm_options_t *o = test_voxm_parse_options(argc, argv);
    
    /* Array dimensions: */
    int32_t NX = 3*o->size; 
    int32_t NY = o->size + o->size/2; 
    int32_t NZ = o->size; 
    
    /* Allocate the voxel array: */
    ppv_size_t sz[ppv_array_NAXES];
    int32_t k;
    for (k = 0; k < ppv_array_NAXES; k++) { sz[k] = 1; }
    sz[0] = NZ;
    sz[1] = NY;
    sz[2] = NX;
    ppv_nbits_t bps = 8;
    ppv_nbits_t bpw = 32;
    ppv_array_t a = ppv_new_array(sz, bps, bpw);
    
    /* Center and radius of array: */
    r3_t ctr = (r3_t){{ 0.5*NX, 0.5*NY, 0.5*NZ }};
    r3_t rad = (r3_t){{ 0.5*NX, 0.5*NY, 0.5*NZ }};
    double fuzzR = 1.5;      /* Half-thickness of fuzzy layer. */

    /* Mark the corners: */
    test_voxm_mark_corners(&a, &ctr, &rad, fuzzR);
    test_voxm_mark_edges(&a, &ctr, &rad, 0, fuzzR);
    test_voxm_mark_edges(&a, &ctr, &rad, 1, fuzzR);
    test_voxm_mark_edges(&a, &ctr, &rad, 2, fuzzR);

    if (strcmp(o->object, "objs") == 0)
      { test_voxm_objs(&a, &ctr, &rad, fuzzR); }
    else if (strcmp(o->object, "tubes") == 0)
      { test_voxm_tubes(&a, &ctr, &rad, fuzzR); }
    else 
      { demand(FALSE, "invalid \"-object\" option"); }
    
    /* Write it out: */
    test_voxm_tomogram_write("-", &a);
    return 0;
  }
    
void test_voxm_tomogram_write(char *outFile, ppv_array_t *a)
  {
    FILE *wr = open_write(outFile, TRUE);

    bool_t plain = FALSE;
    ppv_write_array(wr, a, plain);
    
    fclose(wr); 
  }

test_voxm_options_t *test_voxm_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    test_voxm_options_t *o = (test_voxm_options_t *)malloc(sizeof(test_voxm_options_t)); 
    
    /* Parse keyword parameters: */
    
    /* Size of array: */
    argparser_get_keyword(pp, "-size");
    o->size = (int32_t)argparser_get_next_int(pp, test_voxm_size_MIN, test_voxm_size_MAX);

    /* Object to create: */
    argparser_get_keyword(pp, "-object");
    o->object = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
