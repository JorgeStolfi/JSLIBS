#define PROG_NAME "test_voxb"
#define PROG_DESC "Test the voxel-based 3D modeling routines of {voxb.h}"
#define PROG_VERS "1.0"

#define test_voxb_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-09 01:03:47 by jstolfi */

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
  "  The program writes to standard output a tomogram (3D voxel" \
  " array) containing a model of some test objects.\n" \
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
  "  " test_voxb_C_COPYRIGHT ".\n" \
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
#include <ppv_array_write.h>

#include <voxb_obj.h>
#include <voxb_splat.h>
#include <voxb_splat_tube.h>

#include <test_voxb_objs.h>
#include <test_voxb_tubes.h>
#include <test_voxb_mark.h>
#include <test_voxb_erolate.h>

/* COMMAND-LINE OPTIONS */

typedef struct test_voxb_options_t
  { int32_t size;  /* Size of voxel array. */
    char *object;  /* Name of test object. */
  } test_voxb_options_t;

#define test_voxb_size_MIN (100)
  /* Min tomogram size along any axis. */

#define test_voxb_size_MAX (3000)
  /* Max tomogram size along any axis. */

/* INTERNAL PROTOTYPES */

    /* The unit for all linear dimensions is the voxel side. */ 
    
test_voxb_options_t *test_voxb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {test_voxb_options_t}. */

void test_voxb_tomogram_write(char *outFile, ppv_array_t *A); 

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    test_voxb_options_t *o = test_voxb_parse_options(argc, argv);
    
    ppv_dim_t d = 3;
    
    /* Array dimensions: */
    int32_t NX = 3*o->size; 
    int32_t NY = o->size + o->size/2; 
    int32_t NZ = o->size; 
    
    /* Allocate the voxel array: */
    ppv_size_t sz[d];
    sz[0] = NZ;
    sz[1] = NY;
    sz[2] = NX;
    ppv_sample_t maxsmp = 1;
    ppv_array_t *A = ppv_array_new(d, sz, maxsmp);
    
    /* Center and radius of array: */
    r3_t ctr = (r3_t){{ 0.5*NX, 0.5*NY, 0.5*NZ }};
    r3_t rad = (r3_t){{ 0.5*NX, 0.5*NY, 0.5*NZ }};

    /* Mark the corners: */
    test_voxb_mark_corners(A, &ctr, &rad);
    test_voxb_mark_edges(A, &ctr, &rad, 0);
    test_voxb_mark_edges(A, &ctr, &rad, 1);
    test_voxb_mark_edges(A, &ctr, &rad, 2);

    if (strcmp(o->object, "objs") == 0)
      { test_voxb_objs(A, &ctr, &rad); }
    else if (strcmp(o->object, "tubes") == 0)
      { test_voxb_tubes(A, &ctr, &rad); }
    else  if (strcmp(o->object, "erolate") == 0)
      { double smr = 3.5;  /* Smoothing radius (vx). */
        test_voxb_erolate(A, &ctr, &rad, smr);
      }
    else
      { demand(FALSE, "invalid \"-object\" option"); }
    
    /* Write it out: */
    test_voxb_tomogram_write("-", A);
    return 0;
  }
    
void test_voxb_tomogram_write(char *outFile, ppv_array_t *A)
  {
    FILE *wr = open_write(outFile, TRUE);

    bool_t plain = FALSE;
    ppv_array_write_file(wr, A, plain);
    
    fclose(wr); 
  }

test_voxb_options_t *test_voxb_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    test_voxb_options_t *o = (test_voxb_options_t *)malloc(sizeof(test_voxb_options_t)); 
    
    /* Parse keyword parameters: */
    
    /* Size of array: */
    argparser_get_keyword(pp, "-size");
    o->size = (int32_t)argparser_get_next_int(pp, test_voxb_size_MIN, test_voxb_size_MAX);

    /* Object to create: */
    argparser_get_keyword(pp, "-object");
    o->object = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
