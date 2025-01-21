#define PROG_NAME "test_pst_img_graph"
#define PROG_DESC "checks the {pst_img_graph.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-09 18:09:50 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_pst_img_graph_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -function {ZFUNC} \\\n" \
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
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2024-12-23 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_pst_img_graph_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {pst_img_graph}."

#define PROG_INFO_OPTS \
  ""

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <float_image.h>
#include <float_image_read_gen.h>
#include <image_file_format.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>

#include <pst_img_graph.h>
#include <pst_img_graph_from_maps.h>

typedef struct options_t
  { int32_t DUMMY;  /* Number of function to use. */
  } options_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

pst_img_graph_t* make_graph(bool_t add_diags);
  /* Creates a graph from the slope map file "in/IG.fni" and the
    weight map file "in/IW.fni". */

void write_graph(pst_img_graph_t *g);
  /* Writes {g} to file "out/graph_g.txt". */
  
pst_img_graph_t* read_graph(void);
  /* Reads back a graph from file "out/graph_g.txt". */
  
float_image_t *read_image(char *fname);
  /* Reads a FNI image from file "{fname}". */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test graph: */
    pst_img_graph_t *g = make_graph(FALSE);
    
    write_graph(g);
    pst_img_graph_t *h = read_graph();
    pst_img_graph_equal(g, h);

    /* Cleanup: */
    pst_img_graph_free(g);
    free(o); o = NULL;
    return 0;
  }
    
pst_img_graph_t* make_graph(bool_t add_diags)
  { 
    float_image_t *IG = read_image("in/IG.fni");
    float_image_t *IW = read_image("in/IW.fni");
    pst_img_graph_t *g = pst_img_graph_from_gradient_and_weight_maps(IG, IW, add_diags);
    float_image_free(IG);
    float_image_free(IW);
    return g;
  }
  
float_image_t *read_image(char *fname)
  { uint16_t *maxval = NULL;
    double expoDec, bias;
    bool_t verbose = TRUE;
    float_image_t *img = float_image_read_gen_named
      ( fname, image_file_format_FNI, 0, 1, &maxval, &expoDec, &bias, verbose );
    free(maxval);
    return img;
  }
void write_graph(pst_img_graph_t *g)
  { 
    FILE *wr = open_write("out/graph_g.txt", TRUE);
    pst_img_graph_write(wr, g);
    fclose(wr);
  }
    
pst_img_graph_t* read_graph(void)
  { 
    FILE *rd = open_read("out/graph_g.txt", TRUE);
    pst_img_graph_t *g = pst_img_graph_read(rd);
    fclose(rd);
    return g;
  }

options_t *tb_parse_options(int32_t argc, char **argv)
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

    /* PARSE KEYWORD ARGUMENTS: */
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
