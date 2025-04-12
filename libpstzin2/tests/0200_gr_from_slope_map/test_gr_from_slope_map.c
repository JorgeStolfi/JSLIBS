#define PROG_NAME "test_gr_from_maps"
#define PROG_DESC "checks the {pst_gr.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2025-03-15 14:14:31 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_gr_from_maps_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "  -function {NFUNC} {XFUNC} \\\n" \
  "  -noisy {NOISY} \\\n" \
  "  -size {NX} {NY} \\\n" \
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
  "  " test_gr_from_maps_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {pst_gr_from_slope_map}.\n" \
  "\n" \
  "  The program reads a file \"in/proc/{TEST}-G.fni\" where {TEST}" \
  " is \"{NFUNC}-{XFUNC}-{NOISY}-{NX}x{NY}\", {NFUNC} is formatted as \"%02d\", {NOISY}" \
  " is \"N\" or \"Y\", and {NX,NY} are each formatted as \"%04d\".\n" \
  "\n" \
  "  The file should contain a slope map with size {NX} by {NY}, with" \
  " 2 or 3 channels.  The program builds a graph from that slope map, then" \
  " writes it out to \"out/{TEST}-graph.txt\", and plots it" \
  " to file \"out/{TEST}-graph.eps\"."

#define PROG_INFO_OPTS \
  "  -function {NFUNC} {XFUNC}\n" \
  "    This argument specifies the numeric ID and name of the" \
  " procedural slope map. See {pst_proc_map.h}.\n" \
  "\n" \
  "  -noisy {NOISY}\n" \
  "    This argument specifies if the slope map is clean (\"N\") or" \
  " contaminated by noise (\"Y\").\n" \
  "\n" \
  "  -size {NX} {NY} \\\n" \
  "    This argument specifies the size (cols and rows)" \
  " of the slope map."

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_read_gen.h>
#include <image_file_format.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <argparser.h>

#include <pst_gr.h>
#include <pst_gr_read.h>
#include <pst_gr_write.h>
#include <pst_gr_plot.h>
#include <pst_gr_from_slope_map.h>

typedef struct options_t
  { uint32_t function_num;  /* Number of function to use. */
    char *function_name;    /* Number of function to use. */
    bool_t noisy;           /* True for contaminated version of the slope map. */
    uint32_t NX;  /* Cols in slope map. */
    uint32_t NY;  /* Rows in slope map. */
  } options_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */
  
float_image_t* read_slope_map(char *prefix);
  /* Reads a slope map from the file "in/proc/{prefix}-G.fni". 
    It should have 2 or 3 channels.  If the latter, channel 2
    is interprested as the pixel weight. */

pst_gr_t* make_graph(float_image_t *G, bool_t addDiags);
  /* Creates a graph from the slope map {G}. */

void write_graph(char *prefix, pst_gr_t *gr);
  /* Writes {gr} to file "out/{prefix}-graph.txt". */
 
void plot_graph(char *prefix, pst_gr_t *gr);
  /* Writes an Encapsulated Postscript plot of {gr} to
    file "out/{prefix}-graph.eps". */
 
pst_gr_t* read_graph(char *prefix);
  /* Reads back a graph from file "out/graph_g.txt". */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    char *prefix = jsprintf
      ( "%02d-%s-%c-%04dx%04d", 
        o->function_num, o->function_name, "NY"[o->noisy], o->NX, o->NY
      );
    
    fprintf(stderr, "reading slope map ...\n");
    float_image_t *G = read_slope_map(prefix);
    
    /* Create a test graph: */
    fprintf(stderr, "creating the graph ...\n");
    pst_gr_t *gr = make_graph(G, FALSE);

    fprintf(stderr, "validating the graph ...\n");
    pst_gr_check_consistency(gr);
    
    fprintf(stderr, "writing the graph ...\n");
    write_graph(prefix, gr);
    
    fprintf(stderr, "reading back the graph ...\n");
    pst_gr_t *hr = read_graph(prefix);
    
    fprintf(stderr, "comparing with original ...\n");
    pst_gr_equal(gr, hr);
    
    fprintf(stderr, "plotting ...\n");
    plot_graph(prefix, gr);

    /* Cleanup: */
    float_image_free(G);
    pst_gr_free(gr);
    free(o); o = NULL;
    return 0;
  }
    
pst_gr_t* make_graph(float_image_t *G, bool_t addDiags)
  { pst_gr_t *gr = pst_gr_from_slope_map(G, addDiags);
    return gr;
  }
  
float_image_t *read_slope_map(char *prefix)
  { char *fname = jsprintf("in/proc/%s-G.fni", prefix);
    uint16_t *maxval = NULL;
    double expoDec, bias;
    bool_t yUp = FALSE; /* Irrelevant. for FNI files */
    bool_t verbose = TRUE;
    float_image_t *G = float_image_read_gen_named
      ( fname, image_file_format_FNI, yUp, 0, 1, &maxval, &expoDec, &bias, verbose );
    free(maxval);
    free(fname);
    return G;
  }
  
void write_graph(char *prefix, pst_gr_t *gr)
  { char *fname = jsprintf("out/%s-graph.txt", prefix);
    bool_t verbose = TRUE;
    pst_gr_write_named(fname, gr, verbose);
    free(fname);
  }
    
pst_gr_t* read_graph(char *prefix)
  { char *fname = jsprintf("out/%s-graph.txt", prefix);
    bool_t verbose = TRUE;
    pst_gr_t *gr = pst_gr_read_named(fname, verbose);
    free(fname);
    return gr;
  }
  
void plot_graph(char *prefix, pst_gr_t *gr)
  { char *fname = jsprintf("out/%s-graph.eps", prefix);
    double fontSize = 0.0;     /* Label font size (pt) or 0 to omit labels. */
    double vertexRadius = 1.0; /* Radius of vertices (mm) */
    pst_gr_plot_named(fname, gr, fontSize, vertexRadius); 
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
    
    argparser_get_keyword(pp, "-function");
    o->function_num = (uint32_t)argparser_get_next_int(pp, 0, 99);
    o->function_name = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-noisy");
    o->noisy = argparser_get_next_bool(pp);
    
    argparser_get_keyword(pp, "-size");
    o->NX = (uint32_t)argparser_get_next_int(pp, 1, 4096);
    o->NY = (uint32_t)argparser_get_next_int(pp, 1, 4096);
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
