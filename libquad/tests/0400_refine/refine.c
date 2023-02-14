#define PROG_NAME "refine"
#define PROG_DESC "Refinement of quad-edge maps by quad triangulation."
#define PROG_VERS "1.0"

/* Last edited on 2007-02-05 10:13:34 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2007  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "ALgorithms taken from {MakeShape.m3} by Rober M. Rosi"\
  " and J. Stolfi, ca. 1994."
  
#define PROG_HIST \
  "Converted to C by J. Stolfi in fev/2007."
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -input {NAME_IN} \\\n" \
  "    -order {ORDER} \\\n" \
  "    -output {NAME_OUT}"

#include <oct.h>
#include <bool.h>

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct options_t
  { char *input;
    uint order;
    char *output;
  } options_t;

oct_arc_t calc_triang(oct_arc_t a)
  { Enum(a, glue_patch);
    Enum(a, putwr);
    return Middle(a)
  }

static uint edge_count = 0;

oct_arc_t mk_edge(uint grid_order)
  { frgb_t color1 = (frgb_t){{1.00, 0.95, 0.70}};
    frgb_t color2 = (frgb_t){{1.00, 0.95, 0.70}};
    oct_arc_t a = make_patch(grid_order, edge_count);

    set_primal_properties(
      a, 
      /* vertexRadius */  0.03,
      /* vertexColor */   (frgb_t){{0.0, 0.0, 0.0}};
      /* edgeRadius */    0.01,
      /* edgeColor */     (frgb_t){{0.0, 0.0, 0.0}};
      /* faceColor */     color1,
      /* faceTransp */    (frgb_t){{0.8, 0.8, 0.8}};
    );

    set_primal_properties(
      oct_rot(a),
      /* vertexRadius */  0.02,
      /* vertexColor */   (frgb_t){{1.0, 0.2, 0.0}};
      /* edgeRadius */    0.01,
      /* edgeColor */     (frgb_t){{1.0, 0.2, 0.0}};
      /* faceColor */     color2,
      /* faceTransp */    (frgb_t){{0.8, 0.8, 0.8}};
    );

    vc = oct_org(middle(a));
    set_node_color(vc, (frgb_t){{0.5, 0.1, 0.0});
    set_vertex_radius(vc, 0.02)

    edge_count += 1; 
    return a;
  }
  
 int main(int argc, char **argv)
  { options_t *o = get_options(argc, argv);
    oct_arc_t m = read_map(o->input);
    write_map(o->output, m);
    return 0;
  }
    
oct_arc_t read_map(char *prefix)
  { char *filename = NULL;
    asprintf(&filename, "%s.qe", prefix);
    FILE *rd = open_read(filename, TRUE);
    oct_read(rd, m);
    fclose(rd);
    free(filename);
  }
    
void write_map(char *prefix, oct_arc_t a)
  { char *filename = NULL;
    asprintf(&filename, "%s.qe", prefix);
    FILE *wr = open_write(filename, TRUE);
    oct_write_map(wr, m, NULL);
    fclose(wr);
    free(filename);
  }

void putwr (oct_arc_t a)
  { fprintf(stderr, "edge: "); 
    oct_write_arc(stderr, a, 3);
    fprintf(stderr, "  onext: "); 
    oct_write_arc(stderr, oct_onext(a), 3);
    fprintf(stderr, "\n");
  }

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a quad-edge data structure that" \
  " represents the topology of a 2D map, and writes out" \
  " a triangulation that is a refinement of that map.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -input {NAME_IN}\n" \
  "    Specifies the name of the input data file" \
  " (minus the \".qe\" extension).\n" \
  "\n" \
  "  -order {ORDER}\n" \
  "    Specifies the linear refinement factor to apply.\n" \
  "\n" \
  "  -output {NAME_OUT}\n" \
  "    Specifies the name of the output data file" \
  " (minus the \".qe\" extension).\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ls(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  " PROG_AUTH ".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  " PROG_HIST ".\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

options_t *get_options (int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-input");  
    o->input = argparser_get_next(pp);

    argparser_get_keyword(pp, "-order");                               
    o->grid_order = ParseParams.GetNextInt(1, 20); 

    argparser_get_keyword(pp, "-output");  
    o->output = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
        
    return o;
  }
  
