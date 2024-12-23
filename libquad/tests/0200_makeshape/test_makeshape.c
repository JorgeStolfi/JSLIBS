#define PROG_NAME "makeshape"
#define PROG_DESC "creates some 2D maps using the quad-edge structure"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-22 11:22:01 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2007  State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#define PROG_AUTH \
  "Modula-3 version {MakeShape.m3} created by Rober M. Rosi"\
  " and J. Stolfi, ca. 1994."
  
#define PROG_HIST \
  "Converted to C by J. Stolfi in fev/2007."
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -shape { torus | klein | ... | star }"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes to standard output a quad-edge description of\n" \
  " a simple map selected by name.\n" \
  " to disk.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -shape {NAME}\n" \
  "    Specifies which basic shape is to be built.  The options are:\n" \
  "      " SHAPE_NAMES "\n" \
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
  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <oct.h>
#include <oct_shapes.h>
#include <argparser.h>
#include <vec.h>

/* INTERNAL DEFINITIONS AND PROTOTYPES */

#define SHAPE_NAMES \
  "torus bitorus tritorus klein klein2 klein3 " \
  "projective tetra stick ring cube sausage orange fork star"
  
typedef struct options_t
  { char *shape;  /* Name of map to build. */
  } options_t;

int32_t main(int32_t argc, char **argv);

string_vec_t split_shape_names(char *names);
  /* Splits a single string {names}, containing shape names
    separated by by blanks, into a {string_vec_t} of separate names. */

oct_arc_t make_map(char *name);
  /* Returns the oct-edge representation of the map called {name}. */

void write_map(oct_arc_t a);

options_t *get_options (int32_t argc, char **argv, string_vec_t *shape_name);

uint32_t get_shape_num(char *shp, string_vec_t *shape_name);
  /* Returns the index {i} in {0..NS-1} such that
    {shape_name->e[i]===shp}, or {NS} if not there; where
    {NS==shape_name->ne}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Prepare the table of shape names: */
    string_vec_t shape_name = split_shape_names(SHAPE_NAMES);
    
    /* Parse the command line options: */
    options_t *o = get_options(argc, argv, &shape_name);
    
    /* Build the map and check number-name consistency: */
    oct_arc_t m = make_map(o->shape);
    
    /* Write the map to standard output: */
    write_map(m);
    return 0;
  }

#define eq(X,Y) (0 == strcmp((X), (Y)))

oct_arc_t make_map(char *name)
  { 
    if (eq(name, "torus"))       { return make_torus();     } else
    if (eq(name, "bitorus"))     { return make_bitorus();   } else
    if (eq(name, "tritorus"))    { return make_tritorus();  } else
    if (eq(name, "klein"))       { return make_klein();     } else
    if (eq(name, "klein2"))      { return make_klein2();    } else
    if (eq(name, "klein3"))      { return make_klein3();    } else
    if (eq(name, "projective"))  { return make_projetive(); } else
    if (eq(name, "tetra"))       { return make_tetra();     } else
    if (eq(name, "stick"))       { return make_stick();     } else
    if (eq(name, "ring"))        { return make_ring(5);     } else
    if (eq(name, "cube"))        { return make_cube();      } else
    if (eq(name, "sausage"))     { return make_sausage(3);  } else
    if (eq(name, "orange"))      { return make_fork(7, 0);  } else
    if (eq(name, "fork"))        { return make_fork(3, 2);  } else
    if (eq(name, "star"))        { return make_fork(5, 1);  } else
      { demand(FALSE, "invalid shape number"); return oct_arc_NULL; }
  } 
    
void write_map(oct_arc_t a)
  { oct_arc_vec_t root = oct_arc_vec_new(1);
    root.e[0] = a;
    oct_write_map(stdout, &root, NULL);
    fflush(stdout);
  }

options_t *get_options (int32_t argc, char **argv, string_vec_t *shape_name)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-shape");  
    o->shape = argparser_get_next(pp);
    uint32_t num = get_shape_num(o->shape, shape_name);
    if (num >= shape_name->ne) 
      { argparser_error(pp, "invalid shape name"); }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
        
    return o;
  }

uint32_t get_shape_num(char *shp, string_vec_t *shape_name)
  {
    for (uint32_t i = 0;  i < shape_name->ne; i++)
      { if (0 == strcmp(shp, shape_name->e[i])) { return i; } }
    return shape_name->ne;
  }

string_vec_t split_shape_names(char *names)
  {
    string_vec_t shape_names = string_vec_new(0);
    uint32_t NS = 0;
    char *p = names;
    while (TRUE)
      { while ((*p) == ' ') { p++; } 
        if ((*p) == 0) { break; }
        char *q = p;
        while (((*q) != 0) && ((*q) != ' ')) { q++; }
        uint32_t nc = (uint32_t)(q - p);
        char *s = notnull(malloc(nc + 1), "no mem");
        (void)strncpy(s, p, nc);
        s[nc] = 0;
        string_vec_expand(&shape_names, (vec_index_t)NS);
        shape_names.e[NS] = s; 
        NS++;
        p = q;
      }
     string_vec_trim(&shape_names, NS);
     return shape_names;
  }
