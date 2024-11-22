#define PROG_NAME "test_haf_read_obj"
#define PROG_DESC "tests the routines from {haf_read_obj.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-06-27 14:09:36 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2024  State University of Campinas (UNICAMP)\n\n" jslibs_copyright
  
#define PROG_AUTH \
  "J. Stolfi, 2024."
  
#define PROG_HIST \
  "Created by J. Stolfi in jun/2024."
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ {OBJ_NAME} ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads file \"in/{OBJ_NAME}.obj\".\n" \
  "\n" \
  "OPTIONS\n" \
  "  You have no choice in this matter.\n" \
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
  
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <argparser.h>
#include <vec.h>
#include <jsfile.h>
#include <r3.h>

#include <haf.h>
#include <obj_file.h>

#include <haf_read_obj.h>

typedef struct options_t
  { char *name;   /* Name of object file minus folder and extension. */
  } options_t;

int32_t main(int32_t argc, char **argv);

options_t *get_options (int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Parse the command line options: */
    options_t *o = get_options(argc, argv);
    
    /* Read the input file: */
    char *rd_fname = NULL; char *rd_fname = jsprintf("in/%s.obj", o->name);
    fprintf(stderr, "reading from file %s...\n", rd_fname);
    FILE *rd = open_read(rd_fname, TRUE);
    bool_t verbose = TRUE;
    int32_t nf = -1;
    r3_vec_t V = r3_vec_new(0);
    string_vec_t Vlab = string_vec_new(0);
    haf_arc_vec_t A = haf_arc_vec_new(0);
    int32_t *fleft = NULL;
    int32_t *vorg = NULL;
    haf_edge_id_t eid0 = 4615;
    haf_read_obj_file(rd, eid0, &nf, &V, &Vlab, &A, &fleft, &vorg, verbose);
    
    fclose(rd);
    
    int32_t nv = V.ne;
    int32_t ne = A.ne;
    assert(Vlab.ne == nv);
    
    fprintf(stderr, "got %6d vertices\n", nv);
    fprintf(stderr, "got %6d edges\n", ne);
    fprintf(stderr, "got %6d faces\n", nf);
    
    haf_check_topology(ne, A.e, eid0, verbose);

    /* !!! Should do more tests !!! */
    
    return 0;
  }

options_t *get_options (int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Parse keyword parameters: */

    /* Parse positional arguments: */
    o->name = argparser_get_next(pp);
    
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
        
    return o;
  }
