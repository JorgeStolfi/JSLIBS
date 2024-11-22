#define PROG_NAME "test_obj_file"
#define PROG_DESC "tests the routines from {obj_file_read.h} and {obj_file_write.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-06-27 18:08:54 by stolfi */ 

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
  "  The program reads file \"in/{OBJ_NAME}.obj\" and writes out the same model to.\n" \
  " to \"out/{OBJ_NAME}-wr.obj\".\n" \
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

#include <obj_file.h>
#include <obj_file_read.h>
#include <obj_file_write.h>

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
    obj_file_data_t *D = obj_file_read(rd, verbose);
    fclose(rd);
    
    fprintf(stderr, "found %6d vertices\n", D->V.ne);
    fprintf(stderr, "found %6d texpoints\n", D->T.ne);
    fprintf(stderr, "found %6d normals\n", D->N.ne);
    fprintf(stderr, "found %6d faces\n", D->FV.ne);
    assert(D->VL.ne == D->V.ne);
    
    for (int32_t kf = 0; kf < D->FV.ne; kf++)
      { int32_vec_t *FVk = &(D->FV.e[kf]);
        int32_vec_t *FTk = &(D->FT.e[kf]);
        int32_vec_t *FNk = &(D->FN.e[kf]);
        int32_t nc = FVk->ne;
        fprintf(stderr, "  face %6d has %6d corners:", kf, nc);
        assert(FTk->ne == nc);
        assert(FNk->ne == nc);
        for (int32_t kc = 0; kc < nc; kc++)
          { assert(FVk->e[kc] >= 0);
            fprintf(stderr, " %d", FVk->e[kc] + 1);
            if ((FTk->e[kc] >= 0) || (FNk->e[kc] >= 0))
              { fprintf(stderr, "/");
                if (FTk->e[kc] >= 0) { fprintf(stderr, "%d", FTk->e[kc] + 1); }
                if (FNk->e[kc] >= 0)
                  { fprintf(stderr, "/");
                    fprintf(stderr, "%d", FTk->e[kc] + 1);
                  }
              }
          }
        fprintf(stderr, "\n");
      }

    /* Write it out: */
    char *wr_fname = NULL; char *wr_fname = jsprintf("out/%s-wr.obj", o->name);
    fprintf(stderr, "writing out to file %s...\n", wr_fname);
    FILE *wr = open_write(wr_fname, TRUE);
    int32_t prec = 4;
    obj_file_write(wr, D, prec);
    fclose(wr);
    
    /* Read it back: */
    char *rd2_fname = wr_fname;
    fprintf(stderr, "reading back from file %s...\n", rd2_fname);
    FILE *rd2 = open_read(rd2_fname, TRUE);
    obj_file_data_t *D2 = obj_file_read(rd2, verbose);
    fclose(rd2);
    
    /* Compare: */
    fprintf(stderr, "comparing readback with original...\n");
    double tol = pow(0.1, prec);
    bool_t eq = obj_file_data_compare(D, D2, tol, TRUE);
    if (eq)
      { fprintf(stderr, "compared equal!\n"); }
    else
      { demand(FALSE, "readback differs from original"); }
    
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
