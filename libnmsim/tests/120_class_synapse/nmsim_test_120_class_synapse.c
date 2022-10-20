#define PROG_NAME "nmsim_test_007_synapse_class"
#define PROG_DESC "basic tests of {limnmism} synapse class procedures"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:34:24 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-05"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <nmsim_class_synapse.h>

nmsim_class_synapse_t *nmsim_class_synapse_test_new(void);
  /* Creates an sample {nmsim_class_synapse_t} record. */

void nmsim_class_synapse_test_write(char *fname, nmsim_class_synapse_t *sclass);
  /* Writes the parameters record {sclass} of a synapse class to file {fname}. */

void nmsim_class_synapse_test_read(char *fname, nmsim_class_synapse_t *sclass);
  /* Reads the parameters record of a synapse class from file {fname},
    and compares it to {sclass}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_class_synapse_t *sclass = nmsim_class_synapse_throw();
    nmsim_class_synapse_show(stderr, "class = ", sclass, "\n");
    char *fname = "out/test_synapse_class.txt";
    nmsim_class_synapse_test_write(fname, sclass);
    nmsim_class_synapse_test_read(fname, sclass);
    nmsim_class_synapse_free(sclass);
    return 0;
  }
  
void nmsim_class_synapse_test_write(char *fname, nmsim_class_synapse_t *sclass)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_class_synapse_write(wr, sclass);
    fclose(wr);
  }

void nmsim_class_synapse_test_read(char *fname, nmsim_class_synapse_t *sclass)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_class_synapse_t *sclass_read = nmsim_class_synapse_read(rd);
    fclose(rd);
    nmsim_class_synapse_compare(sclass_read, sclass);
    nmsim_class_synapse_free(sclass_read);
  }

