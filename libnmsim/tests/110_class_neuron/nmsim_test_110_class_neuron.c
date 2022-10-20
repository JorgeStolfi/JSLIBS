#define PROG_NAME "nmsim_test_110_class_neuron"
#define PROG_DESC "basic tests of {limnmism} neuron class tools"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:34:30 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-03"
  
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

#include <nmsim_class_neuron.h>

nmsim_class_neuron_t *nmsim_class_neuron_test_new(void);
  /* Creates an sample {nmsim_class_neuron_t} record. */

void nmsim_class_neuron_test_write(char *fname, nmsim_class_neuron_t *nclass, double timeStep);
  /* Writes the neuron class parameters record {nclass} to file {fname},
    assuming the given {timeStep}. */

void nmsim_class_neuron_test_read(char *fname, nmsim_class_neuron_t *nclass, double timeStep);
  /* Reads a neuron class parameters record from file {fname},
    assuming the given {timeStep}, and compares it to {nclass}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    double timeStep = 1.0; /* Simulation time step (ms). */
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_throw();
    nmsim_class_neuron_show(stderr, "class = ", nclass, "\n");
    char *fname = "out/test_neuron_nclass.txt";
    nmsim_class_neuron_test_write(fname, nclass, timeStep);
    nmsim_class_neuron_test_read(fname, nclass, timeStep);
    nmsim_class_neuron_free(nclass);
    return 0;
  }
  
void nmsim_class_neuron_test_write(char *fname, nmsim_class_neuron_t *nclass, double timeStep)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_class_neuron_write(wr, nclass, timeStep);
    fclose(wr);
  }

void nmsim_class_neuron_test_read(char *fname, nmsim_class_neuron_t *nclass, double timeStep)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_class_neuron_t *nclass_read = nmsim_class_neuron_read(rd, timeStep);
    fclose(rd);
    nmsim_class_neuron_compare(nclass_read, nclass);
    nmsim_class_neuron_free(nclass_read);
  }

