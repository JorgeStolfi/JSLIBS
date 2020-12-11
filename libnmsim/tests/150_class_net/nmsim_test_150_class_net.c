#define PROG_NAME "nmsim_test_150_class_net"
#define PROG_DESC "basic tests of {limnmism} class-level network description"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-11 13:34:31 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-10"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_class_net_throw.h>

void nmsim_class_net_test(int32_t nnc, int32_t nsc);
  /* Tests creation, reading, and writing of a class-level
    network description with {nnc} neuron classes and {nsc} synapse 
    classes. The parameter {nnc} must be positive. */

void nmsim_class_net_test_write(char *fname, nmsim_class_net_t *cnet, double timeStep);
  /* Tests {nmsim_class_net_write} with the class-level network description
    {cnet}, creating a file {fname}. */

void nmsim_class_net_test_read(char *fname, nmsim_class_net_t *cnet, double timeStep);
  /* Reads a class-level network description from file
    {fname}, and compares it to the given {cnet}.*/

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_class_net_test(1,0);
    nmsim_class_net_test(1,1);
    nmsim_class_net_test(2,3);
    return 0;
  }
    
void nmsim_class_net_test(int32_t nnc, int32_t nsc)
  {
    /* Create the network: */
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    /* Test write and read: */
    char *fname = NULL;
    asprintf(&fname, "out/test_class_net_%04dnc_%04dsc.txt", nnc, nsc);
    double timeStep = 1.0;
    nmsim_class_net_test_write(fname, cnet, timeStep);
    nmsim_class_net_test_read(fname, cnet, timeStep);
    free(fname);
  }

void nmsim_class_net_test_write(char *fname, nmsim_class_net_t *cnet, double timeStep)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_class_net_write(wr, cnet, timeStep);
    fclose(wr);
  }

void nmsim_class_net_test_read(char *fname, nmsim_class_net_t *cnet, double timeStep)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_class_net_t *cnet_read = nmsim_class_net_read(rd,timeStep);
    fclose(rd);

    nmsim_class_net_compare(cnet_read, cnet);
    nmsim_class_net_free(cnet_read);
  }

