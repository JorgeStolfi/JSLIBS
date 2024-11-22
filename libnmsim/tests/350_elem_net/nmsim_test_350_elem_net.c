#define PROG_NAME "nmsim_test_350_elem_net"
#define PROG_DESC "basic tests of {limnmism} element-level network procedures"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 05:29:30 by stolfi */ 

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
#include <vec.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_class_net_throw.h>

#include <nmsim_group_net.h>
#include <nmsim_group_net_throw.h>

#include <nmsim_elem_net.h>
#include <nmsim_elem_net_throw.h>

void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne, 
    int32_t nse
  );
  /* Tests a net neuron with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups, {nsg} neuron groups, {nne} neurons, and {nse}
    synapses, all generated at random. 
    
    The network must have at least one neuron, therefore {nnc,nng,nne}
    must be positive. If it has some synapses ({nse > 0}), it must have
    at least one synaptic bundle ({nsg > 0}). If it has a synaptic
    bundle ({nsg > 0}), it must have at least one synapse class ({nsc >
    0}). */

void nmsim_elem_net_test_write(char *fname, nmsim_elem_net_t *enet, double timeStep);
  /* Tests {nmsim_elem_net_write} with the elem-level network description
    {enet}, creating a file {fname}. */
    
void nmsim_elem_net_test_read(char *fname, nmsim_elem_net_t *enet, double timeStep);
  /* Reads a elem-level network description from file
    {fname}, and compares it to the given {enet}.*/

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_elem_net_test(1,0, 1,0, 1,0);
    nmsim_elem_net_test(2,3, 1,0, 1,0);
    nmsim_elem_net_test(2,3, 3,5, 3,0);
    nmsim_elem_net_test(5,10, 5,7, 4,6);
    nmsim_elem_net_test(5,10, 50,100, 10,20);
    return 0;
  }
    
void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne, 
    int32_t nse
  )
  {
    /* Assemble the file name {fname}: */
    char *fname = jsprintf(
        "out/test_elem_net_%04dnc_%04dsc_%04dng_%04dsg_%04dne_%04dse.txt", 
        nnc, nsc, nng, nsg, nne, nse
      );

    /* Create a random network: */
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    nmsim_group_net_t *gnet = nmsim_group_net_throw(cnet, nng, nsg, nne, nse);
    nmsim_elem_net_t *enet = nmsim_elem_net_throw(gnet);
    
    /* Test write and read: */
    double timeStep = 1.0;
    nmsim_elem_net_test_write(fname, enet, timeStep);
    nmsim_elem_net_test_read(fname, enet, timeStep);
    free(fname);
  }

void nmsim_elem_net_test_write(char *fname, nmsim_elem_net_t *enet, double timeStep)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_elem_net_write(wr, enet, timeStep);
    fclose(wr);
  }

void nmsim_elem_net_test_read(char *fname, nmsim_elem_net_t *enet, double timeStep)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_elem_net_t *enet_read = nmsim_elem_net_read(rd,timeStep);
    fclose(rd);

    nmsim_elem_net_compare(enet_read, enet);
    nmsim_elem_net_free(enet_read);
  }
         
