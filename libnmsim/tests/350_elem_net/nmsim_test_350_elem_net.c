#define PROG_NAME "nmsim_test_350_elem_net"
#define PROG_DESC "basic tests of {limnmism} element-level network procedures"
#define PROG_VERS "1.0"

/* Last edited on 2019-06-18 11:30:40 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-05"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
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
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>

void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne_g_min,
    int32_t nne_g_max
  );
  /* Tests a net neuron with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups and {nsg} neuron groups, all generated at random. 
    
    If {nne_g_max > 0}, the neuron groups will be expanded to have
    between {nne_g_min} and {nne_g_max} neurons each; and each synapse group will be expanded
    into about {K*nne_g_max} synapses between those neurons, where {K}
    is the average in-degree of the post-synaptic neurons, as
    specified in the attributes of the synaptic group. See
    {nmsim_elem_net_expand_groups} for details. */

void nmsim_elem_net_test_write(char *fname, nmsim_elem_net_t *enet, double timeStep);
  /* Tests {nmsim_elem_net_write} with the elem-level network description
    {enet}, creating a file {fname}. */
    
void nmsim_elem_net_test_read(char *fname, nmsim_elem_net_t *enet, double timeStep);
  /* Reads a elem-level network description from file
    {fname}, and compares it to the given {enet}.*/

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_elem_net_test(0,0, 0,0, 0,0);
    nmsim_elem_net_test(2,3, 0,0, 0,0);
    nmsim_elem_net_test(2,3, 4,5, 0,0);
    nmsim_elem_net_test(5,10, 5,7, 4,6);
    nmsim_elem_net_test(5,10, 50,100, 10,20);
    return 0;
  }
    
void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne_g_min, 
    int32_t nne_g_max
  )
  {
    /* Assemble the file name {fname}: */
    char *fname = NULL;
    asprintf
      ( &fname, 
        "out/test_elem_net_%04dnc_%04dsc_%04dng_%04dsg_%04d-%04dnpg.txt", 
        nnc, nsc, nng, nsg, nne_g_min, nne_g_max
      );

    /* Create a random network: */
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    double K_max = 0.8*((double)nne_g_max);
    double K_min = 0.5*K_max;
    nmsim_group_net_t *gnet = 
      nmsim_group_net_throw(cnet, nng, nsg, nne_g_min, nne_g_max, K_min, K_max);
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
         
