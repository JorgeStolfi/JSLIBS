#define PROG_NAME "nmsim_test_360_elem_net_group_stats"
#define PROG_DESC "basic tests of {limnmism} synapse group statistics tols"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:49:43 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2020  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2020-12-13"
  
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
#include <jsmath.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_class_net.h>
#include <nmsim_class_net_throw.h>
#include <nmsim_group_net.h>
#include <nmsim_group_net_throw.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_throw.h>

#include <nmsim_elem_net_group_stats.h>

void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne, 
    int32_t nse
  );
  /* Tests synapse group statistics on an element-level
    network with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups, {nsg} neuron groups, {nne} neurons, and {nse}
    synapses, all generated at random. 
    
    The network must have at least one neuron, therefore {nnc,nng,nne}
    must be positive. If it has some synapses ({nse > 0}), it must have
    at least one synaptic bundle ({nsg > 0}). If it has a synaptic
    bundle ({nsg > 0}), it must have at least one synapse class ({nsc >
    0}).
    
    The statistics are written to a file called "out/test_stats_{STUFF}.txt". */

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
    /* Create a random network: */
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    nmsim_group_net_t *gnet = nmsim_group_net_throw(cnet, nng, nsg, nne, nse);
    nmsim_elem_net_t *enet = nmsim_elem_net_throw(gnet);
    
    /* Assemble the file name {fname}: */
    char *prefix = NULL;
    asprintf
      ( &prefix, 
        "out/test_stats_%04dnc_%04dsc_%04dng_%04dsg_%04dne_%04dse", 
        nnc, nsc, nng, nsg, nne, nse
      );

    /* Write out the network for reference: */
    char *fname_net = NULL;
    asprintf(&fname_net, "%s_net.txt", prefix);
    FILE *wr_net = open_write(fname_net, TRUE);
    double timeStep = 1.0;
    nmsim_elem_net_write(wr_net, enet, timeStep);
    fclose(wr_net);
    free(fname_net);
    
    /* Print the synapse group statistics: */
    char *fname_stats = NULL;
    asprintf(&fname_stats, "%s_stats.txt", prefix);
    FILE *wr_stats = open_write(fname_stats, TRUE);
    nmsim_elem_net_group_stats_print_all(wr_stats, enet, "  | ", " |\n");
    fclose(wr_stats);
    free(fname_stats);
  }

         
