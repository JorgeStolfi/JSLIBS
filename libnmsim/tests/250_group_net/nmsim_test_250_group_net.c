#define PROG_NAME "nmsim_test_250_group_net"
#define PROG_DESC "basic tests of {limnmism} group-level network description"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:33:59 by stolfi */ 

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

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_class_net_throw.h>

#include <nmsim_group_net.h>
#include <nmsim_group_net_throw.h>

void nmsim_group_net_test
  ( int32_t nnc, 
    int32_t nsc, 
    int32_t nng,
    int32_t nsg
  );
  /* Tests a neuron net with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups (populations), and {nsg} synapse groups (bundles). */

void nmsim_group_net_test_write(char *fname, nmsim_group_net_t *gnet, double timeStep);
  /* Tests {nmsim_group_net_write} with the grpup-level network description
    {gnet}, creating a file {fname}. */
    
void nmsim_group_net_test_read
  ( char *fname, 
    nmsim_group_net_t *gnet, 
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_elem_neuron_count_t nse_g_max, 
    double timeStep
  );
  /* Reads a group-level network description from file
    {fname}, and compares it to the given {gnet}.*/

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_group_net_test(1,0,1,0);
    nmsim_group_net_test(2,3,1,0);
    nmsim_group_net_test(2,3,1,1);
    nmsim_group_net_test(5,10,50,100);
    return 0;
  }
    
void nmsim_group_net_test
  ( int32_t nnc, 
    int32_t nsc, 
    int32_t nng,
    int32_t nsg
  )
  {
    assert((nnc >= 1) && (nng >= 1));
    assert((nsc >= 0) && (nsg >= 0));
    if (nsc == 0) { assert(nsg == 0); } 
    /* Create the network: */
    nmsim_elem_neuron_count_t nne_min = nng;
    nmsim_elem_neuron_count_t nne_max = 4*nng;
    nmsim_elem_synapse_count_t nse_min = 0;
    nmsim_elem_synapse_count_t nse_max = 10*nsg;
    for (int32_t i = 0; i < 10; i++) 
      { nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
        nmsim_elem_neuron_count_t nne = (nmsim_elem_neuron_count_t)int64_abrandom(nne_min, nne_max);
        nmsim_elem_synapse_count_t nse = (nmsim_elem_synapse_count_t)int64_abrandom(nse_min, nse_max);
        nmsim_group_net_t *gnet = nmsim_group_net_throw(cnet, nng, nsg, nne, nse);
        /* Test write and read: */
        char *fname = jsprintf("out/test_group_net_%03d_%04dnc_%04dsc_%06dng_%06dsg.txt", i, nnc, nsc, nng, nsg);
        double timeStep = 1.0;
        nmsim_group_net_test_write(fname, gnet, timeStep);
        nmsim_group_net_test_read(fname, gnet, nne_max, nse_max, timeStep);
        free(fname);
        nmsim_group_net_free(gnet);
      }
  }

void nmsim_group_net_test_write(char *fname, nmsim_group_net_t *gnet, double timeStep)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_group_net_write(wr, gnet, timeStep);
    fclose(wr);
  }

void nmsim_group_net_test_read
  ( char *fname, 
    nmsim_group_net_t *gnet, 
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_elem_neuron_count_t nse_g_max, 
    double timeStep
  )
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_group_net_t *gnet_read = nmsim_group_net_read(rd, nne_g_max, nse_g_max, timeStep);
    fclose(rd);

    nmsim_group_net_compare(gnet_read, gnet);
    nmsim_group_net_free(gnet_read);
  }
