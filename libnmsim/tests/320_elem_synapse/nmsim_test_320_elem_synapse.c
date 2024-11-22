#define PROG_NAME "nmsim_test_320_elem_synapse"
#define PROG_DESC "basic tests of {limnmism} individual synapse attributes"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:33:52 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-06"
  
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
#include <jsmath.h>

#include <nmsim_elem_synapse.h>

void nmsim_elem_synapse_test_write(char *fname, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn);
  /* Writes the synapse attributes record {syn} to file {fname},
    assuming index {ise}. */

void nmsim_elem_synapse_test_read(char *fname, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn);
  /* Reads a synapse index and attributes record from file {fname}, and compares
    them to {ise} and {syn}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_group_synapse_count_t nsg = 50;  /* Assumed number of synapse classes. */
    nmsim_elem_neuron_count_t nne = 500; /* Assumed number of neurons. */
    
    for (int32_t i = 0; i < 10; i++)
      { nmsim_group_synapse_ix_t isg_max = (nmsim_group_synapse_ix_t)imin(3*i, nsg-1); 
        nmsim_group_synapse_ix_t isg_min = (nmsim_group_synapse_ix_t)(i < 3 ? 0 : isg_max); 
        nmsim_elem_neuron_ix_t ine_max = nne-1;
        nmsim_elem_neuron_ix_t ine_min = 0;
        double W_avg = -5.0;
        double W_dev = 2.0;
        nmsim_elem_synapse_t syn = nmsim_elem_synapse_throw
          ( isg_min, isg_max, ine_min, ine_max, ine_min, ine_max, W_avg, W_dev );
        nmsim_elem_synapse_show(stderr, "synapse = ", &syn, "\n");
        char *fname = jsprintf("out/test_synapse_elem_%03d.txt", i);
        nmsim_elem_synapse_ix_t ise = 17;
        nmsim_elem_synapse_test_write(fname, ise, &syn);
        nmsim_elem_synapse_test_read(fname, ise, &syn);
        free(fname);
      }
    return 0;
  }
  
void nmsim_elem_synapse_test_write(char *fname, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_elem_synapse_write(wr, ise, syn);
    fprintf(wr, "\n");
    fclose(wr);
  }

void nmsim_elem_synapse_test_read(char *fname, nmsim_elem_synapse_ix_t ise, nmsim_elem_synapse_t *syn)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_group_synapse_ix_t isg_max = syn->isg;
    nmsim_elem_neuron_ix_t ine_max = (nmsim_elem_neuron_ix_t)imax(syn->ine_pre, syn->ine_pos);
    nmsim_elem_synapse_t syn_read = nmsim_elem_synapse_read(rd, ise, isg_max, ine_max);
    fclose(rd);
    
    nmsim_elem_synapse_compare(&syn_read, syn);
  }

