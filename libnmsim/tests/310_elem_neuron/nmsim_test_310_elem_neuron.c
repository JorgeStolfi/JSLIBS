#define PROG_NAME "nmsim_test_310_elem_neuron"
#define PROG_DESC "basic tests of {limnmism} individual neuron attributes"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-11 18:46:25 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-11"
  
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
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <nmsim_elem_neuron.h>

void nmsim_elem_neuron_test_write(char *fname, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu);
  /* Writes the neuron attributes record {neu} to file {fname},
    assuming index {ine}. */

void nmsim_elem_neuron_test_read(char *fname, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu);
  /* Reads a neuron index and attributes record from file {fname}, and compares
    them to {ine} and {neu}. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_group_neuron_count_t nng = 418;
    for(int32_t i = 0; i < 10; i++)
      { nmsim_elem_neuron_t neu = nmsim_elem_neuron_throw(nng - 1);
        nmsim_elem_neuron_show(stderr, "neuron = ", &neu, "\n");
        char *fname = NULL;
        asprintf(&fname, "out/test_neuron_elem_%03d.txt", i);
        nmsim_elem_neuron_ix_t ine = 17;
        nmsim_elem_neuron_test_write(fname, ine, &neu);
        nmsim_elem_neuron_test_read(fname, ine, &neu);
        free(fname);
      }
    return 0;
  }
  
void nmsim_elem_neuron_test_write(char *fname, nmsim_elem_neuron_ix_t ine, nmsim_elem_neuron_t *neu)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_elem_neuron_write(wr, ine, neu);
    fprintf(wr, "\n");
    fclose(wr);
  }

void nmsim_elem_neuron_test_read
  ( char *fname, 
    nmsim_elem_neuron_ix_t ine, 
    nmsim_elem_neuron_t *neu
  )
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_group_neuron_ix_t ing_max = neu->ing; 
    nmsim_elem_synapse_count_t nse_out_max = neu->nse_out;
    nmsim_elem_neuron_t neu_read = nmsim_elem_neuron_read(rd, ine, ing_max,nse_out_max);
    fclose(rd);
    
    nmsim_elem_neuron_compare(&neu_read, neu);
  }

