#define PROG_NAME "nmsim_test_210_group_neuron"
#define PROG_DESC "basic tests of {limnmism} neuron population attributes"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:34:10 by stolfi */ 

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

#include <nmsim_group_neuron.h>

void nmsim_group_neuron_test_write(char *fname, nmsim_group_neuron_ix_t ing, nmsim_group_neuron_t *ngrp);
  /* Writes the neuron population attributes record {ngrp} to file {fname},
    assuming index {ing}. */

void nmsim_group_neuron_test_read
  ( char *fname, 
    nmsim_group_neuron_ix_t ing,
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_group_synapse_count_t nsg_out_max,
    nmsim_group_neuron_t *ngrp
  );
  /* Reads a neuron population index and attributes record from file {fname}, and compares
    them to {ing} and {ngrp}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_class_neuron_ix_t inc_max = 418;
    nmsim_elem_neuron_count_t nne_g_min = 1;
    nmsim_elem_neuron_count_t nne_g_max = 20;
    nmsim_group_synapse_count_t nsg_out_min = 0;
    nmsim_group_synapse_count_t nsg_out_max = 30;
    for (int32_t i = 0; i < 10;i++)
      { nmsim_group_neuron_t ngrp = nmsim_group_neuron_throw
          ( inc_max, nne_g_min, nne_g_max, nsg_out_min, nsg_out_max);
        char *fname = NULL;
        asprintf(&fname, "out/test_neuron_group_%03d.txt", i);
        nmsim_group_neuron_show(stderr, "neuron group = ", &ngrp, "\n");
        nmsim_group_neuron_ix_t ing = 17;
        nmsim_group_neuron_test_write(fname, ing, &ngrp);
        nmsim_group_neuron_test_read(fname, ing, nne_g_max, nsg_out_max, &ngrp);
        free(fname);
      }
    return 0;
  }
  
void nmsim_group_neuron_test_write(char *fname, nmsim_group_neuron_ix_t ing, nmsim_group_neuron_t *ngrp)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_group_neuron_write(wr, ing, ngrp);
    fprintf(wr, "\n");
    fclose(wr);
  }

void nmsim_group_neuron_test_read
  ( char *fname, 
    nmsim_group_neuron_ix_t ing,
    nmsim_elem_neuron_count_t nne_g_max,
    nmsim_group_synapse_count_t nsg_out_max,
    nmsim_group_neuron_t *ngrp
  )
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_class_neuron_ix_t inc_max = ngrp->inc;
    nmsim_group_neuron_t ngrp_read = nmsim_group_neuron_read
      ( rd, ing, inc_max, nne_g_max, nsg_out_max );
    fclose(rd);
    
    nmsim_group_neuron_compare(&ngrp_read, ngrp);
  }

