#define PROG_NAME "nmsim_test_220_group_synapse"
#define PROG_DESC "basic tests of {limnmism} synaptic band attributes"
#define PROG_VERS "1.0"

/* Last edited on 2019-02-28 20:17:18 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-06"
  
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
#include <jsmath.h>

#include <nmsim_group_synapse.h>

void nmsim_group_synapse_test_write(char *fname, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp);
  /* Writes the synapse group attributes {*sgrp} to file {fname},
    assuming index {isg}. */

void nmsim_group_synapse_test_read(char *fname, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp);
  /* Reads a synapse group index and attributes record from file {fname}, and compares
    them to {isg} and {*sgrp}. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_class_neuron_count_t nsc = 50;  /* Assumed number of synapse classes. */
    nmsim_group_neuron_count_t nng = 500; /* Assumed number of neuron groups. */
    nmsim_group_neuron_ix_t ing_pre = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
    nmsim_group_neuron_ix_t ing_pos = (nmsim_group_neuron_ix_t)int64_abrandom(0, nng-1);
    double K_min = 5.0;
    double K_max = 10.0;
    nmsim_class_synapse_ix_t isc_max = nsc - 1; /* Max class index. */
    nmsim_group_synapse_t sgrp = nmsim_group_synapse_throw(isc_max, ing_pre, ing_pos, K_min, K_max);
    nmsim_group_synapse_show(stderr, "synaps group = ", &sgrp, "\n");
    char *fname = "out/test_synapse_group.txt";
    nmsim_group_synapse_ix_t isg = 17;
    nmsim_group_synapse_test_write(fname, isg, &sgrp);
    nmsim_group_synapse_test_read(fname, isg, &sgrp);
    return 0;
  }
  
void nmsim_group_synapse_test_write(char *fname, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp)
  { FILE *wr = open_write(fname, TRUE);
    nmsim_group_synapse_write(wr, isg, sgrp);
    fprintf(wr, "\n");
    fclose(wr);
  }

void nmsim_group_synapse_test_read(char *fname, nmsim_group_synapse_ix_t isg, nmsim_group_synapse_t *sgrp)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_class_synapse_ix_t isc_max = sgrp->isc;
    nmsim_group_neuron_ix_t ing_max = (nmsim_group_neuron_ix_t)imax(sgrp->ing_pre, sgrp->ing_pos);
    nmsim_group_synapse_t sgrp_read = nmsim_group_synapse_read(rd, isg, isc_max, ing_max);
    fclose(rd);
    
    nmsim_group_synapse_compare(&sgrp_read, sgrp);
  }

