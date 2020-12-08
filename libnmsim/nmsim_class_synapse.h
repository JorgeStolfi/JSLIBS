#ifndef nmsim_class_synapse_H
#define nmsim_class_synapse_H

/* Common parameters for a class of synapses. */ 
/* Last edited on 2020-12-06 18:56:01 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_class_synapse_t
  { double W_avg; /* Average resting synapse strength. */
    double W_dev; /* Deviation of resting synapse strengths. */
  } nmsim_class_synapse_t;
  /* Parameters of a set of synapse.  See {nmsim_class_synapse_INFO} 
    for the meaning of the parameters. */
    
#define nmsim_class_synapse_INFO \
  "The parameters {W_avg} and {W_dev} are the /average/ and" \
  " /deviation/ of the /resting strength/ of the synapses in the" \
  " class.  Namely, if {W_dev} is zero, then the actual resting strength {W} of every" \
  " synapse of this class is {W_avg}.  Otherwise, the {W} values of" \
  " those synapses are samples from a random distribution with mean {W_avg} and deviation {W_dev}," \
  " such that {W} has the same sign as {W_avg} and {|W|} has a log-normal" \
  " distribution. Note that {W} always has the" \
  " same sign as {W_avg}, even when {W_dev} is greater than {|W_avg|}.\n" \
  "\n" \
  "  If {W_avg} and {W_dev} are zero, then all synapses have strength" \
  " zero (as if they did not exist)."

typedef nmsim_class_synapse_t *nmsim_class_synapse_ref_t;

nmsim_class_synapse_t *nmsim_class_synapse_new(double W_avg, double W_dev);
  /* Allocates a new synapse class record with the given fields. */

void nmsim_class_synapse_free(nmsim_class_synapse_t *sclass);
  /* Releases the storage of {*sclass}. */

/* INPUT/OUTPUT */

void nmsim_class_synapse_write(FILE *wr, nmsim_class_synapse_t *sclass);
  /* Writes the synapse class parameter record {p} to file {wr}, as described by
    {nmsim_class_synapse_read_INFO} below. */

nmsim_class_synapse_t* nmsim_class_synapse_read(FILE *rd);
  /* Reads the parameters of a synapse class
    from file {rd}, in the format described 
    by {nmsim_class_synapse_read_INFO} below. */
    
#define nmsim_class_synapse_read_INFO \
  "A synapse class description begins with a line \n" \
  "    begin " nmsim_class_synapse_FILE_TYPE " (format of " nmsim_class_synapse_VERSION ")\n" \
  "followed by lines in the format \"{NAME} = {VALUE}\",\n" \
  " where {NAME} is \"W_avg\" and \"W_dev\", in that order, and" \
  " a closing line \"end " nmsim_class_synapse_FILE_TYPE "\"." \
  "\n" \
  "  There must be at least one space before and after the \"=\" in the" \
  " parameter lines. The {VALUE} is a" \
  " decimal fraction.  The values \"W_avg\" and \"W_dev\" are" \
  " in millivolts.\n" \
  "\n" \
  "  " nmsim_class_synapse_INFO
    
#define nmsim_class_synapse_FILE_TYPE "nmsim_class_synapse"
    
#define nmsim_class_synapse_VERSION "2019-01-10"

/* TESTING AND DEBUGGING */

void nmsim_class_synapse_show(FILE *wr, char *pref, nmsim_class_synapse_t *sclass, char *suff);
  /* Writes the parameters {sclass} as a compact line, with prefix {pref} 
    and suffix {suff}. */
 
nmsim_class_synapse_t* nmsim_class_synapse_throw(void);
  /* Generates a random synapse class record, for testing purposes. */

void nmsim_class_synapse_compare(nmsim_class_synapse_t *sclass_read, nmsim_class_synapse_t *sclass_orig);
  /* Compares the synapse class parameters {sclass_read} read from a file
    with the expected class parameters {sclass_orig}.  Aborts if any 
    parameter does not match (within roundoff tolerance). */

vec_typedef(nmsim_class_synapse_ref_vec_t,nmsim_class_synapse_ref_vec,nmsim_class_synapse_ref_t);

#endif
