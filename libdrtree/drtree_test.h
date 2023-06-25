/* Tools for testing {drtree} fundions. */
#ifndef drtree_test_H
#define drtree_test_H
/* Last edited on 2023-06-16 18:33:57 by stolfi */

#define drtree_test_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 
#include <argparser.h> 

#include <drtree.h> 

typedef struct drtree_test_options_t
  { int32_t nIndivs;  /* Number of individuals in population. */
    int32_t nRoots;   /* Number of root individuals. */
    int32_t tStart;   /* Nominal start time of evolution. */
    int32_t tStop;    /* Nominal final time of evolution. */
    bool_t orphans;   /* If false, no roots allowed after {tStart}. */
    int32_t ageMax;   /* Max individual age. */
    int32_t nchMax;   /* Max num of children. */
    int32_t tRef;     /* A reference time for lineage highlighting. */
    char *outPrefix;  /* Prefix for output files. */
  } drtree_test_options_t;
  /* The command line options.  See {drtree_test_options_HELP} 
    and {drtree_test_options_INFO} below.  */

#define drtree_test_num_indivs_MAX 1000
#define drtree_test_num_times_MAX 300
#define drtree_test_num_children_MAX 10
  /* Limits for command line options.  */

drtree_test_options_t *drtree_test_parse_options(argparser_t *pp);
  /* Parses the command line arguments. */
 
#define drtree_test_options_HELP \
  "    -nIndivs {nIndivs} \\\n" \
  "    -nRoots {nRoots} \\\n" \
  "    -tStart {tStart} \\\n" \
  "    -tStop {tStop} \\\n" \
  "    -orphans { \"T\" | \"F\" } \\\n" \
  "    -ageMax {ageMax} \\\n" \
  "    -nchMax {nchMax} \\\n" \
  "    -tRef {tRef} \\\n" \
  "    -outPrefix {outPrefix}"
  
#define drtree_test_options_INFO \
  "  -nIndivs {nIndivs}\n" \
  "    This mandatory argument specifies the number of individuals to generate in for test.\n" \
  "\n" \
  "  -nRoots {nRoots}\n" \
  "    This mandatory argument specifies the desired number of root individuals.  Must be at most {nIndivs}.\n" \
  "\n" \
  "  -tStart {tStart}\n" \
  "  -tStop {tStop}\n" \
  "    These mandatory arguments specify the initial and final time of the simulation. The values may be negative.\n" \
  "\n" \
  "  -orphans { \"T\" | \"F\" }\n" \
  "    This mandatory argument specifies whether root (parent-less) individuals may be born after {tStart}.\n" \
  "\n" \
  "  -ageMax {ageMax}\n" \
  "    This mandatory argument specifies the maximum age of an individual..\n" \
  "\n" \
  "  -nchMax {nchMax}\n" \
  "    This mandatory argument specifies the maximum number of children per individual.\n" \
  "\n" \
  "  -tRef {tRef}\n" \
  "    This mandatory argument specifies a reference time in the range {tStart..tStop}.  The" \
  " lineages that are surviving at this time will be highlighted somehow.\n" \
  "\n" \
  "  -outPrefix {outPrefix}\n" \
  "    This mandatory argument specifies a common prefix for all output file names."

void drtree_test_create_individuals
  ( int32_t ni,       /* Number of individuals in population. */
    int32_t nRoots,   /* Number of root individuals. */
    int32_t tStart,   /* Nominal start time of evolution. */
    int32_t tStop,    /* Nominal final time of evolution. */
    bool_t orphans,   /* If false, no roots allowed after {tStart}. */
    int32_t ageMax,   /* Max individual age. */
    int32_t nchMax,   /* Max num of children. */
    drtree_node_t dt[]
  );
  /* Fills {dt[0..ni-1]} with the description of {ni} individuals
    with random life spans and parenting, as specified 
    by the parameters {nRoots,tStart,tStop,orphans,ageMax,nchMax}.
    See [drtree_test_options_INFO} for their meanings.
    
    The life spans of the individuals will usually cover the whole interval
    {o->tStart..o->tStop} but will not be contained in it. */

#endif

