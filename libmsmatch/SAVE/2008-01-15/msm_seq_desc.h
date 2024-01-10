#ifndef msm_seq_desc_H
#define msm_seq_desc_H

/* Abstract sequence descriptors. */
/* Last edited on 2008-01-12 12:11:16 by stolfi */

#define msm_seq_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>

#include <vec.h>

#include <msm_basic.h>

typedef int msm_seq_id_t; 
  /* A {msm_seq_id_t} is an internal identifier for a sequence, 
     presently an index into a table of sequences. */

#define msm_seq_id_none -1
  /* A {msm_seq_id_t} value that means "no sequence". */

typedef struct msm_seq_desc_t 
  { msm_seq_id_t id;  /* Internal identifier of the sequence. */
    char *name;       /* External identifier of the sequence. */
    int level;        /* Filtering level. */
    int nbas;         /* Number of fundamental indices in the sequence. */
    bool_t circ;      /* Tells whether the sequence is circular. */
  } msm_seq_desc_t;
  /* Essential parameters of a sequence to be matched. In this
    library, the fields {id}, {name}, and {level} are used only for
    I/O and consistency checks.  
    
    Typically, the {id} and/or the {name} fields identify the original
    sequence, and {level} indicates the number of coarsening
    (filtering and downsampling) steps applied to the original
    sequences (so that {level} is zero for the original sequence).
    
    The sequence is assumed to have {nbas} distinct elements. If
    {circ} is FALSE, the sequence is assumed to be open; the only
    valid indices are {0..nbas-1}. If {circ} is TRUE, the sequence is
    circular; the valid indices are all integers, taken modulo
    {nbas}. */

msm_seq_desc_t msm_seq_desc_make
  ( msm_seq_id_t id,  /* Internal identifier of the sequence. */
    char *name,       /* External identifier of the sequence. */
    int level,        /* Filtering level. */
    int nbas,         /* Number of fundamental indices in the sequence. */
    bool_t circ       /* Tells whether the sequence is circular. */
  );
  /* Assembles an {msm_seq_desc_t} from the given fields. */

bool_t msm_seq_desc_same_seq(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die);
  /* Returns TRUE if the equence descriptors {sa} and {sb} refer to
    the same sequence (that is, if they have the same {name}, {id},
    and {level} fields). Otherwise, if {die} is true, fails with an
    error message; if {die} is FALSE, returns FALSE silently. The
    {name} strings are compared with {strcmp}. */

bool_t msm_seq_desc_equal(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die);
  /* Returns TRUE if sequence descriptors {sa} and {sb} are identical
    (refer to the same sequence and have the same {nbas}).
    Otherwise, if {die} is true, fails with an error message; if {die}
    is FALSE, returns FALSE silently. The {name} strings are compared
    with {strcmp}. */

/* INDEX MAPPING BETWEEN SUCCESSIVE SCALES OF RESOLUTION 

  The procedures in this section provide index mapping between a
  sequence {sFine} with {nFine} distinct elements and a filtered and
  downsampled version {sCoarse} of it, with {nCoarse} elements.
  
  These procedures assume that the filter performs a discrete
  convolution with a symmetric weight table of size {nwtb}, and
  downsamples the result with a stride of 2.
  
  Thus, if the sequence is circular, then {nCoarse == (nFine+1)/2}; if
  the sequence is open, then {nCoarse == (nFine-(nwtb-1)+1)/2}.
  
  In any case, the procedures assume that element {ix} of {sCoarse} is
  a weighted average of elements of {sFine}, centered on element {2*ix
  + (nwtb-1)/2} of the latter.
  
  The mappings between the two index ranges are monotonic
  non-decreasing: i.e. increasing one index does not decrease the
  other. Moreover, each index {ix} of {sCoarse} maps to two distinct
  elements in {sFine}. However, if {sFine} is circular and {nFine} is
  odd, then element {nCoarse-1} corresponds to {nFine-1} alone. In
  this case, there is a `hiccup' in the mapping at every index that is
  a multiple */

int msm_seq_desc_filtered_size(int nFine, bool_t circ, int nwtb);
  /* Computes the number of distinct elements in the filtered and downsampled 
    version of a sequence with {nFine} elements and circularity {circ}. */

int msm_seq_desc_map_index_to_coarser(int ixFine, int nFine, int nCoarse, int nwtb);
  /* Maps the index {ixFine} of an element in sequence {sFine} to the
    index {ixCoarse} of the approximately corresponding element in
    the filtered and downsampled sequence {sCoarse}. */

int msm_seq_desc_map_index_to_finer(int ixCoarse, int nCoarse, int nFine, int nwtb);
  /* Maps the index {ixCoarse} of an element in the filtered
    sequence {sCoarse} to the index {ixFine} of the approximately
    corresponding element in the original sequence {sFine}.
    
    The procedure is monotonic non-decreasing: i.e. increasing
    {ixCoarse} does not decrease the result.  Usually the procedure
    maps two distinct If {nFine} is odd,
    */

/* RANDOM SEGMENTS */

void msm_seq_desc_throw_segment(int ns, bool_t circ, int len, int *ini, int *fin);
  /* Selects initial and final sample indices {*ini,*fin} for a random 
    segment with {len} samples, taken from a sequence with {ns} samples
    and circularity {circ}. */ 

/* PRINTOUT */

void msm_seq_desc_write
  ( FILE *wr, 
    char *pre, 
    msm_seq_desc_t *s, 
    int idsize, 
    int namesize, 
    int nssize,
    char *suf
  );
  /* Writes to {wr} the sequence descriptor {s}, in a format
    compatible with {msm_seq_desc_read}; namely,
    
     "{pre} {s.id} {s.name} {s.level} {s.nbas} {s.circ} {suf}". 
      
     The delimiters {pre} and {suf} will be omitted if NULL. The fields
    {s.id}, {s.name}, and {s.nbas} will be blank-padded to {idSize},
    {nameSize}, and {nbasSize} characters, respectively. The field
    {s.circ} will be printed as "T" or "F". */

msm_seq_desc_t msm_seq_desc_read(FILE *rd, char *pre, char *suf);
  /* Parses from {rd} a sequence descriptor, in the format 
  
      "{pre} {id} {name} {level} {nbas} {circ} {suf}". 
      
    If the arguments {pre} and {suf} are not NULL, they must be
    present in the input. In general, consecutive fields should be
    separated by one or more spaces. The {name} field and the
    {pre,suf} arguments must not contain embedded spaces or other
    formatting chars. The {circ} field should be "T","t","1" for TRUE,
    "F","f","0" for FALSE. */

#endif
