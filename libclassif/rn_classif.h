/* rn_classif.h --- tools for vector classifiers. */
/* Last edited on 2024-12-21 11:22:38 by stolfi */

#ifndef rn_classif_H
#define rn_classif_H

#include <stdio.h>
#include <math.h>

#include <r2.h>
#include <frgb.h>
#include <jsrandom.h>
  
/* VECTOR CLASSIFICATION PROBLEMS */

/* A (discrete) /vector classification function/ is a function that
  assigns to vectors from some space {\Reals^NA} a /class index/ {cl}
  in some range {1..NC}.

  The components of the vector are called /attributes/ and the subset
  of {\Reals^NA} that maps to a given class {cl} is called
  the /domain/ of that class.  

  The function need not be defined on the whole {\Reals^NA}. Here the
  class index 0 is used to indicate `undefined class' or `invalid
  vector' or `do not care'. The region that maps to class 0 is not
  considered a class domain but the /background/, complement of all 
  class domains.
  
  A (discrete) /vector classification problem/ is a vector
  classification function with a probability distribution over its
  valid classes. */
  
typedef int rn_classif_labeler_t(int NA, int NC, double p[]);
  /* Type of a vector classification function from {\Reals^NA} 
    with {NC} valid classes.
    
    A function of this type should return the index of a class in
    {1..NC} whose domain contains the given vector {p[0..NA-1]}. If
    {p} is outside all domains it should return 0. If {p} is close to
    a domain boundary the result is undefined. (However it should
    preferably err in favor of domains with smaller dimension, so as
    to make them thick enough for visualization.) 
    
    Each combination of {NA} and {NC} defines a different
    classification function, even with the same {rn_classif_labeler_t}
    procedure. The procedure may fail for certain combinations of {NA}
    and {NC}. */
    
typedef void rn_classif_thrower_t(int i, int NA, int NC, double p[], int *classP);    
  /* Type of a sampler procedure associated to some vector classification
    function from {\Reals^NA} with {NC} valid classes. 
    
    A procedure of this type should select a valid class domain {cl} and generate a random
    vector in {\Reals^NA} in that domain.  It should store the 
    attributes into {p[0..NA-1]} and the class {cl} into {*classP}.
    
    The attributes of a generated sample depend on the integer {i} and on the
    current state of the number generator {drandom()}, which must have
    been properly seeded by the client. Usually each attribute is in
    the range {U=[-1 _ +1]}. The procedure may choose the class domain
    to be sampled as a periodic function of {i}, so that if {i} ranges
    over {1..N} the samples will be more unformly distributed among the
    classes. 
    
    Each combination of {NA} and {NC} defines a different sampling
    function, even with the same {rn_classif_thrower_t} procedure.
    The procedure may fail for certain combinations of {NA} and
    {NC}. */
  
typedef struct rn_classif_problem_t 
  { int NA;                    /* Total number of attributes per sample. */
    int NC;                    /* Number of classes. Classes are numbered {1..NC}. */
    rn_classif_labeler_t *lab; /* Labelling procedure */
    rn_classif_thrower_t *gen; /* Generation procedure */
  } rn_classif_problem_t;
  /* A Discrete classification problem. */

void rn_classif_class_count(int NS, int class[], int NC, int num[]); 
  /* Stores into {num[cl]} the number of elements of {class[0..NS-1]}
    with value {cl}, where {cl} ranges in {0..NC}. Note that {num} must
    have {NC+1} elements. */

void rn_classif_class_count_print(FILE *wr, int NC, int num[]);
  /* Prints the class count table {num[0..NC]}. */

void rn_classif_sample_print(FILE *wr, char *pre, int NA, double p[], char *sep, char *suf);
  /* Writes to {wr} the attributes {p[0..NA-1]}, preceded by {pre},
    separated by {sep}, and terminated by {suf}. */
 
typedef struct rn_classif_dataset_t 
  { int NS;        /* Number of samples. Samples are numbered {0..NS-1}. */
    int NA;        /* Total number of attributes per sample. */
    double **smp;  /* {smp[i][t]} is attribute {t} of sample {i}. */
  } rn_classif_dataset_t;
  /* A dataset for a vector classifier, cosisting of {NS} sample
    vectors, each with {NA} attributes. */
    
rn_classif_dataset_t *rn_classif_dataset_new(int NS, int NA);
  /* Allocates the storage for a dataset {D} with {NS} samples and {NA} attributes per sample.
     Allocates the vector {D.smp} but not the individual attribute vectors. */

void rn_classif_dataset_free(rn_classif_dataset_t *D);
  /* Deallocates all storage of {D} including the header {*D} itself. */

void rn_classif_dataset_label(rn_classif_dataset_t *D, int NAL, int NCL, rn_classif_labeler_t *lab, int class[]);
  /* Classifies each sample in {D} with the labeler {lab} into classes {0..NCL},
    and stores the result in {class[0..D.NS-1]}. The classifier is given the first {NAL}
    attributes of each sample, so one must have {NAL <= D.NA}. */
 
void rn_classif_dataset_write(FILE *wr, rn_classif_dataset_t *D, int cl, int class[]);
  /* Writes {D} to {WR} in a human- and machine-readable format. If
    {class} is not NULL, then it must have {D.NS} elements, and the
    procedure writes only the samples {i} such that {class[i] =
    cl}. */
  
void rn_classif_compare(int NS, int cold[], int cnew[], int *nsP, int *nfP);
  /* The vectors {cold} and {cnew} should have {NS} elements each. The procedure
    computes the number {ns} of `successes' (indices {i} in {0..NS-1}
    such that {cold[i] == cnew[i]} and the number {nf} of `failures'
    (all the other cases). These counts are returned in {*nsP} and {*nfP}. */
  
void rn_classif_cross_matrix_build(int NS, int cold[], int NCold, int cnew[], int NCnew, int **nccP);
  /* The vectors {cold} and {cnew} should have {NS} elements each, in
    the ranges {0..NCold} and {0..NCnew}, respectively. The procedure builds a
    /cross-classification matrix/ {ncc} where each entry
    {ncc[old,new]} is the number of positions {i} with {cold[i] = old}
    and {cnew[i] = new}. The `non-class' 0 is counted like any valid
    class, so the array has {NCnew+1} columns and {NCold+1} rows.
    
    The array is linearized by rows, so {ncc[old,new]} is actually in
    {ncc[old*(NCnew+1)+new]}. The array {ncc} is allocated by the
    procedure and its address is returned in {*nccP}. */

void rn_classif_cross_matrix_print(FILE *wr, int NCold, int NCnew, int *ncc, bool_t score);
  /* Prints to {wr} the cross-classification matrix {ncc} which is
    assumed to have {NCold} rows (`old' classes) and {NCnew} columns
    (`new' classes).  If {score} is true prints also the number of elements in the diagonal
    (`successes') and off the diagonal (`failures'). */

typedef double rn_classif_pq_dist_t(int NA, double p[], double q[]);
typedef double rn_classif_pi_dist_t(int NA, double p[], int i);
typedef double rn_classif_ij_dist_t(int i, int j);
  /* Types of distance functions for samples and elements of sample sets.
  
    A function of this type should returns some measure of the distance(not
    necessarily a metric in the mathematical sense) between the two given
    samples.  The samples may be given as two vectors {p[0..NA-1],q[0..NA-1]}
    or as a vector {p[0..NA-1]} and an index {i} into some sample list,
    or two indices {i,j} into two such list (or the same list). */

/* ITERATIVE IMPROVEMENT ("TRAINING") */

void rn_classif_find_nearest_neighbor_index(int NS, int NA, double p[], rn_classif_pi_dist_t *dist, int *iminP, double *dminP);
  /* Finds in dataset of {NS} samples the sample whose first {NA}
    attributes are closest to {p[0..NA-1]} in the metric {dist}.
    Returns the index of said element in {*iminP} and the distance in
    {*distP}. */

int rn_classif_find_nearest_in_dataset(rn_classif_dataset_t *M, double HM[], int NA, double p[], rn_classif_pq_dist_t *dist);
  /* Finds the index {i} such that the first {NA} attributes of sample
    {M.smp[i]} are closest to {p[0..NA-1]}, in the metric {dist}.
    Requires {NA<=M.NA}. If {HM} is not null it must be a vector of
    {M.NS} handicaps; the distance {dist} is then replced by
    {DIST(p,M.smp[i]) = max(dist(p,M.smp[i]),HM[i])} */

void rn_classif_nn_label_dataset
  ( rn_classif_dataset_t *M,    /* Model samples. */
    double HM[],                /* {HM[i]} is the handicap of sample {M.smp[i]}. */
    int classM[],               /* {classM[i]} is the class of sample {M.smp[i]}. */
    rn_classif_dataset_t *C,    /* Samples to be classified. */
    rn_classif_pq_dist_t *dist, /* Sample distance function. */
    int classC[]                /* (OUT) {classC{j]} is the class assigned to {C.smp[j]}. */
  );
  /* Calls {i = rn_classif_find_nearest_in_dataset(M,HM,M.NA,C.smp[j],dist)}
    for each sample {C.smp[j]} and sets {classC[j]} to {classM[i]}.
    The vectors {HM} and {classM} must have {M.NS} elements.
    The vector {classC} must have {C.NS} elements. */

typedef void rn_classif_handicap_proc_t(rn_classif_dataset_t *M, int classM[], double HM[]);
  /* Type of a procedure that computes handicaps {HM[0..M.NS-1]} for a 1NN
    classifier {(M,HM,classM)}, given {M} and {classM}. */

void rn_classif_nn_improve
  ( rn_classif_dataset_t *M,    /* (IN/OUT) Model samples. */
    double HM[],                /* (IN/OUT) {HM[i]} is the handicap of sample {M.smp[i]}. */
    int classM[],               /* (IN/OUT) {classM[i]} is the given class of sample {M.smp[i]}. */
    rn_classif_dataset_t *R,    /* (IN/OUT) Training samples. */
    int classR[],               /* (IN/OUT) {classR[i]} is the given class of sample {R.smp[i]}. */
    int classX[],               /* (OUT) {classX[i]} is the class of {R.smp[i]} assigned by {M}. */
    rn_classif_pq_dist_t *dist, /* Sample distance function. */
    rn_classif_handicap_proc_t *hproc, /* Procedure that computes handicaps for {M}. */
    int maxIters,               /* Maximum iterations. */
    bool_t sameClass,           /* TRUE only swaps samples of the same class. */
    bool_t verbose              /* TRUE to show the progress at each iteration. */
  );
  /* Try to improve the NN model {M} by replacing some of its elements
    by elements of {R}. 
    
    The goal is to make the labeling of {R.smp} by the 1NN classifier
    {(M,HM,classM)} be as close as possible to the given classification
    {classR}. Each tentative improvement consists in swapping some
    sample of {M} by some sample of {R} (usually one that was
    misclassified by {M}) together with their given classifications in
    {classM} and {classR}. If the swap does not succeed it is undone,
    so the final result is never worse than the original one. This
    basic step is repeated at most {maxIter} times.
    
    If {sameClass} is TRUE then the swap is restricted to samples of
    the same class, unless {M} has no samples with the same class as
    the chosen {R} sample.
    
    The procedure {hproc} is called to compute the handicap vector {HM} 
    upon entry and every time the model {M} is changed. 
    
    On output, the sets {M,R} and their classifications
    {classM,classR} may have been partly swapped and rearranged. The
    vector {HM} will contain the handicaps computed for the final {M}
    by {hproc}. The vector {classX} will contain the labeling of {R}
    implied by {(M,HM,classM)}. */
 
void rn_classif_nn_print_handicaps(FILE *wr, char *title, int NS, double H[]);
  /* Writes the sample handicap vector {H[0..NS-1]} to the file {wr},
    with the given {title}. */
    
frgb_t *rn_classif_pick_class_colors(int NC);
  /* Returns an array of {NC+1} colors for classes
    {0..NC}.  Color 0 is always white, the rest
    has increasing brightness and variable hues. */

#endif
