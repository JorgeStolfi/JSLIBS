#define PROG_NAME "test_nnclassif"
#define PROG_DESC "test nearest neihbor classifiers"
#define PROG_VERS "1.0"

#define make_test_classif_data_C_COPYRIGHT \
  "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2020-10-11 03:22:30 by jstolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -problem {PROBLEM} \\\n" \
  "    -seed {SEED} \\\n" \
  "    -samples {NS_MODEL} {NS_REFINE} {NS_EVAL} \\\n" \
  "    -prefix {PREFIX} \\\n" \
  "    [ -attributes {NAR} {NAI} ] \\\n" \
  "    [ -classes {NC} ] \\\n" \
  "    [ -noise {SIGMA} ] \\\n" \
  "    [ -opf { N | Y } ] \\\n" \
  "    [ -refine {REFINE_ITERS} ] \\\n" \
  "    [ -swap { same | full } ] \\\n" \
  "    [ -grid { N | Y } ] \\\n" \
  "    [ -verify { N | Y } ] \\\n" \
  "    [ -image {IMG_SIZE} {SUBSMP} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESCRIPTION "\n" \
  "\n" \
  "OUTPUTS\n" \
  PROG_INFO_OUTPUTS "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTIONS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  The San Francisco Exploratorium.\n" \
  "\n" \
  "SMELL ALSO\n" \
  "  The roses along the road.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2010-05-24 by J. Stolfi: created program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " make_test_classif_data_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESCRIPTION \
  "  Tests the nearest-neighbor (1NN) vector classifier, with a model of {NS_MODEL} samples" \
  " for a /discrete classification function/ specified by {PROBLEM}, {NAR} and {NC}.\n" \
  "\n" \
  "  Each sample is a vector of attributes.  Each attribute is a real number.  The total number of" \
  " attributes per sample is {NA = NAR + NAI} where {NAR} and {NAI} are" \
  " user-specifiable numbers.  The first {NAR} attributes are \"relevant\" to the test problem;" \
  " the the remaining {NAI} attributes are \"irrelevant\" random numbers, uniformly distributed" \
  " in {U = [-1 _ +1]}, appended to them.\n" \
  "\n" \
  "  Considering only the {NAR} relevant attributes, the classification function" \
  " defined by {PROBLEM,NAR,NC} define a space of /valid" \
  " attribute vectors/ contained in {U^NAR} that is partitioned {NC}" \
  " disjoint /class domains/ numbered {1..NC}.  The number" \
  " and shape of the domains, and the random sampling probability distribution" \
  " in each domain depends on the.  For some {PROBLEM}s" \
  " the number of attributes {NAR} and/or classes {NC} may be restricted or fixed.  The class domains" \
  " may not cover the whole of {U^NAR} and may or may not touch each other. Some" \
  " domains may be zero-measure sets with dimension smaller than {NAR}.\n" \
  "\n" \
  "  The program" \
  " generates {NS} samples in these domains (where {NS} is either {NS_MODEL}, {NS_REFINE}," \
  " or {NS_EVAL}), and assigns to each sample" \
  " a /nominal class/ which is the index of the containing domain.  The samples" \
  " may be generated at random or may be taken from a regular grid.\n" \
  "\n" \
  "  When sampling on a grid, the number {NS} is implicitly rounded up to the next" \
  " perfect power {NG^NAR}, and the problem is sampled at a regular grid of points" \
  " with {NAR} samples along each axis and spanning the hypercube {U^NAR}.  Only" \
  " those samples that fall inside one of" \
  " the class domains will be written out.  In that case the actual number of" \
  " samples in the output may be substantially larger or smaller than {NS}.  Note" \
  " that the arrangement of samples is regular only in {U^NAR}, not in {U^NA}\n" \
  "\n" \
  "  After generating each sample (whether randomly or from a grid) and" \
  " appending the {NAI} irrelevant attributes, the program" \
  " adds to each attribute (relevant or irrelevant)" \
  " a Gaussian noise with mean 0 and deviation {SIGMA}.  The noise is truncated to {4*SIGMA}" \
  " so that each final attribute is in the range {V = [-VMAX _ +VMAX]}" \
  " where {VMAX = 1 + 4*SIGMA}.  The nominal class of the sample is retained," \
  " even though the perturbation may cause the" \
  " attribute vector to fall outside the domain of that class and" \
  " possibly into the domain of a different class.\n" \
  "\n" \
  "  The program generates three datasets, a /model set/ {M} with {NS_MODEL} samples," \
  " a /refinement set/ {R} with {NS_REFINE} samples, and an /evaluation set/ {E} with" \
  " {NS_EVAL} samples.  The model set {M}, always random, is used with the nearest" \
  " neighbor algorithm to classify other samples.  The evaluation set {E}, either random or grid-like, is used" \
  " to check the accuracy of the 1NN classifier with the model.  The refinement set {R}, also" \
  " random, is used as a source of alternative samples that can be iteratively swapped with those of {M} to" \
  " improve the performance of the 1NN classifier.\n" \
  "\n" \
  "  The 1NN algorithm has three parameters: the list of representative samples {M}, a /class" \
  " table/ {classM}, and a /handicap table/.  For each {s} in {M}," \
  " {classM(s)} is a positive class index" \
  " and {HM(s)} is a non-negaive real number.  The algorithm itself defines another classification function" \
  " for {Reals^NA}, namely it assigns to each sample vector {t} in {Reals^NA} the class" \
  " {classM(s)} of the sample {s} in {M} that minimizes {max(H(s),dist(t,s))}, where {dist}" \
  " is some metric (currently the Euclidian one).  The handicaps {H(m)} are either zero," \
  " or are computed by the program from {M} and {classM} by the OPF algorithm" \
  " (Spina, Falcão and Suzuki, 2009)."
  
#define PROG_INFO_OUTPUTS \
  "  If \"-image\" is specified, the program will write a set of PPM images" \
  " files called \"{PREFIX}-{TAG}-C.ppm\" for various {TAG} ad {KIND} strings.  showing the\n" \
  " classes of the ideal classifier defined by {PROBLEM,NAR,NC} (with {TAG} = \"I-tru\"), and of the 1NN" \
  " classifiers {(D,HD,classD)} where {D} is either {M} ({TAG} = \"M-ini\") or {R} ({TAG} = \"R-ini\")." \
  " In any clase {classD} is the labeling of {D} samples defined by the ideal" \
  " classifier {PROBLEM,NAR,NC}, and {HD} is a set of handicaps" \
  " computed from {D} and {classD} according to the \"-opf\"" \
  " parameter.  If {NS_REFINE} is positive it also generates images of the 1NN classifier with" \
  " the model set {M'} after refinement ({TAG} = \"M-ref\").  In these images, each class domain will" \
  " be painted with a distinct nonwhite color, with brightness" \
  " increasing from 1 to {NC}.  The background (class 0) is painted white.\n" \
  "\n" \
  "   For each of the above 1NN classifiers, the program will also produce an" \
  " image file called \"{PREFIX}-{TAG}-D.ppm\"" \
  " showing the discrepancy between its classification and that of the" \
  " deal classifier {PROBLEM,NAR,NC}.  Points where the two classifiers agree, or where" \
  " the ideal classifier returns 0, are shown in white; points where they" \
  " disagree are shown in red.\n" \
  "\n" \
  "   The number of relevant attributes {NAR} must be 2, and the X and Y axis are" \
  " attribute 0 and attribute 1, respectively. Each image will span the augmented" \
  " square {V^2}.  In all these images, if \"-noise\" is requested, each sample is perturbed by" \
  " the noise, and then the position of the original sample is painted with the color" \
  " of the domain that contains the perturbed point.  If {NAI} is positive, the" \
  " irrelevant attributes will set to random values in the set {D}, as explained" \
  " above, and to zeros in the samples being imaged.  See the \"-image\" option for further details."

#define PROG_INFO_OPTIONS \
  "  -problem {PROBLEM}\n" \
  "    This mandatory argument specifies the ideal class domains and" \
  " the sampling probability distributions.\n" \
  "\n" \
  "  -samples {NS_MODEL} {NS_REFINE} {NS_EVAL}\n" \
  "    This mandatory argument specifies the number of samples to generate" \
  " in each of the three sets.  If {NS_REFINE} is zero then no refinement" \
  " will be tried.\n" \
  "\n" \
  "  -seed {SEED}\n" \
  "    This mandatory argument specifies a seed for the random number generator.\n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    This mandatory argument specifies the prefix for all output filenames.\n" \
  "\n" \
  "  -attributes {NAR} {NAI}\n" \
  "    This optional argument specifies the number of attributes in" \
  " each sample, namely {NAR} relevant ones and {NAI} irrelevant ones.  The allowed" \
  "  values of {NAR} depend on the {PROBLEM}.  If {NAR} is zero" \
  " (the default) the number of relevant attributes depends on the {PROBLEM}.  If {NAI} is zero" \
  " (the default) there are no irrelevant attributes.\n" \
  "\n" \
  "  -classes {NC}\n" \
  "    This optional argument specifies the number of classes into which" \
  " the space of all attribute vectors is partitioned.  If {NC} is zero" \
  " (the default) the number of classes is selected by the {PROBLEM}.  The" \
  " allowed range of {NC} depends on the {PROBLEM}.\n" \
  "\n" \
  "  -noise {SIGMA}\n" \
  "    This optional argument specifies the standard deviation" \
  " of the Gaussian perturbation to be added to each attribute" \
  " after sampling.  Note that the noise may move some samples" \
  " into the domain of a different class.  The default is {SIGMA=0}" \
  " (no perturbation).\n" \
  "\n" \
  "  -opf { N | Y }\n" \
  "    This optional argument requests the use of handicaps" \
  " when computing distances to the representantive samples" \
  " in the classifer.  The handicaps are" \
  " computed by the OPF algorithm   The default is \"N\" (no OPF handicaps).\n" \
  "\n" \
  "  -refine {REFINE_ITERS}\n" \
  "    This optional argument directs the program to attempt an iterative" \
  " improevment of the model set {M} by swapping some of its elements with" \
  " those of the refinement set {R} so as to improve the accuracy of {M} when classifying {R}.  At most" \
  " {REFINE_ITERS} iterations will be performed.  The default" \
  " {REFINE_ITERS} depends on the size of {M} and {R}, with the \"full\" option.\n" \
  "\n" \
  "  -swap { same | full }\n" \
  "    This optional argument specifies whether the refinement iteration should swap only" \
  " elements of the same class if possible (\"same\") or may swap " \
  " elements of distinct classes may be swapped as well (\"full\").  The default is \"full\".\n" \
  "\n" \
  "  -verify { N | Y }\n" \
  "    This optional argument specifies whether the procedure should" \
  " check that each randomly generated sample belongs to the domain of its nominal" \
  " class.  Errors detected before the" \
  " random perturbation are fatal; errors after the perturbation are" \
  " counted and reported to {stderr}.  The default is \"N\" (no checking).\n" \
  "\n" \
  "  -grid { N | Y }\n" \
  "    This optional argument specifies whether the procedure should" \
  " generate the evaluation set {E} at a regular grid instead of" \
  " randomly.  The default is \"Y\" (grid sampling).\n" \
  "\n" \
  "  -image {SIZE} {SUBSMP}\n" \
  "    This optional argument directs the program to write PPM" \
  " image files called \"{PREFIX}-{TAG}-{KIND}.ppm\" showing the domain" \
  " classes.  This option can be used" \
  " only when {NAR=2}. The images will have {SIZE} columns and {SIZE}" \
  " rows.  The program will generate {SUBSMP^2} samples" \
  " inside each pixel, and average their colors to obtain the pixel" \
  " color."

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <float_image.h>
#include <r2.h>
#include <rn.h>
#include <rn_classif.h>
#include <rn_classif_opf.h>
#include <rn_classif_uint16_image.h>
#include <rn_classif_test.h>
  
#define MAX_ATTRIBS 128
  /* Max num of relevant attributes (param safety only) */

#define MAX_CLASSES 128
  /* Max num of classes (param safety only) */

#define MAX_SAMPLES 100000000
  /* Max num of samples to generate (param safety only) */

#define MAX_SEED (~0LU)
  /* Max seed value (param safety only) */

#define MAX_NOISE 10.0
  /* Max variance of post-sampling noise (param safety only) */

#define MAX_IMAGE_SIZE 1024
  /* Max rows and columns in image (param safety only) */

#define MAX_SUBSMP 5
  /* Max pixel subsampling order for image (param safety only) */

#define MAX_REFINE_ITERS 100000
  /* Max refinement iterations (param safety only) */

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { char *problem;      /* Name of classification problem. */
    int attributes_rel; /* Number of relevant attributes desired. */
    int attributes_irr; /* Number of irrelevant attributes to append. */
    int classes;        /* Number of classes desired. */
    int samples_model;  /* Number of samples in NN model. */
    int samples_refine; /* Number of samples in refinement (retraining) set. */
    int samples_eval;   /* Number of samples in evaluation set. */
    uint32_t seed;      /* Seed for randomizer. */
    double noise;       /* Deviation of post-sampling noise. */
    bool_t opf;         /* TRUE to use OPF handicaps in 1NN classifier. */
    bool_t grid;        /* TRUE to use grid sampling to obtain the {E} set. */
    bool_t verify;      /* TRUE to verify the generated samples. */
    int refine_iters;   /* Number of refinement iterations. */
    bool_t swap_same;   /* TRUE if refinement should swap only in same class. */
    char *prefix;       /* The ouput file prefix. */
    /* Image output options: */
    int image_size;     /* Height and width of PPM image, or 0 if none. */
    int image_subsmp;   /* Will generate {subsmp*subsmp} samples in each pixel. */
  } options_t;

/* CANNED PROBLEMS */

typedef enum
  { problem_kind_SATURN, /* The "saturn" ("B9") dataset from Papa et al.(2009); {NC = 2}. */
    problem_kind_PETALS, /* The "petals" ("B10") dataset from Papa et al.(2009); {NC = 4}. */
    problem_kind_VESSEL, /* The "boat" ("B11") dataset from Papa et al.(2009); {NC = 3}. */
    problem_kind_MBALLS, /* Each class {2..NC} comprises {(NC-1)^(NAR-1)} balls; class 1 is background. */
    problem_kind_SHELLS, /* Classes {2..NC} are concentric shells; class 1 is background and center. */
    problem_kind_STAOLI, /* Classes {2..3} are a fat and thin box; class 1 is background. */
    problem_kind_NUMBER  /* Number of valid problem kinds. Must be last. */
 } problem_kind_t;
 /* Numeric problem kinds.  The valid kinds are {0..problem_kind_NUMBER-1} */

char *problem_kind_name[problem_kind_NUMBER] = 
  { [problem_kind_SATURN] = "saturn", 
    [problem_kind_PETALS] = "petals", 
    [problem_kind_VESSEL] = "vessel", 
    [problem_kind_MBALLS] = "mballs", 
    [problem_kind_SHELLS] = "shells", 
    [problem_kind_STAOLI] = "staoli" 
  };
  /* External names of the problem kinds. */
  
typedef void problem_def_proc_t(int *NAR, int *NC);
  /* Type of a problem definition procedure.  It sets {*NC}
    and {*NAR} to suitable defaults if any of them are zero, then 
    performs range checking. */

/* INTERNAL PROTOTYPES */

typedef rn_classif_dataset_t dataset_t;
typedef rn_classif_problem_t problem_t;

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

problem_t *get_problem(char *name, int NAR, int NC);
  /* Returns a {problem_t} given the problems {name} and the desired
    numbers of attributes {NAR} and classes {NC}.  Fails if {NAR,NC}
    are not a valid combination for that problem. */
    
void generate_dataset(problem_t *P, int NAI, int NS, bool_t grid, options_t *o, dataset_t **DP, int **classDP);
  /* Returns a dataset containing {NS} samples of the classification
    problem {P}, augmented with {NAI} irrelevant attributes, according to options {o}. */

void generate_raw_dataset_random(problem_t *P, int NA, int NS, uint32_t seed, bool_t verify, dataset_t **DP, int **classDP);
void generate_raw_dataset_grid(problem_t *P, int NA, int NS, uint32_t seed, dataset_t **DP, int **classDP);
  /* Returns a {dataset_t} containing random or grid samples of the classification
    problem {P} according to options {o}.  Only the {P->NAR} relevant attributes
    of each sample are set.  Does not add any noise. */
   
void dataset_stats_print(FILE *wr, int NS, int NC, int classD[]);
  /* Prints statistics on number and percentage of samples per 
    class in the dataset {D}. */

void output_dataset(dataset_t *D, int classD[], int NC, options_t *o);
  /* Writes the dataset as {NC} files "{o->prefix}-c{NNN}.dat" where
    {CCC} is each class in 3 digit format. See {rn_classif_dataset_write}
    for the format. */
    
void verify_dataset(problem_t *P, dataset_t *D, int classD[], char *nameD);
  /* Compares the classification {classD} for the samples of {D}
    with the classification specified by {P.lab} with {P.NA} atributes
    and {P.NC} classes.  Assumes that {classD} has only classes in {0..P.NC} too.
    If there are discrepancies, prints a cross-classification matrix with the 
    datset {nameD} in the title.*/

void compute_handicaps(dataset_t *M, int classM[], bool_t opf, bool_t verbose, double HM[]);
  /* Computes the 1NN handicaps for the dataset {M} with classes {classM}.  If {opf} is true 
    uses the OPF algorithm, else sets all handicaps to zero.  */

int evaluate_nn_classifier(int NC, dataset_t *M, int classM[], double HM[], char *nameM, dataset_t *E, int classE[], char *nameE);
  /* Uses the 1NN classifier {(M,classM,HM)} to classify each sample in the dataset {E},
    ans compares the result with the classification {classE}.  Prints a cross-classification table
    and returns the number of discrepancies {nf}, with the dataset {nameM,nameE} in the title. */

void refine_nn_model(int NC, dataset_t *M, int classM[], double HM[], dataset_t *R, int classR[], options_t *o, bool_t verbose);
  /* Tries to improve the 1NN classifier {(M,HM,classM)} by swapping some elements with those of {R}. */

void output_ideal_image(problem_t *P, char *tag, options_t *o);
  /* Writes PPM file "{o->prefix}-{tag}-C.ppm" showing the class domains of the problem {P}. */

void output_dataset_image(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o);
  /* Writes PPM file "{o->prefix}-{tag}-P.ppm" showing the points of the dataset {D}. 
  !!! Should also write "{o->prefix}-{tag}-H.ppm" showing the handicaps as circles. */

void output_nn_images(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o);
  /* Writes PPM file "{o->prefix}-{tag}-C.ppm" showing the class domains of the 1NN classifier
    {(D,classD,HD)}. Then writes another image file "{o->prefix}-{tag}-D.pgm" that shows where that classification
    differs from the classification of {P}. */

void output_nn_hand_cmp_image(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o);
  /* Writes PPM file "{o->prefix}-{tag}-X.ppm" showing the effect of the handicaps {HD}
    on the 1NN classifier; that is, comparing the 1NN classifier {(D,HD,classD)} with {(D,NULL,classD)}. */

void output_image(char *prefix, char *tag, char *kind, uint16_image_t *img);
  /* Writes the image {img} to file "{prefix}-{tag}-{kind}.ppm". */

void add_random_noise(dataset_t *D, double sigma, uint32_t seed);
  /* Add the post-sampling noise.  The perturbations are a fixed 
    random function of the seed, indepednent of the random
    values used in sampling. */

void append_irrelevant_attributes(dataset_t *D, problem_t *P, uint32_t seed);
  /* Appends {D.NA - P.NAR} random values uniform in {U} to each 
    sample in {D}.  The values are a fixed 
    random function of the seed, independent of those used in 
    sampling and noise. */

void get_grid_indices(int g, int NG, int NAR, int ix[]);
  /* The integer {g} must be in the range {0..NG^NAR-1}. Decomposes
     {g} into its {NAR} digits {ix[0..NAR-1]} in base {NG}. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    problem_t *P = get_problem(o->problem, o->attributes_rel, o->classes);
    if (o->image_size > 0) 
      { output_ideal_image(P, "I-tru", o); }
    
    /* Generate the model dataset (random): */
    dataset_t *M = NULL;
    int *classM = NULL;
    assert(o->samples_model > 0);
    generate_dataset(P, o->attributes_irr, o->samples_model, FALSE, o, &M, &classM);
    if (o->verify) { verify_dataset(P, M, classM, "M"); }
    double *HM = notnull(malloc(M->NS*sizeof(double)), "no mem");
    compute_handicaps(M, classM, o->opf, FALSE, HM);
    if (o->image_size > 0) 
      { output_dataset_image(P, M, classM, HM, "M-ini", o);
        output_nn_images(P, M, classM, HM,   "M-ini", o);
        if (o->opf) { output_nn_hand_cmp_image(P, M, classM, HM,   "M-ini", o); }
      }
     
    /* Generate the evaluation dataset (grid or random): */
    dataset_t *E = NULL;
    int *classE = NULL;
    if (o->samples_eval > 0)
      { generate_dataset(P, o->attributes_irr, o->samples_eval, o->grid, o, &E, &classE);
        if (o->verify) { verify_dataset(P, E, classE, "E"); }
        evaluate_nn_classifier(P->NC, M, classM, HM, "M-ini", E, classE, "E");
        if (o->image_size > 0) { output_dataset_image(P, E, classE, NULL, "E-ini", o); }
      }
    
    /* Output images of the ideal and the unrefined 1NN classifiers: */
    
    
    if (o->samples_refine > 0)
      { /* Generate the refinement dataset (random): */
        dataset_t *R;
        int *classR;
        generate_dataset(P, o->attributes_irr, o->samples_refine, FALSE, o, &R, &classR);
        if (o->verify) { verify_dataset(P, R, classR, "R"); }
        if (o->image_size > 0) 
          { output_dataset_image(P, R, classR, NULL, "R-ini", o); }

        /* Refine model and output images of the result: */
        refine_nn_model(P->NC, M, classM, HM, R, classR, o, TRUE); 
        
        if (E != NULL)
          { evaluate_nn_classifier(P->NC, M, classM, HM, "M-ref", E, classE, "E"); }

        /* Output images of the refined 1NN classifier: */
        if (o->image_size > 0) 
          { output_nn_images(P, M, classM, HM,   "M-ref", o);
            output_dataset_image(P, M, classM, HM, "M-ref", o);
            if (o->opf) { output_nn_hand_cmp_image(P, M, classM, HM,   "M-ref", o); }
          }
      }
    
    return 0;
  }
  
void compute_handicaps(dataset_t *M, int classM[], bool_t opf, bool_t verbose, double HM[])
  {
    if (opf)
      { rn_classif_opf_compute_handicaps(M, classM, rn_dist, HM, verbose); }
    else
      { int i; for (i = 0; i < M->NS; i++) { HM[i] = 0; } }
  }

void generate_dataset(problem_t *P, int NAI, int NS, bool_t grid, options_t *o, dataset_t **DP, int **classDP)
  { 
    int NAR = P->NA;
    int NA = NAR + NAI;
    if (grid)
      { generate_raw_dataset_grid(P, NA, NS, o->seed, DP, classDP); }
    else
      { generate_raw_dataset_random(P, NA, NS, o->seed, o->verify, DP, classDP); }
    dataset_t *D = (*DP);   
    int *classD = (*classDP);
    dataset_stats_print(stderr, D->NS, P->NC, classD);
    append_irrelevant_attributes(D, P, o->seed);
    if (o->noise > 0) { add_random_noise(D, o->noise, o->seed); } 
  }

void dataset_stats_print(FILE *wr, int NS, int NC, int classD[])
  {
    int *num = notnull(malloc((NC+1)*sizeof(int)), "no mem");
    rn_classif_class_count(NS, classD, NC, num); 
    rn_classif_class_count_print(wr, NC, num);
    free(num);
  }

void output_dataset(dataset_t *D, int classD[], int NC, options_t *o)
  { int cl;
    for (cl = 1; cl <= NC; cl++)
      { char *fname = NULL;
        asprintf(&fname, "%s-c%03d.dat", o->prefix, cl);
        FILE *wr = open_write(fname, TRUE);
        rn_classif_dataset_write(wr, D, cl, classD);
        fclose(wr);
        free(fname);
      }
  }

void verify_dataset(problem_t *P, dataset_t *D, int classD[], char *nameD)
  { int *classDX = notnull(malloc(D->NS*sizeof(int)), "no mem");
    rn_classif_dataset_label(D, P->NA, P->NC, P->lab, classDX);
    int ns, nf;
    rn_classif_compare(D->NS, classD, classDX, &ns, &nf);
    if (nf != 0)
      { fprintf(stderr, "\n");
        fprintf(stderr, "=== cross-classification of ideal vs 1NN(%s) on set %s ===\n", nameD, nameD);
        int *ncc;
        rn_classif_cross_matrix_build(D->NS, classD, P->NC, classDX, P->NC, &ncc);
        rn_classif_cross_matrix_print(stderr, P->NC, P->NC, ncc, TRUE);
        free(ncc);
        fprintf(stderr, "\n");
      }
    free(classDX);
  }

int evaluate_nn_classifier(int NC, dataset_t *M, int classM[], double HM[], char *nameM, dataset_t *E, int classE[], char *nameE)
  {
    int *classEX = notnull(malloc(E->NS*sizeof(int)), "no mem");
    rn_classif_nn_label_dataset(M, HM, classM, E, rn_dist, classEX);
    
    fprintf(stderr, "=== cross-classification of 1-NN(%s) on evaluation set %s ===\n", nameM, nameE);
    int *ncc;
    rn_classif_cross_matrix_build(E->NS, classE, NC, classEX, NC, &ncc);
    rn_classif_cross_matrix_print(stderr, NC, NC, ncc, TRUE);
    free(ncc);
    fprintf(stderr, "\n");
    
    int ns, nf;
    rn_classif_compare(E->NS, classE, classEX, &ns, &nf);
    free(classEX);
    return nf;
  }
  
void refine_nn_model(int NC, dataset_t *M, int classM[], double HM[], dataset_t *R, int classR[], options_t *o, bool_t verbose)
  {
    auto void hproc(dataset_t *D, int classD[], double HD[]);
    
    int *classRX = notnull(malloc(R->NS*sizeof(int)), "no mem");
    int maxIters = o->refine_iters;
    if (maxIters <= 0) { maxIters = (int)imin(3*M->NS, 3*R->NS); }
    rn_classif_nn_improve(M, HM, classM, R, classR, classRX, rn_dist, hproc, maxIters, o->swap_same, verbose);
    
    /* Local procedure implementation: */
    
    void hproc(dataset_t *D, int classD[], double HD[])
      { compute_handicaps(D, classD, o->opf, FALSE, HD); }
  }

void output_ideal_image(problem_t *P, char *tag, options_t *o)
  {
    r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
    double HV = 1.0 + 4*o->noise;  /* Half-extent of imaged area. */
    srandom(o->seed + 314159);
    
    fprintf(stderr, "computing the ideal classifier image '%s-C' ...\n", tag);
    uint16_image_t *imgC = rn_classif_uint16_image
      ( P->NA,  
        P->NC, P->lab,  
        0, NULL,  
        0, NULL, NULL, NULL, 
        &ctr, HV, o->image_size, o->image_subsmp, o->noise
      );
    output_image(o->prefix, tag, "C", imgC);
    uint16_image_free(imgC);
  }

void output_dataset_image(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o)
  { 
    r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
    double HV = 1.0 + 4*o->noise;  /* Half-extent of imaged area. */
    srandom(o->seed + 314159);
    
    fprintf(stderr, "computing the dataset image '%s-P' ...\n", tag);
    uint16_image_t *imgC = 
      rn_classif_uint16_image
      ( P->NA,  
        0, NULL, 
        P->NC, P->lab,  
        P->NC, D, HD, classD,  
        &ctr, HV, o->image_size, o->image_subsmp, o->noise
      );
    output_image(o->prefix, tag, "P", imgC);
    uint16_image_free(imgC);
  }

void output_nn_images(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o)
  { 
    auto int evalDx(int NAX, int NCX, double px[]);
      /* The 1NN classifier {(D,classD,HD)}, with {px} already expanded
        from {NA} attributes to {D.NA} attributes.  Requires {NCX==P.NC}, {NAX==D.NA}; 
        assumes that {classD} ranges over {0..P.NC}. */

    auto int evalD(int NAX, int NCX, double p[]);
      /* The 1NN classifier {(D,classD,HD)}, that expands {p} with random {U} values
        from {NA} attributes to {D.NA} attributes.  Requires {NCX==P.NC}, {NAX==P.NA}
        and {P.NA <= D.NA}; assumes that {classD} ranges over {0..P.NC}. */

    auto int diffPD(int NAX, int NCX, double p[]);
      /* A labeler that compares the classifier {P} on {p[0..NA-1}
        with the 1NN {evalD} classifier above. If {P} assigns class 0,
        returns 0. Otherwise, if the two classes agree, returns 2.
        Otherwise returns 1. Requires {NCX==2} and {NAX==P.NA}. */
    
    r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
    double HV = 1.0 + 4*o->noise;  /* Half-extent of imaged area. */
    
    fprintf(stderr, "computing the 1NN classifier image '%s-C' ...\n", tag);
    srandom(o->seed + 314159);
    uint16_image_t *imgC = 
      rn_classif_uint16_image
        ( P->NA,  
          P->NC, evalD,  
          P->NC, P->lab,  
          P->NC, D, HD, classD,  
          &ctr, HV, o->image_size, o->image_subsmp, o->noise
        );
    output_image(o->prefix, tag, "C", imgC);
    uint16_image_free(imgC);
    
    fprintf(stderr, "computing the ideal/1NN classifier discrepancy image '%s-D' ...\n", tag);
    srandom(o->seed + 314159);
    uint16_image_t *imgD = 
      rn_classif_uint16_image
        ( P->NA,  
          2, diffPD,  
          0, NULL,  
          0, D, HD, NULL,  
          &ctr, HV, o->image_size, o->image_subsmp, o->noise
        );
    output_image(o->prefix, tag, "D", imgD);
    uint16_image_free(imgD);

    int evalDx(int NAX, int NCX, double px[])
      { assert(NAX == D->NA);
        int imin = rn_classif_find_nearest_in_dataset(D, HD, D->NA, px, rn_dist);
        assert((imin >= 0) && (imin <= D->NS));
        int cl = classD[imin];
        assert((cl >= 0) && (cl <= P->NC));
        return cl;
      }
    
    int evalD(int NAX, int NCX, double p[])
      { assert(NAX == P->NA);
        assert(NAX <= D->NA);
        /* Expand {p[0..NA-1]} to {px[0..D.NA-1]} by appending random values: */
        double px[D->NA];
        int t;  for (t = 0; t < D->NA; t++) { px[t] = (t < D->NA ? p[t] : 2*drandom() - 1); }
        /* Classify {px} with {(D,HD,classD)}: */
        return evalDx(D->NA, NCX, px);
      }
    
    int diffPD(int NAX, int NCX, double p[])
      { assert(NAX == P->NA);
        assert(NCX == 2);
        int clP = P->lab(NAX, P->NC, p);
        if (clP == 0)
          { return 0;}
        else
          { int clD = evalD(NAX, P->NC, p);
            return (clD == clP ? 2 : 1);
          }
      }
  }

void output_nn_hand_cmp_image(problem_t *P, dataset_t *D, int classD[], double HD[], char *tag, options_t *o)
  { 
    auto int evalDx(int NAX, int NCX, double HX[], double px[]);
      /* The 1NN classifier {(D,classD,HD)}, with {px} already
        expanded from {NA} attributes to {D.NA} attributes. The
        handicap tbale HX may be NULL. Requires {NCX==P.NC},
        {NAX==D.NA}; assumes that {classD} ranges over {0..P.NC}. */

    auto int diffH(int NAX, int NCX, double p[]);
      /* A labeler that compares the 1NN classifier {(D,HD,classD)}
        with {(D,NULL,classD)}. If {P} returns 0, returns 0.
        Otherwise, if the handicaps were an iprovement, returns 3. If
        they were a disimprovement, returns 1. If they made no
        difference, returns 2. If {Requires {NCX==3} and
        {NAX==P.NA}. */
    
    r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
    double HV = 1.0 + 4*o->noise;  /* Half-extent of imaged area. */
    srandom(o->seed + 314159);
    
    fprintf(stderr, "computing the 1NN classifier with/sans H discrepancy image '%s-G' ...\n", tag);
    srandom(o->seed + 314159);
    uint16_image_t *imgG = 
      rn_classif_uint16_image
      ( P->NA,  
        3, diffH, 
        P->NC, P->lab,  
        0, D, HD, NULL,  
        &ctr, HV, o->image_size, o->image_subsmp, o->noise
      );
    output_image(o->prefix, tag, "G", imgG);
    uint16_image_free(imgG);

    int evalDx(int NAX, int NCX, double HX[], double px[])
      { assert(NAX == D->NA);
        int imin = rn_classif_find_nearest_in_dataset(D, HX, D->NA, px, rn_dist);
        assert((imin >= 0) && (imin <= D->NS));
        int cl = classD[imin];
        assert((cl >= 0) && (cl <= P->NC));
        return cl;
      }
    
    int diffH(int NAX, int NCX, double p[])
      { assert(NAX == P->NA);
        assert(NCX == 3);
        int clP = P->lab(NAX, P->NC, p);
        if (clP == 0)
          { return 0;}
        else
          { /* Expand {p[0..NA-1]} to {px[0..D.NA-1]} by appending random {U} values: */
            double px[D->NA];
            int t;  for (t = 0; t < D->NA; t++) { px[t] = (t < D->NA ? p[t] : 2*drandom() - 1); }
            /* Classify {px} with and without {HD}: */
            int cl_with_H = evalDx(D->NA, P->NC, HD, px);
            int cl_sans_H = evalDx(D->NA, P->NC, NULL, px);
            /* Compare classifications: */
            if (cl_with_H == cl_sans_H)
              { if (cl_with_H == clP)
                  { /* All three agree: */ return 2; }
                else
                  { /* With/sans agree but both are wrong: */ return 2; }
              }
            else
              { if (cl_with_H == clP)
                  { /* Handicaps improved: */  return 3; }
                else if (cl_sans_H == clP)
                  { /* Handicaps made it worse: */  return 1; }
                else
                  { /* Handicap changed but did not improve: */ return 2; }
              }
          }
      }
  }

void output_image(char *prefix, char *tag, char *kind, uint16_image_t *img)
  {
    char *fname = NULL;
    asprintf(&fname, "%s-%s-%s.%s", prefix, tag, kind, (img->chns == 1 ? "pgm" : "ppm"));
    uint16_image_write_pnm_named(fname, img, FALSE, TRUE);
    free(fname);
  }

void append_irrelevant_attributes(dataset_t *D, problem_t *P, uint32_t seed)
  {
    int NS = D->NS;
    int NAR = P->NA;
    int NA = D->NA;
    /* Append the irrelevant attributes {NAR..NA-1}: */
    /* Must do this separately so that the essential ones are fixed for a given {seed}. */
    srandom(seed + 4615);
    int i, t;
    for (i = 0; i < NS; i++)
      { double *attri = D->smp[i];
        for (t = NAR; t < NA; t++) { attri[t] = 2*drandom() - 1; }
      }
  }

void add_random_noise(dataset_t *D, double sigma, uint32_t seed)
  {
    int NS = D->NS;
    int NA = D->NA;
    srandom(seed + 19501129);
    int i, t;
    for (i = 0; i < NS; i++)
      { double *attri = D->smp[i];
        for (t = 0; t < NA; t++) 
          { /* Add to {attri[t]} a Gaussian deviate truncated to {4*sigma} */
            double dit = sigma*fmax(-4.0, fmin(+4.0, dgaussrand()));
            attri[t] += dit;
          }
      }
  }
     
void generate_raw_dataset_random(problem_t *P, int NA, int NS, uint32_t seed, bool_t verify, dataset_t **DP, int **classDP)
  {
    int NC = P->NC;
    int NAR = P->NA;
    dataset_t *D = rn_classif_dataset_new(NS, NA);
    int *classD = notnull(malloc(NS*sizeof(int)), "no mem");

    /*  Generate samples: */
    srandom(seed);
    int nerr = 0;
    int i;
    for (i = 0; i < NS; i++)
      { double *attri = notnull(malloc(NA*sizeof(double)), "no mem");
        D->smp[i] = attri;
        P->gen(i, NAR, NC, attri, &(classD[i]));
        if (verify) 
          { int ver = P->lab(NAR, NC, attri);
            if (classD[i] != ver)
              { /* Since sample is unperturbed,this is a fatal error or very bad luck: */
                fprintf(stderr, "** inconsistent class %d != %d for sample %d:\n", classD[i], ver, i);
                rn_classif_sample_print(stderr, "  ", NAR, attri, " ", "\n");
                nerr++;
                if (nerr > 100) { fprintf(stderr, "** too many errors\n"); exit(1); }
              }
          }
      }
    (*DP) = D;
    (*classDP) = classD;
  }
  
void generate_raw_dataset_grid(problem_t *P, int NA, int NS, uint32_t seed, dataset_t **DP, int **classP)
  {
    int NC = P->NC;
    int NAR = P->NA;
    
    /* Compute the grid order {NG} such that {NG^NAR >= NS}: */
    int NG = (int)ceil(pow(NS, 1.0/NAR));
    int NT = ipow(NG, NAR); /* Total number of tentative samples. */
    assert(NT >= NS);
    demand(NT <= MAX_SAMPLES, "too many samples in grid");
    
    /* Allocate sample array with maximum size: */
    int *classD = notnull(malloc(NT*sizeof(int)), "no mem");
    dataset_t *D = rn_classif_dataset_new(NT, NA);
    
    /* Each sample is identified by a tuple of {NAR} indices in {0..NG-1}: */
    int ix[NAR];
    
    /* Generate grid samples and store those that have nonzero class. */
    /* Those samples are stored in postions {0..NOK-1} of {attr} and {classD}. */
    NS = 0;
    srandom(seed);
    double attr[NA]; /* Temporary attribute vector. */
    int g;
    for (g = 0; g < NT; g++)
      { /* Break {g} down into {NAR} grid indices: */
        get_grid_indices(g, NG, NAR, ix);
        /* Convert the indices into attributes in {U}: */
        int t;
        for (t = 0; t < NAR; t++) { attr[t] = 2*(ix[t] + 0.5)/NG - 1; }
        /* Now use the problem's classifier to obtain the classD: */
        int cl = P->lab(NAR, NC, attr);
        if (cl != 0)
          { /* Sample is inside some domain, store it: */
            double *attrn = notnull(malloc(NA*sizeof(double)), "no mem");
            for (t = 0; t < NAR; t++) { attrn[t] = attr[t]; }
            D->smp[NS] = attrn;
            classD[NS] = cl;
            NS++;
          }
      }
      
    /* Trim excess storage and return: */
    D->NS = NS;
    D->smp = notnull(realloc(D->smp, NS*sizeof(double*)), "no mem");
    classD = notnull(realloc(classD, NS*sizeof(int)), "no mem");
    (*DP) = D;
    (*classP) = classD;
  }
  
void get_grid_indices(int g, int NG, int NAR, int ix[])
  {
    int k;
    int tmp = g;
    for (k = 0; k < NAR; k++) { ix[k] = tmp % NG; tmp = tmp/NG; } 
  }
  
problem_t *get_problem(char *name, int NAR, int NC)
  {
    problem_kind_t kind;
    for (kind = 0; kind < problem_kind_NUMBER; kind++) 
      { if (strcmp(name, problem_kind_name[kind]) == 0)
          { /* Found the problem: */
            /* Get the problem-specific procedure {def}: */
            problem_def_proc_t *chk = NULL;
            rn_classif_labeler_t *lab = NULL;
            rn_classif_thrower_t *gen = NULL;
            switch(kind)
              {
                case problem_kind_SATURN: 
                  chk = &rn_classif_test_check_saturn;
                  lab = &rn_classif_test_label_saturn;
                  gen = &rn_classif_test_throw_saturn;
                  break;
                case problem_kind_PETALS: 
                  chk = &rn_classif_test_check_petals;
                  lab = &rn_classif_test_label_petals;
                  gen = &rn_classif_test_throw_petals;
                  break;
                case problem_kind_VESSEL: 
                  chk = &rn_classif_test_check_vessel;
                  lab = &rn_classif_test_label_vessel;
                  gen = &rn_classif_test_throw_vessel;
                  break;
                case problem_kind_MBALLS: 
                  chk = &rn_classif_test_check_mballs;
                  lab = &rn_classif_test_label_mballs;
                  gen = &rn_classif_test_throw_mballs;
                  break;
                case problem_kind_SHELLS: 
                  chk = &rn_classif_test_check_shells;
                  lab = &rn_classif_test_label_shells;
                  gen = &rn_classif_test_throw_shells;
                  break;
                case problem_kind_STAOLI: 
                  chk = &rn_classif_test_check_staoli;
                  lab = &rn_classif_test_label_staoli;
                  gen = &rn_classif_test_throw_staoli;
                  break;
                default: assert(FALSE);
              }
            problem_t *P = notnull(malloc(sizeof(problem_t)), "no mem");
            (*P) = (problem_t) { .NA = NAR, .NC = NC, .lab = lab, .gen = gen };
            chk(&(P->NA), &(P->NC));
            return P;
          }
      }
    demand(FALSE, "unrecognized problem kind");
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-problem");
    o->problem = argparser_get_next(pp);

    argparser_get_keyword(pp, "-samples");
    o->samples_model = (int)argparser_get_next_int(pp, 1, MAX_SAMPLES);
    o->samples_refine = (int)argparser_get_next_int(pp, 0, MAX_SAMPLES);
    o->samples_eval = (int)argparser_get_next_int(pp, 0, MAX_SAMPLES);

    argparser_get_keyword(pp, "-seed");
    o->seed = (uint32_t)argparser_get_next_uint(pp, 1, MAX_SEED);
    
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp);
    
    if (argparser_keyword_present(pp, "-attributes"))
      { o->attributes_rel = (int)argparser_get_next_int(pp, 0, MAX_ATTRIBS);
        o->attributes_irr = (int)argparser_get_next_int(pp, 0, MAX_ATTRIBS - o->attributes_rel); 
      }
    else
      { o->attributes_rel = 0; o->attributes_irr = 0; }

    if (argparser_keyword_present(pp, "-classes"))
      { o->classes = (int)argparser_get_next_int(pp, 0, MAX_CLASSES); }
    else
      { o->classes = 0; }

    if (argparser_keyword_present(pp, "-noise"))
      { o->noise = argparser_get_next_double(pp, 0, MAX_NOISE); }
    else
      { o->noise = 0; }

    if (argparser_keyword_present(pp, "-verify"))
      { o->verify = argparser_get_next_bool(pp); }
    else
      { o->verify = FALSE; }

    if (argparser_keyword_present(pp, "-opf"))
      { o->opf = argparser_get_next_bool(pp); }
    else
      { o->opf = FALSE; }

    if (argparser_keyword_present(pp, "-grid"))
      { o->grid = argparser_get_next_bool(pp); }
    else
      { o->grid = FALSE; }

    if (argparser_keyword_present(pp, "-refine"))
      { o->refine_iters = (int)argparser_get_next_int(pp, 0, MAX_REFINE_ITERS); }
    else
      { o->refine_iters = 0; }

    if (argparser_keyword_present(pp, "-swap"))
      { if (argparser_keyword_present_next(pp, "same"))
          { o->swap_same = TRUE; }
        else if (argparser_keyword_present_next(pp, "full"))
          { o->swap_same = FALSE; }
        else   
          { argparser_error(pp, "expecting \"same\" or \"full\""); }
      }
    else
      { o->swap_same = FALSE; }

    if (argparser_keyword_present(pp, "-image"))
      { o->image_size = (int)argparser_get_next_int(pp, 0, MAX_IMAGE_SIZE);
        o->image_subsmp = (int)argparser_get_next_int(pp, 0, MAX_SUBSMP);
      }
    else
      { o->image_size = 0;
        o->image_subsmp = 0;
      }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
