/* rn_classif_test.h --- generate random points for testing classifiers. */
/* Last edited on 2024-12-05 10:24:32 by stolfi */

#ifndef rn_classif_test_H
#define rn_classif_test_H

#include <stdio.h>
#include <math.h>

#include <r2.h>
  
/* 
  The procedures in this interface are intended to generate random
  samples for testing point classification methods, such as Nearest
  Neighbor.  
  
  For each values of the parameters {NC} and {NA}, each procedure
  {rn_classif_test_label_{XXX}} or {} defines {NC} pairwise disjoint /class
  domains/ contained in the hypercube {U^NA=[-1_+1]^NA} each
  identified with a /class index/ in the range {1..NC}. For each
  problem there is a set of allowed values of {NA} and {NC}, a
  specific set of class domains (that usually depend on {NC}), and a
  specific sampling probability distribution inside each domain. The
  domains may not cover the whole of {U^NA} and may or may not touch
  each other. Some domains are . zero-measure sets with dimension
  smaller than {NA}.

  The procedure {rn_classif_test_label_{XXX}} finds the
  class domain that contains the given point {*p}, and returns
  the index of that domain.  If {*p} is not in any domain but 
  is within a small distance of some class domain, may return that class.
  If {*p} is too far from any class domain, it returns.
  
  The procedure {rn_classif_test_throw_{XXX}} generates samples
  from the same domains. It first computes a class index {*class} as some
  periodic function of the parameter {i}, and choses a point {*p} randomly inside the
  domain with index {*class}.  The input values of {*p} and {*class} are ignored.

  In some cases, the class index is simply {(i % NC) + 1}. In many
  cases, the mapping of {i} to {*class} is chosen so that each domain
  gets chosen with a frequency approximately proportional to its measure.

  The companion procedure {rn_classif_test_check_{XXX}} checks whether
  the parameters {*NA} and {*NC} are valid, and provides suitable
  defaults if any of them are zero. */
  
void rn_classif_test_check_saturn(int *NA, int *NC);
int rn_classif_test_label_saturn(int NA, int NC, double p[]);
void rn_classif_test_throw_saturn(int i, int NA, int NC, double p[], int *class);
  /* Two thin elliptical rings, interscting on 4 places. The number of
    attributes {NA} must be 2, the number of classes {NC} must be 2, and
    {*class} alternates between 1 and 2. The domains are unidimensional;
    domain 1 is an ellipse, domain 2 a circle, both concentric and
    intertwined. This distribution approximates the "saturn" ("B9")
    dataset from Papa et al.(2009) */
    
void rn_classif_test_check_petals(int *NA, int *NC);
int rn_classif_test_label_petals(int NA, int NC, double p[]);
void rn_classif_test_throw_petals(int i, int NA, int NC, double p[], int *class);
  /* A `flower' with four teardrop-shaped `petals'. The number of
    attributes {NA} must be 2, the number of classes {NC} must be 4, and
    {*class} cycles between 1 and 4. The domains are teardrop-shaped
    regions radiating from the origin, narrow end inwards. The sampling
    is concentrated mostly towards the edges of each domain. The petals
    meet at the origin as four narrow angles. This distribution
    approximates the "petals" ("B10") dataset from Papa et al.(2009) */
    
void rn_classif_test_check_vessel(int *NA, int *NC);
int rn_classif_test_label_vessel(int NA, int NC, double p[]);
void rn_classif_test_throw_vessel(int i, int NA, int NC, double p[], int *class);
  /* An elliptical ribbon enclosing two boxes. The number of attributes
    {NA} must be 2, the number of classes {NC} must be 3, and {*class}
    cycles between 1 and 3. Domain 1 is unidimensional, along an almost
    circular ellipse. Domains 2 and 3 are two rectangles contained in
    the ellipse, sampled mostly near their edges. This distribution was
    called "boat" ("B11") in Papa et al. */
    
void rn_classif_test_check_mballs(int *NA, int *NC);
int rn_classif_test_label_mballs(int NA, int NC, double p[]);
void rn_classif_test_throw_mballs(int i, int NA, int NC, double p[], int *class);
  /* An arrangement of balls of alternating classes. The number of
    attributes{NA} may be any positive integer and the number of classes
    {NC} must be 2 or more. The class domains comprise an array of
    {(NC-1)^NA} balls nested inside {U^NA}. Each domain {2..NC} consists
    of {(NC-1)^(NA-1)} of those balls, interleaved with the balls of
    other domains. Domain 1 is {U^NA} minus those balls. */
    
void rn_classif_test_check_shells(int *NA, int *NC);
int rn_classif_test_label_shells(int NA, int NC, double p[]);
void rn_classif_test_throw_shells(int i, int NA, int NC, double p[], int *class);
  /* Concentric shells. The number of attributes{NA} may be any positive
    integer and the number of classes {NC} must be 2 or more. Domains
    {2..NC} are spherical shells centered on the origin. Domain 1
    consists of the ball inside the innermost shell and of the outermost
    shell. All domains have the same measure. */
    
void rn_classif_test_check_staoli(int *NA, int *NC);
int rn_classif_test_label_staoli(int NA, int NC, double p[]);
void rn_classif_test_throw_staoli(int i, int NA, int NC, double p[], int *class);
  /* A fat box and a thin box, like Stanlie and Ollie (or Penn and
    Teller). The number of attributes{NA} may be any positive integer
    and the number of classes {NC} must be 3. Domain 2 is a regular
    hypercube nested inside {U^NA}. Domain 3 is a plate with with one
    side {a} and the remaining sides {4*a}, disjoint and separated from
    domain 2. Domain 1 is the remaider of {U^NA}. Unlike their
    namesakes, classes 2 and 3 have the same measure. */

#endif
