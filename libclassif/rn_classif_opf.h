/* rn_classif_opf.h --- tools for the optimum path forest classifier. */
/* Last edited on 2024-12-05 10:24:27 by stolfi */

#ifndef rn_classif_opf_H
#define rn_classif_opf_H

#include <stdio.h>
#include <math.h>

#include <r2.h>
#include <rn_classif.h>
  
void rn_classif_opf_compute_handicaps
  ( rn_classif_dataset_t *M,    /* Model samples. */
    int classM[],               /* {classM[i]} is the given class of {M.smp[i]}. */
    rn_classif_pq_dist_t *dist, /* Sample distance function. */
    double H[],                 /* (OUT) {H[i]} is the OPF handicap of sample {M.smp[i]}. */
    bool_t verbose              /* TRUE for diagnostic messages. */
  );
  /* Computes the OPF handicaps {H[0..M.NS-1]} for the 1NN 
    classifier {(M,H,classM)} as described in J. P. Papa, A. X. Falcao, 
    C. T. N. Suzuki "Supervised Classification ..." (IJIST, 2009). */

#endif
