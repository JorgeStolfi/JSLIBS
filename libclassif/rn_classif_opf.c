/* See rn_classif_opf.h. */
/* Last edited on 2024-12-05 10:24:24 by stolfi */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rn_classif.h>
#include <rn_classif_opf.h>
#include <opf.h>
#include <mst.h>

void rn_classif_opf_compute_handicaps
  ( rn_classif_dataset_t *M,    /* Model samples. */
    int classM[],               /* {classM[i]} is the given class of {M.smp[i]}. */
    rn_classif_pq_dist_t *dist, /* Sample distance function. */
    double H[],                 /* (OUT) {H[i]} is the OPF handicap of sample {M.smp[i]}. */
    bool_t verbose              /* TRUE for diagnostic messages. */
  )
  {
    int NS = M->NS;
    int NA = M->NA;

    auto double ijdist_full(int i, int j);
      /* Returns {dist(M.smp[i],M.smp[j])}, except that fudges 0 to a tiny value if {i!=j}. */
    
    auto double ijdist_homo(int i, int j);
      /* If {i} and {j} have the same class, returns {ijdist_full(i,j)}, else {+INF}. */
    
    /* Output vectors for {opf_build_complete}: */
    int *P = notnull(malloc(NS*sizeof(int)), "no mem"); /* Predecessor map. */
    int *R = notnull(malloc(NS*sizeof(int)), "no mem"); /* Root map. */

    /* Build a MST of the complete graph {K} on {M} with the {dist} arc cost. */
    int i;
    for (i = 0; i < NS; i++) { H[i] = +INF; }
    H[0] = 0;
    mst_build_complete(NS, ijdist_full, P, H, verbose);
    if (verbose) { rn_classif_nn_print_handicaps(stderr, "mst", NS, H); }

    /* Build an OPF for the graph {G} on {M} with same-class edges only. */
    /* Restrict the seeds to the vertices with MST neighbors in a different class. */
    for (i = 0; i < NS; i++) { H[i] = +INF; }
    for (i = 0; i < NS; i++)
      { int j = P[i];
        if ((j != i) && (classM[i] != classM[j])) { H[i] = 0; H[j] = 0; }
      }
    /* Build the OPF of {M}, with {fmax} path cost, within each class only: */
    opf_build_complete(NS, H, ijdist_homo, fmax, P, R, verbose);
    if (verbose) { rn_classif_nn_print_handicaps(stderr, "opf", NS, H); }
    free(P);
    free(R);

    /* Local procedure implementations */
    
    double ijdist_full(int i, int j)
      { double d = dist(NA, M->smp[i], M->smp[j]);
        assert((! isnan(d)) && (d >= 0));
        if ((i != j) && (d == 0)) { d = 1.0e-200; }
        return d;
      }
    
    double ijdist_homo(int i, int j)
      { if (classM[i] != classM[j])
          { return +INF;}
        else
          { return ijdist_full(i,j); }
      }
  }
  
void rn_classif_nn_print_handicaps(FILE *wr, char *title, int NS, double H[])
  { 
    fprintf(wr, "%s handicaps:\n", title); 
    int i;
    for (i = 0; i < NS; i++) {
      fprintf(wr, "  %5d %8.6f\n", i, H[i]); 
    }
    fprintf(wr, "\n");
  }

