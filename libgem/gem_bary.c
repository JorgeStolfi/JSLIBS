/* See {gem_bary.h}  */
/* Last edited on 2024-12-05 10:25:56 by stolfi */

#include <affirm.h>

#include <gem.h>

#include <gem_bary.h>

#define gem_bary_INI_NODES 500  
  /* Initial size of node stacks. */

void gem_bary_splice(gem_ref_t a, gem_ref_t b, int k)
  {
    gem_ref_vec_t nodesA = gem_ref_vec_new(gem_bary_INI_NODES);
    gem_ref_vec_t nodesB = gem_ref_vec_new(gem_bary_INI_NODES);
    int nnA = 0, nnB = 0;

    if (k == 0)
      { gem_splice(a, b, 1); }
    else
      { gem_traverse(a, k-1, &nodesA, &nnA);
        gem_traverse(b, k-1, &nodesB, &nnB);
        demand(nnA == nnB, "faces do not match");

        int i;
        for (i = 0; i < nnA; i++) { gem_splice(nodesA.e[i], nodesB.e[i], k+1); }
      }
  }
