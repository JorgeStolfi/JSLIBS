/* See dg_basis.h */
/* Last edited on 2009-05-18 14:21:41 by stolfi */

#include <math.h>

#include <affirm.h> 
#include <vec.h> 

#include <mdg_grid.h> 
#include <dg_tree.h>
#include <dg_locus.h>
#include <dg_tent.h> 

#include <dg_tent_basis.h> 

/* INTERNAL PROTOTYPES */

void dg_tent_basis_gather_templates
  ( mdg_dim_t d,                    /* Grid dimension. */
    udg_pulse_family_t fam[],       /* Tent family. */
    int nmt,                       /* Number of mother tents in family. */
    int *ntpp,                     /* (OUT) number of distinct templates. */
    mdg_grid_size_t rsz[],          /* (OUT) {rsz[k*d + i]} is size of template {k} along axis {i}. */
    dg_tent_mother_index_t ftix[], /* (OUT) smallest mother tent index with given template. */
    dg_tent_mother_index_t ntix[]  /* (OUT) next mother tent index with same template. */
  );
  /* Collects the distinct support templates for all {nmt} mother
    tents of tent family {fam[0..d-1]}.
    
    The function will store into {*ntpp} the number {ntp} of distinct
    templates found.
    
    For each {k} in {0..ntp-1}, the procedure will store template
    number {k} into {rsz[k*d + i]} for {i} in {0..d-1}. It will also
    store into {ftix[k]} the index of the first mother tent with
    that support template. 
    
    Finally, for each mother tent index {tix}, it will store into
    {ntix[tix]} the index of the next mother tent with same
    support template as mother tent {tix}; or
    {dg_tent_mother_index_NONE} if {tix} is the last mother tent with
    that support. */

/* IMPLEMENTATIONS */

dg_tent_vec_t dg_tent_basis_get
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[], 
    dg_tree_t G,
    bool_t superior
  )
  {
    /* Allocate the basis, initially empty: */
    dg_tent_vec_t B = dg_tent_vec_new(0);
    int nB = 0;
    
    /* Compute the number of distinct tent indices: */
    int nmt = dg_tent_mother_index_max(d, fam) + 1;
    
    /* Get all distinct templates and their tent indices. */
    int ntp;                   /* Number of ditinct templates (at most {nmt}). */
    mdg_grid_size_t rsz[nmt*d];        /* {rsz[d*k + i]} is the size of template {k} along axis {i}. */
    dg_tent_mother_index_t ftix[nmt]; /* {ftix[k]} is min mother tent index with template {k}. */
    dg_tent_mother_index_t ntix[nmt]; /* {ntix[tix]} is next mother tent w/ same template as {tix}. */
    dg_tent_basis_gather_templates(d, fam, nmt, &ntp, rsz, ftix, ntix);
 
    int k;
    for (k = 0; k < ntp; k++)
      { /* Get hold of template number {k}: */
        mdg_grid_size_t *rszk = &(rsz[k*d]);
        /* Enumerate all instances of template {rszk[0..d-1]} in the tree, and output the tents: */

        auto void process_template_instance(dg_tree_vec_t *N, bool_t is_leaf);
          /* Given an instance {N} of the template {rszk} in {G},
            appends to {B.e[0..nB-1]} one copy of each mother tent from family {fam}
            that has support size {rszk} --- scaled, translated, and folded so that 
            their support is {N}.
            
            The procedure does nothing, however, if those tents would
            create linear dependencies.  Which tents get eliminated 
            depends on the {superior} flag. */

        void process_template_instance(dg_tree_vec_t *N, bool_t is_leaf)
          {
            /* Decide whether this instance is appropriate: */
            bool_t take;
            if (superior)
              { take = dg_tree_pack_is_even(N); }
            else
              { take = dg_tree_pack_is_leaf(N); }
            if (take)
              { /* Enumerate all tents with this template: */
                dg_tent_mother_index_t tix = ftix[k];
                while (tix != dg_tent_mother_index_NONE)
                  { /* Output a tent with mother index {tix} and position {N[0]}: */
                    dg_tent_vec_expand(&B, nB);
                    B.e[nB] = (dg_tent_t){ N->e[0]->index, tix };
                    nB++;
                    tix = ntix[tix];
                  }
              }
          }
        
        dg_tree_pack_enum(d, G, rszk, process_template_instance);
      }
      
    dg_tent_vec_trim(&B, nB);
    return B;
  }

void dg_tent_basis_gather_templates
  ( mdg_dim_t d,                    /* Grid dimension. */
    udg_pulse_family_t fam[],       /* Tent family. */
    int nmt,                       /* Number of mother tents in family. */
    int *ntpp,                     /* (OUT) number of distinct templates. */
    mdg_grid_size_t rsz[],          /* (OUT) {rsz[k*d + i]} is size of template {k} along axis {i}. */
    dg_tent_mother_index_t ftix[], /* (OUT) smallest mother tent index with given template. */
    dg_tent_mother_index_t ntix[]  /* (OUT) next mother tent index with same template. */
  )
  {
    int tix;
    int ntp = 0;
    int k, i;
    for (tix = 0; tix < nmt; tix++)
      { /* Get mother pulse indices {pix[0..d-1]} from mother tent index {tix}: */
        udg_pulse_mother_index_t pix[d];
        dg_tent_mother_index_unpack(d, fam, tix, pix);

        /* Get the template of mother tent {tix}: */
        mdg_grid_size_t msz[d];
        dg_tent_mother_supp_count(d, fam, pix, msz);
        
        /* Check if this template has occurred before: */
        bool_t matched = FALSE;
        for (k = 0; (k < ntp) && (! matched); k++)
          { /* Check if {msz[0..d-1]} matches {rsz[k*d..k*d+d-1]}: */
            mdg_grid_size_t *rszk = &(rsz[k*d]);
            matched = TRUE;
            for (i = 0; i < d; i++)  { matched &= (msz[i] == rszk[i]); }
            if (matched) { /* Exit loop on {k}: */ break; }
          }
        if (! matched)
          { /* Save this new template, set {k} to its index: */
            k = ntp; ntp++;
            mdg_grid_size_t *rszk = &(rsz[k*d]);
            for (i = 0; i < d; i++) { rszk[i] = msz[i]; }
            ftix[k] = dg_tent_mother_index_NONE; /* For now. */            
          }
        /* Associate template index {k} to the tent index {tix}: */
        ntix[tix] = ftix[k];
        ftix[k] = tix;
      }
    /* Reverse the lists: */
    for (k = 0; k < ntp; k++)
      { dg_tent_mother_index_t atix = ftix[k]; 
        dg_tent_mother_index_t rtix = dg_tent_mother_index_NONE;
        while (atix != dg_tent_mother_index_NONE)
          { tix = ntix[atix]; 
            ntix[atix] = rtix; 
            rtix = atix; atix = tix; 
          }
        ftix[k] = rtix;
      }

    /* Return the number of templates found; */
    (*ntpp) = ntp;
  }
