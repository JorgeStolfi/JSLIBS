/* See pst_integrate.h */
/* Last edited on 2025-01-18 12:35:09 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <jsfile.h>
#include <r2.h>
#include <jswsize.h>
#include <rn.h>

#include <pst_basic.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_interpolate.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_integrate.h>

/* INTERNAL PROTOTYPES */

void pst_integrate_set_equations(pst_imgsys_t *S, float_image_t *G, float_image_t *W, bool_t verbose);
  /* Builds the equation list {S.eq} for the integration system {S}, given
    the gradient map {G} and its weight map {W}.*/
   
#define MAX_COEFFS pst_imgsys_MAX_COEFFS

#define FLUFF (1.0e-140)
  /* Practically zero. */

/* IMPLEMENTATIONS */

pst_imgsys_t* pst_integrate_build_system(float_image_t *G, bool_t verbose) 
  {
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "gradient map must have 3 channels");
    
    /* Get the size of the system: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    
    pst_imgsys_t *S = pst_imgsys_new_grid(NX_Z, NY_Z);
    
    /* Gather the equations {S->eq[0..S->N-1]} of the system: */
    pst_integrate_set_equations(S, G, W, verbose);
    
    return S;
  }

void pst_integrate_set_equations(pst_imgsys_t *S, float_image_t *G, bool_t verbose)
  {
    /* Get/check the sizes of the slope maps: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    assert(NC_G == 3);
    
    /* Get size of the height map: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    assert(S->N == NX_Z*NY_Z);
    
    /* Conceptually, we define four axial quadratic mismatch terms
      {qpo[x,y](z), qmo[x,y](z), qop[x,y](z), qom[x,y](z)} for each
      corner of the pixel grid, that depend on the height map {z}.
      Each term is the square of the difference between two 
      estimates of the height increment along one of the edges incident to {(x,y)}:
      one obtained by interpolation of the slope map {G}, and the 
      other obtained by numerical difference of the (unknown) heights at the two
      endpoints of that edge. 
      
        { qpo[x,y](z) = (dpo[x,y] - (-z[x,y]  + z[x+1,y]))^2 }

        { qmo[x,y](z) = (dmo[x,y] - (-z[x,y]  + z[x-1,y]))^2 }

        { qop[x,y](z) = (dop[x,y] - (-z[x,y]  + z[x,y+1]))^2 }

        { qom[x,y](z) = (dom[x,y] - (-z[x,y]  + z[x,y-1]))^2 }

      where {dpo[x,y], dmo[x,y], dop[x,y], dom[x,y]} are the
      interpolated slopes at the edge midpoints {(x+1/2,y), (x-1/2,y),
      (x,y+1/2), (x,y-1/2)}, respectively. These slopes have
      interpolated reliability weights {wpo[x,y], wmo[x,y], wop[x,y],
      wom[x,y]}, respectively, which are defined to be the weights of
      the corresponding axial mismatch terms.
      
      We also define four skew mismatch terms
      {qmm[x,y](z), qmp[x,y](z), qpm[x,y](z), qpp[x,y](z)} 
      
        { qmm[x,y](z) = (dmm[x,y] - (-z[x,y]  + z[x-1,y-1]))^2 }
      
        { qmp[x,y](z) = (dmp[x,y] - (-z[x,y]  + z[x-1,y+1]))^2 }
      
        { qpm[x,y](z) = (dpm[x,y] - (-z[x,y]  + z[x+1,y-1]))^2 }
      
        { qpp[x,y](z) = (dpp[x,y] - (-z[x,y]  + z[x+1,y+1]))^2 }
      
      where {dmm[x,y], dmp[x,y], dpm[x,y], dpp[x,y]} are the dot product
      of the gradient interpolated at
      the skew edge midpoints {(x-1/2,y-1/2), (x-1/2,y+1/2), (x+1/2,y-1/2),
      (x+1/2,y+1/2)}, respectively, obtained from the slope maps,
      and the corresponding edge displacement vectors {(-1,-1), (-1,+1), (+1,-1), (+1,+1)}.
      These differences have reliability weights {wmm[x,y], wmp[x,y],
      wpm[x,y], wpp[x,y]}, respectively, which are defined to be the
      weights of the corresponding skew mismatch terms.
      
      For each height variable {z[x,y]} we select from these eight mismatch
      terms a subset of at most four terms. If an axial term has non-zero 
      weight, it is included in the system; otherwise one of the adjacent
      skew terms is included.      
      
      We want to minimize the sum {Q(z)} of all selected terms
      {qrs[x,y]}, each weighted by the corresponding weight {wrs[x,y]};
      where {rs} is {po}, {mo}, etc..
      
      This is a classical least squares problem. Minimizing {Q} is
      equivalent to finding the {z} variables that satisfy the
      equilibrium equations {dQ(z)/d(z[x,y]) = 0} for all {x,y}. Each
      term of {Q} contributes two terms for these equations. For
      instance, the term {wmm[x,y]*qmm[x,y]} contributes
      
        { d(wmm[x,y]*qmm[x,y])/d(z[x,y]) = 2*wmm[x,y]*(dmm[x,y] + z[x,y]) } 
        
      to the equation {x,y}, and 
      
        { d(wmm[x,y]*qmm[x,y])/d(z[x-1,y-1]) = 2*wmm[x,y]*(dmm[x,y] - z[x-1,y-1]) }
      
      to the equation {x-1,y-1}.
    */

    auto void append_equation_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy);
      /* Tries to add to equation {eqk} a term corresponding to a grid
        edge or diagonal that starts at {(x,y)} and ends at
        {(x+ux,y+uy)}, where {ux,uy} are in {-1..+1} (but not both
        zero). If the weight of that edge turns out to be zero, does
        nothing. Otherwise appends the term for the height value
        {z[x+ux,y+uy]}, increments term 0 (expected to be the height value
        {z[x,y]}), increments {eq->rhs,eq->wtot}. */

    /* Make sure that we can build the system: */
    assert(MAX_COEFFS >= 5);
    
    uint32_t N = S->N;
    
    uint32_t N_hors = 0; /* Equations that do not correspond to any height map pixel. */
    uint32_t N_good = 0; /* Equations that got at least one term with nonzero weight. */
    uint32_t N_null = 0; /* Equations that correspond to pixels but are in a zero-weight hole. */
    
    for (int32_t k = 0; k < N; k++)
      { uint32_t uidk = (uint32_t)k;  /* Index of main unknow of equation {k}. */
        pst_imgsys_equation_t *eqk = &(S->eq[k]);
        
        /* Initialize equation: */
        eqk->nt = 1;
        eqk->uid[0] = uidk;
        eqk->cf[0] = 0;
        eqk->wtot = 0;
        eqk->rhs = 0;
        
        /* Get the indices of height map pixel correspoinding to {uidk}: */
        int32_t x = S->col[uidk];
        int32_t y = S->row[uidk];
        assert((x == -1) == (y == -1));
        if (x == -1)
          { /* Height Value {uidk} is not a height map pixel: */
            N_hors++;
          }
        else
          { /* The height value is pixel {x,y} of the height map. Try to add the axial terms: */
            append_equation_term(eqk, x, y, -1, 00);
            append_equation_term(eqk, x, y, +1, 00);
            append_equation_term(eqk, x, y, 00, -1);
            append_equation_term(eqk, x, y, 00, +1);
            /* If some of them were zero, try to add the diagonal terms: */
            if (eqk->nt < MAX_COEFFS) { append_equation_term(eqk, x, y, -1, -1); }
            if (eqk->nt < MAX_COEFFS) { append_equation_term(eqk, x, y, -1, +1); }
            if (eqk->nt < MAX_COEFFS) { append_equation_term(eqk, x, y, +1, -1); }
            if (eqk->nt < MAX_COEFFS) { append_equation_term(eqk, x, y, +1, +1); }
            if (eqk->nt == 1)
              { /* No equations for this height value: */ N_null++; }
            else
              { N_good++; }
          }
      }
    assert(N_hors + N_good + N_null == N);
    if (verbose)
      { fprintf(stderr, "%5d system variables are not height map pixels\n", N_hors);
        fprintf(stderr, "%5d height map pixels with null equations\n", N_null);
        fprintf(stderr, "%5d non-null equations found\n", N_good);
      }

    return;
       
    void append_equation_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy)
      { 
        /* Check if we got enough equations: */
        if (eqk->nt >= MAX_COEFFS) { return; }
        
        /* Check if both height values {Z[x,y]} and {Z[x',y']} are in the height map: */
        if ((x < 0) || (x >= S->NX)) { return; }
        if ((y < 0) || (y >= S->NY)) { return; }
        
        int32_t x1 = x + ux, y1 = y + uy;
        if ((x1 < 0) || (x1 >= S->NX)) { return; }
        if ((y1 < 0) || (y1 >= S->NY)) { return; }
        
        /* Check if {x',y'} is in the system: */
        int32_t uid1 = x1 + y1*S->NX;
        assert((uid1 >= 0) && (uid1 < N));
      
        /* Get the estimated height difference {d} and its weight {w}: */
        double d, w;
        pst_slope_map_get_edge_data(G, W, x, y, ux, uy, &d, &w);
        if (fabs(w) >= FLUFF)
          { /* Add term to equation */
            uint32_t nt = eqk->nt;
            eqk->uid[nt] = (uint32_t)uid1; 
            eqk->cf[nt] = -w; 
            eqk->cf[0] += w;
            eqk->rhs += -w*d; 
            eqk->wtot += w; 
            nt++;
            eqk->nt = nt;
          }
      }
  }
   
void pst_integrate_fill_holes(pst_imgsys_t *S, int32_t NX_Z, int32_t NY_Z)
  { uint32_t N = S->N; 

    auto void add_term(pst_imgsys_equation_t *eq, int32_t x1, int32_t y1);
      /* If the variable {x1,y1} exists, adds to {eq} a term with coefficient 1
        that pulls the variable of {eq} towards the mean of its neighbors. */
    
    for (uint32_t k = 0; k < N; k++)
      { uint32_t uidk = k;
        pst_imgsys_equation_t *eqk = &(S->eq[k]);
        if (pst_imgsys_equation_is_null(uidk, eqk, N))
          { assert(eqk->nt == 1);
            assert(eqk->uid[0] == uidk);
            /* Get neighbors if variable {uidk}: */
            int32_t x = S->col[uidk];
            int32_t y = S->row[uidk];
            assert((x >= 0) && (x < NX_Z));
            assert((y >= 0) && (y < NY_Z));
            add_term(eqk, x-1,y);
            add_term(eqk, x+1,y);
            add_term(eqk, x,y-1);
            add_term(eqk, x,y+1);
          }
        eqk->wtot = 1.0e-5;
      }      
    return;
    
    void add_term(pst_imgsys_equation_t *eq, int32_t x1, int32_t y1)
      { if ((x1 < 0) || (x1 >= NX_Z)) { return; }
        if ((y1 < 0) || (y1 >= NY_Z)) { return; }
        /* Variable {x1,y1} exists: */
        int32_t uid1 = x1 + y1*NX_Z;
        assert((uid1 >= 0) && (uid1 < S->N));
        assert((uid1 >= 0) && (uid1 < S->N));
        assert((eq->nt >= 1) && (eq->nt < MAX_COEFFS));
        int32_t j = (int32_t)eq->nt;
        eq->uid[j] = (uint32_t)uid1;
        eq->cf[j] = -1.0;
        eq->cf[0] += 1.0;
        (eq->nt)++;
      }
  }

