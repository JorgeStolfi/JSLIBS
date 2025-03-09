/* See pst_integrate.h */
/* Last edited on 2025-03-05 15:02:24 by stolfi */

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

   
#define MAX_COEFFS pst_imgsys_MAX_COEFFS

#define FLUFF (1.0e-140)
  /* Practically zero. */

/* IMPLEMENTATIONS */

pst_imgsys_t* pst_integrate_build_system
  ( float_image_t *G,
    float_image_t *H,
    double hintsWeight,
    int32_t indent,
    bool_t verbose
  ) 
  {
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand((NC_G == 2) || (NC_G == 3), "gradient map must have 2 or 3 channels");
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
   
    int32_t NC_H = 0;
    if (H != NULL)
      { /* Get the size of the system: */
        float_image_get_size(H, &NC_H, NULL, NULL);
        demand((NC_H == 1) || (NC_H == 2), "gradient map must have 1 or 2 channels");
        float_image_check_size(H, 0, NX_Z, NY_Z, "hints map {H} has wrong size");
      }
    
    pst_imgsys_t *S = pst_imgsys_new_grid(NX_Z, NY_Z);
    assert(S->N == NX_Z*NY_Z);
    
    /* Conceptually, we define four axial quadratic mismatch terms
      {qpo[x,y](Z), qmo[x,y](Z), qop[x,y](Z), qom[x,y](Z)} for each
      corner of the pixel grid, that depend on the height map {Z}.
      Each term is the square of the difference between two 
      estimates of the height increment along one of the edges incident to {(x,y)}:
      one obtained by interpolation of the slope map {G}, and the 
      other obtained by numerical difference of the (unknown) heights at the two
      endpoints of that edge. 
      
        { qpo[x,y](Z) = (dpo[x,y] - (-Z[x,y]  + Z[x+1,y]))^2 }

        { qmo[x,y](Z) = (dmo[x,y] - (-Z[x,y]  + Z[x-1,y]))^2 }

        { qop[x,y](Z) = (dop[x,y] - (-Z[x,y]  + Z[x,y+1]))^2 }

        { qom[x,y](Z) = (dom[x,y] - (-Z[x,y]  + Z[x,y-1]))^2 }

      where {dpo[x,y], dmo[x,y], dop[x,y], dom[x,y]} are the
      interpolated slopes at the edge midpoints {(x+1/2,y), (x-1/2,y),
      (x,y+1/2), (x,y-1/2)}, respectively. (We omit the channel index 0
      in {Z[0,x,y]} for clarity.) These slopes have interpolated
      reliability weights {wpo[x,y], wmo[x,y], wop[x,y], wom[x,y]},
      respectively, which are defined to be the weights of the
      corresponding axial mismatch terms.
      
      We also define four skew mismatch terms
      {qmm[x,y](Z), qmp[x,y](Z), qpm[x,y](Z), qpp[x,y](Z)} 
      
        { qmm[x,y](Z) = (dmm[x,y] - (-Z[x,y]  + Z[x-1,y-1]))^2 }
      
        { qmp[x,y](Z) = (dmp[x,y] - (-Z[x,y]  + Z[x-1,y+1]))^2 }
      
        { qpm[x,y](Z) = (dpm[x,y] - (-Z[x,y]  + Z[x+1,y-1]))^2 }
      
        { qpp[x,y](Z) = (dpp[x,y] - (-Z[x,y]  + Z[x+1,y+1]))^2 }
      
      where {dmm[x,y], dmp[x,y], dpm[x,y], dpp[x,y]} are the dot product
      of the gradient interpolated at the skew edge midpoints
      {(x-1/2,y-1/2), (x-1/2,y+1/2), (x+1/2,y-1/2), (x+1/2,y+1/2)},
      respectively, obtained from the slope maps, and the corresponding
      edge displacement vectors {(-1,-1), (-1,+1), (+1,-1), (+1,+1)}.
      These differences have reliability weights {wmm[x,y], wmp[x,y],
      wpm[x,y], wpp[x,y]}, respectively, which are defined to be the
      weights of the corresponding skew mismatch terms.
      
      For each height variable {Z[x,y]} we select from these eight mismatch
      terms a subset of at most four terms. If an axial term has non-zero 
      weight, it is included in the system; otherwise one of the adjacent
      skew terms is included.   
      
      If {H} is not null, we also add a central mismatch term 
      
        { qoo[x,y](Z) = (H[0,x,y] - Z[x,y])^2 }
        
      with weight {woo[x,y] = hintsWeight*H[1,x,y]}. However, if {H} is
      null or the weight {hintsWeight*H[1,x,y]} is zero, we add a term 
      that corresponds to assuming {H[0,x,y]} is zero with
      a very small positive weight.
      
      We want to minimize the sum {Q(Z)} of all selected terms
      {qrs[x,y]}, each weighted by the corresponding weight {wrs[x,y]};
      where {rs} is {po}, {mo}, etc..
      
      This is a classical least squares problem. Minimizing {Q} is
      equivalent to finding the {Z} variables that satisfy the
      equilibrium equations {dQ(Z)/d(Z[x,y]) = 0} for all {x,y}. Each
      term of {Q} contributes two terms for these equations. For
      instance, the term {wmm[x,y]*qmm[x,y]} contributes
      
        { d(wmm[x,y]*qmm[x,y])/d(Z[x,y]) = 2*wmm[x,y]*(dmm[x,y] + Z[x,y]) } 
        
      to the equation {x,y}, and 
      
        { d(wmm[x,y]*qmm[x,y])/d(Z[x-1,y-1]) = 2*wmm[x,y]*(dmm[x,y] - Z[x-1,y-1]) }
      
      to the equation {x-1,y-1}.  The term {woo[x,y]*doo[x,y]} contributes
      
        { d(woo[x,y]*qoo[x,y])/d(Z[x,y]) = 2*woo[x,y]*(H[0,x,y] - Z[x,y]) }
        
      to that same equation.
      
      The equation for {Z[x,y]} is obtained by adding all these partial equations
      with the indicated weights, leaving out the common factor '2*'. */

    auto bool_t append_edge_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy);
      /* Tries to add to equation {eqk} a term corresponding to a grid
        edge or diagonal that starts at {(x,y)} and ends at
        {(x+ux,y+uy)}, where {ux,uy} are in {-1..+1} (but not both
        zero). If the weight of that edge turns out to be zero, does
        nothing. Otherwise appends the term for the height value
        {Z[x+ux,y+uy]}, increments term 0 (expected to be the height value
        {Z[x,y]}), increments {eq->rhs,eq->wtot}.
        
        Returns {TRUE} iff the edge term was succesfully added. */

    auto void append_hints_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y);
      /* Tries to add to equation {eqk} a term corresponding to the
        equation {Z[0,x,y] == H[0,x,y]} with weight {hintsWeight*H[1,x,y]}. If
        {hintsWeight*H[1,x,y]} is zero, pretends that {H[0,x,y] = 0} and {hintsWeight*H[1,x,y]}
        is a very small weight. Increments the coefficient of term 0 of
        the equation, as well as {eq->rhs,eq->wtot}. */

    auto void append_fudge_term(pst_imgsys_equation_t *eqk);
      /* Adds to to equation {eqk} a term corresponding to the
        equation {Z[0,x,y] == 0} with weight 1.0.  Must be called only if {eqk.wtot}
        and {eqk.nt} are zero. */

    /* Make sure that we can build the system: */
    assert(MAX_COEFFS >= 5);
    
    bool_t ignoreSlopes = FALSE; /* Normally {FALSE}; {TRUE} for debugging. */
    
    uint32_t N = S->N;
    
    uint32_t N_hors = 0; /* Equations that do not correspond to any height map pixel. */
    uint32_t N_cook = 0; /* Equations that ended up {0==0} and had to be cooked. */
    
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
          { /* The height value is pixel {x,y} of the height map. */
            if (! ignoreSlopes)
              { /* Try to add the axial terms: */
                bool_t axm = append_edge_term(eqk, x, y, -1, 00);
                bool_t axp = append_edge_term(eqk, x, y, +1, 00);
                bool_t aym = append_edge_term(eqk, x, y, 00, -1);
                bool_t ayp = append_edge_term(eqk, x, y, 00, +1);
                /* If two adjacent , try to add the diagonal terms: */
                if ((! axm) || (! aym)) { (void)append_edge_term(eqk, x, y, -1, -1); }
                if ((! axm) || (! ayp)) { (void)append_edge_term(eqk, x, y, -1, +1); }
                if ((! axp) || (! aym)) { (void)append_edge_term(eqk, x, y, +1, -1); }
                if ((! axp) || (! ayp)) { (void)append_edge_term(eqk, x, y, +1, +1); }
              }
            /* This must come last: */
            append_hints_term(eqk, x, y);
            /* Check final equation weight: */
            if (eqk->wtot == 0)
              { /* Add a fudge term pulling the variable to 0: */
                append_fudge_term(eqk);
                N_cook++;
              }
          }
      }

    if (verbose)
      { fprintf(stderr, "%*s%5d system variables are not height map pixels\n", indent, "", N_hors);
        fprintf(stderr, "%*s%5d null equations cooked to pull to zero\n", indent, "", N_cook);
      }
      
    /* Normalize {wtot} fields: */
    double nt_exp = (ignoreSlopes ? 0 : 4) + 1.0;
    for (int32_t k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(S->eq[k]);
        eqk->wtot /= nt_exp;
      }

    return S;
       
    bool_t append_edge_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy)
      { 
        bool_t debug = FALSE;
        
        /* Check if we got enough equations: */
        if (eqk->nt >= MAX_COEFFS) { return FALSE; }
      
        /* Term 0 must be there, even if with zero coef: */
        assert(eqk->nt >= 1);
        
        /* Check if both height values {Z[x,y]} and {Z[x',y']} are in the height map: */
        if ((x < 0) || (x >= S->NX)) { return FALSE; }
        if ((y < 0) || (y >= S->NY)) { return FALSE; }
        
        int32_t x1 = x + ux, y1 = y + uy;
        if ((x1 < 0) || (x1 >= S->NX)) { return FALSE; }
        if ((y1 < 0) || (y1 >= S->NY)) { return FALSE; }
        
        /* Check if {x',y'} is in the system: */
        int32_t uid1 = x1 + y1*S->NX;
        assert((uid1 >= 0) && (uid1 < N));
      
        /* Get the estimated height difference {vD} and its weight {wD}: */
        double vD, wD;
        pst_slope_map_get_edge_data(G, x, y, ux, uy, &vD, &wD);
        if (debug && ((x == 0) || (x == S->NX-1) || (y == 0) || (y == S->NY-1)))
          { fprintf(stderr, "%*s  @ {append_edge_term} [%d,%d] - [%d,%d] %+10.6f %8.4f\n", indent, "", x, y, x+ux, y+uy, vD, wD); }
        assert(isfinite(wD) && (wD >= 0));
        if (fabs(wD) >= FLUFF)
          { /* Add term to equation */
            assert(isfinite(vD));
            uint32_t nt = eqk->nt;
            eqk->uid[nt] = (uint32_t)uid1; 
            eqk->cf[nt] = -wD; 
            eqk->cf[0] += wD;
            eqk->rhs += -wD*vD; 
            eqk->wtot += wD; 
            nt++;
            eqk->nt = nt;
            return TRUE;
          }
        else
          { return FALSE; }
      }
       
    void append_hints_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y)
      { 
        /* Check if height values {Z[x,y]} is in the height map: */
        if ((x < 0) || (x >= S->NX)) { return; }
        
        /* Check if {x,y} is in the system: */
        int32_t uid = x + y*S->NX;
        assert((uid >= 0) && (uid < N));
      
        /* Term 0 must be there, even if with zero coef: */
        assert(eqk->nt >= 1);
        assert(eqk->uid[0] == uid);
        
        /* Get the hint (external estimated height) {vH} and its weight {wH}: */
        double wH = ((H == NULL) ? 0.0 : (NC_H < 2 ? 1.0 : float_image_get_sample(H, 1, x, y)));
        demand(isfinite(wH) && (wH >= 0), "invalid hint weight");
        double vH = NAN;  /* Hint value, or 0 if undefined. */
        if (fabs(wH) >= FLUFF)
          { vH = ((H == NULL) || (wH == 0) ? 0.0 : float_image_get_sample(H, 0, x, y));
            if (! isfinite(vH)) { wH = 0; vH = NAN; }
          }
        /* Adjust the weight: */
        wH = wH * hintsWeight;
        if (! isfinite(wH)) { wH = 0.0; }
        if (fabs(wH) >= FLUFF)
          { /* Add term {Z=H} to equation */
            eqk->cf[0] += wH;
            eqk->rhs += wH*vH; 
            eqk->wtot += wH; 
          }
      }
      
    void append_fudge_term(pst_imgsys_equation_t *eqk)
      {
        /* Term 0 must be there, but with zero coef: */
        assert(eqk->nt == 1);
        assert(eqk->cf[0] == 0);
        assert(eqk->wtot == 0);
        assert(eqk->rhs == 0);
        
        double vC = 0.0;  /* Default value. */
        double wC = 1.0;  /* Irrelevant. */
        eqk->cf[0] += wC;
        eqk->rhs += vC; 
        eqk->wtot += wC; 
      }
  }
