/* See pst_slope_map.h */
/* Last edited on 2024-12-23 07:08:22 by stolfi */

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
#include <pst_interpolate.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_slope_map.h>

/* INTERNAL PROTOTYPES */

void pst_slope_map_extract_system_weight_image(pst_imgsys_t *S, float_image_t *OW);
  /* Stores into {OW} the equation weights {.wtot} of the system {S}. */ 

pst_imgsys_equation_t *pst_slope_map_get_eqs(float_image_t *G, float_image_t *GW, bool_t fillHoles);
  /* Builds the equation list {eq} for the integration system, given
    the gradient map {G} and its weight map {GW}. If {fillHoles} is
    false, the equations for isolated vertices will have zero {wtot}
    and will set the height to zero. If {fillHoles} is true, those
    equations will set the heights to the average of their neighbors,
    with an arbitrary nonzero {wtot}. */

/* IMPLEMENTATIONS */

r2_t pst_slope_map_get_pixel(float_image_t *G, int32_t x, int32_t y)
  { demand(G->sz[0] == 2, "wrong slope map depth");
    r2_t grd;
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { grd.c[axis] = (*p); p += G->st[0]; }
    return grd;
  }

void pst_slope_map_set_pixel(float_image_t *G, int32_t x, int32_t y, r2_t *grd)
  { demand(G->sz[0] == 2, "wrong slope map depth");
    float *p = float_image_get_sample_address(G, 0, x, y);
    for (uint32_t axis = 0; axis < 2; axis++)
      { (*p) = (float)(grd->c[axis]); p += G->st[0]; }
  }

void pst_slope_map_to_depth_map_recursive
  ( float_image_t *IG, 
    float_image_t *IW, 
    uint32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportIter,
    pst_slope_map_report_proc_t *reportData,
    pst_imgsys_report_proc_t *reportSys,
    pst_height_map_report_proc_t *reportHeights
  )
  {
    /* Get image size and check compatibility: */
    demand(IG->sz[0] == 2, "wrong {IG} channels");
    uint32_t NX_G = (uint32_t)IG->sz[1];
    uint32_t NY_G = (uint32_t)IG->sz[2];
    
    if (IW != NULL)
      { demand(IW->sz[0] == 1, "wrong {IW} channels");
        demand(IW->sz[1] == NX_G, "wrong {IW} cols");
        demand(IW->sz[2] == NY_G, "wrong {IW} rows");
      }
      
    /* Allocate output images: */
    uint32_t NX_Z = NX_G + 1;
    uint32_t NY_Z = NY_G + 1;
    float_image_t *OZ = float_image_new(1, (int32_t)NX_Z, (int32_t)NY_Z);  
    float_image_t *OW = ((IW == NULL) || (OWP == NULL) ? NULL : float_image_new(1, (int32_t)NX_Z, (int32_t)NY_Z));
    if (OWP != NULL) { (*OWP) = OW; }

    uint32_t indent = 2*level+2; /* Indentation for messages. */

    if (verbose)
      { fprintf(stderr, "%*sEntering level %d with G size %d×%d ...\n", indent, "", level, NX_G, NY_G); }

    if (reportData != NULL) { reportData(level, IG, IW); }

    /* Decide whether to use multi-scale: */
    bool_t trivial = ((NX_G == 1) && (NY_G == 1));
    
    /* Get the initial guess {OZ} for the heights: */
    if (trivial)
      { /* Solution for the 1x1 system: */
        /* Get slopes: */
        float dZdX = float_image_get_sample(IG, 0, 0, 0); 
        float dZdY = float_image_get_sample(IG, 1, 0, 0); 
        /* Assume an affine function that is 0 at center, with the given gradient: */
        double z10 = 0.5*(dZdX - dZdY);
        double z11 = 0.5*(dZdX + dZdY);
        float_image_set_sample(OZ, 0, 0, 0, (float)-z11);
        float_image_set_sample(OZ, 0, 0, 1, (float)-z10);
        float_image_set_sample(OZ, 0, 1, 0, (float)+z10);
        float_image_set_sample(OZ, 0, 1, 1, (float)+z11);
        if (OW != NULL) { float_image_fill_channel(OW, 0, 1.0); }
      }
    else
      { /* Initialize {OZ} with the solution to a coarser problem. */
      
        /* Shrink slope maps and weight map by half: */
        if (verbose) { fprintf(stderr, "%*sShrinking slope and weight maps ...\n", indent, ""); }
        float_image_t *SIG; 
        float_image_t *SIW;
        pst_slope_and_weight_map_shrink(IG,IW, &SIG,&SIW);
        
        /* Compute the half-scale height map: */
        float_image_t *SOZ, *SOW;
        pst_slope_map_to_depth_map_recursive
          ( SIG, SIW, level+1, 
            2*maxIter, convTol/2, topoSort,
            &SOZ, &SOW,
            verbose, reportIter,
            reportData, reportSys, reportHeights
          );
        
        /* Expand the computed height map to double size: */
        if (verbose) { fprintf(stderr, "%*sExpanding height map to %d×%d ...\n", indent, "", NX_Z, NY_Z); }
        pst_height_map_expand(SOZ, SOW, OZ, OW);
          
        /* Free the working storage: */
        float_image_free(SIG);
        if (SIW != NULL) { float_image_free(SIW); }
        float_image_free(SOZ);
        if (SOW != NULL) { float_image_free(SOW); }
      }
    (*OZP) = OZ;

    bool_t solve_sys = ! trivial;
    if (solve_sys)
      { /* Build the linear system: */
        uint32_t NP_Z = NX_Z*NY_Z;
        if (verbose) { fprintf(stderr, "%*sBuilding linear system for %d pixels ...\n", indent, "", NP_Z); }
        bool_t full = FALSE; /* Should be parameter. FALSE means exclude indeterminate pixels. */
        pst_imgsys_t *S = pst_slope_map_build_integration_system(IG, IW, full);
        if (OW != NULL) { pst_slope_map_extract_system_weight_image(S, OW); }
        if (verbose) { fprintf(stderr, "%*sSystem has %d equations and %d unknowns\n", indent, "", S->N, S->N); }
        if (reportSys != NULL) { reportSys(level, S); }

        /* Solve the system for the corner heights: */
        if (verbose) { fprintf(stderr, "%*sSolving the system ...\n", indent, ""); }
        bool_t para = FALSE; /* Should be parameter. 1 means parallel execution, 0 sequential. */
        bool_t szero = TRUE; /* Should be parameter. 1 means adjust sum to zero, 0 let it float. */
        uint32_t *ord = NULL;
        if (topoSort) { ord = pst_imgsys_sort_equations(S); }
        pst_slope_map_solve_system(S, OZ, ord, maxIter, convTol, para, szero, verbose, level, reportIter, reportHeights);
        pst_imgsys_free(S);
        if (ord != NULL) { free(ord); }
      }
    else
      { /* The initial guess is the solution: */
        if (reportSys != NULL) { reportSys(level, NULL); }
        if (reportHeights != NULL) { reportHeights(level, 0, 0.0, TRUE, OZ); }
      }

    if (verbose)
      { fprintf
          ( stderr, ("%*sLeaving level %d (OZ = %" int64_d_fmt "×%" int64_d_fmt ") ...\n"), 
            indent, "", level, OZ->sz[1], OZ->sz[2]
          );
      }
  }

void pst_slope_and_weight_map_shrink
  ( float_image_t *IG, 
    float_image_t *IW, 
    float_image_t **SG, 
    float_image_t **SW
  )
  { 
    uint32_t NC = (uint32_t)IG->sz[0];
    uint32_t NXI = (uint32_t)IG->sz[1]; uint32_t NXJ = (NXI+1)/2;
    uint32_t NYI = (uint32_t)IG->sz[2]; uint32_t NYJ = (NYI+1)/2;
    float_image_t *JG = float_image_new((int32_t)NC, (int32_t)NXJ, (int32_t)NYJ);
    float_image_t *OW = NULL;
    
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert(IW->sz[1] == NXI);
        assert(IW->sz[2] == NYI);
	OW = float_image_new(1, (int32_t)NXJ, (int32_t)NYJ);
      }

    for (int32_t jy = 0; jy < NYJ; jy++)
      { for (int32_t jx = 0; jx < NXJ; jx++)
          { int32_t ix = 2*jx, iy = 2*jy;
	    float wmin = INF;
	    for (uint32_t c = 0; c < 2; c++)
	    {
	      double da,wa;
	      pst_interpolate_two_samples(IG, IW, (int32_t)c, ix,iy,  ix+1,iy+1, &da,&wa);
	      double db,wb;
 	      pst_interpolate_two_samples(IG, IW, (int32_t)c, ix,iy+1, ix+1,iy,  &db,&wb);
	      
	      float d = (float)( wa+wb > 0 ? (wa*da + wb*db)/(wa + wb) : (da+db)/2.0);
	      float_image_set_sample(JG, (int32_t)c, jx, jy, d);
	      float w = (float)(4/(1/wa + 1/wb));
	      if(w < wmin) wmin = w;
	    }
   
	    if (OW != NULL) { float_image_set_sample(OW, 0, jx, jy, wmin);}
          
          }
      }
    *SG = JG;
    *SW = OW;
  }

void pst_slope_map_solve_system
  ( pst_imgsys_t *S, 
    float_image_t *OZ,
    uint32_t ord[], 
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose, 
    uint32_t level, 
    uint32_t reportIter, 
    pst_height_map_report_proc_t *reportHeights
  )
  {
    uint32_t indent = 2*level+2;
    
    auto void reportSol(uint32_t iter, double change, bool_t final, uint32_t N, double Z[]);
      /* A procedure that is called by the Gauss-Seidel solver at each iteration.
         When appropriate, it copies {Z} into the image {OZ} and calls {reportHeights}. 
         This happens when {final} is TRUE or when {reportIter != 0} and {iter}
         is a multiple of {reportIter}. */
    
    void reportSol(uint32_t iter, double change, bool_t final, uint32_t N, double Z[])
      { bool_t doit = final || ((reportHeights != NULL) && (reportIter != 0) && (iter % reportIter == 0));
        if (doit) 
          { pst_slope_map_copy_sol_vec_to_height_map(S, Z, OZ);
            if (reportHeights != NULL) { reportHeights(level, iter, change, final, OZ); }
          }
      }

    double *Z = rn_alloc(S->N);
    pst_slope_map_copy_height_map_to_sol_vec(S, OZ, Z);
    pst_imgsys_solve(S, Z, ord, maxIter, convTol, para, szero, verbose, indent, reportSol);
    /* {reportSol} must have copied the final solution into {OZ}. */
    free(Z);
  }

void pst_slope_map_copy_height_map_to_sol_vec(pst_imgsys_t *S, float_image_t *IZ, double VZ[])
  {
    uint32_t N = S->N;
    for (uint32_t k = 0; k < N; k++)
      { uint32_t x = S->col[k];
        uint32_t y = S->row[k];
        VZ[k] = float_image_get_sample(IZ, 0, (int32_t)x, (int32_t)y);
      }
  }

void pst_slope_map_copy_sol_vec_to_height_map(pst_imgsys_t *S, double VZ[], float_image_t *IZ)
  {
    uint32_t NX = (uint32_t)IZ->sz[1];
    uint32_t NY = (uint32_t)IZ->sz[2];
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { int32_t k = S->ix[x + (int32_t)NX*y];
            float_image_set_sample(IZ, 0, x, y, (float)(k < 0 ? 0.0 : VZ[k]));
          }
      }
  }

pst_imgsys_t* pst_slope_map_build_integration_system
  ( float_image_t *G, 
    float_image_t *GW, 
    bool_t full
  ) 
  {
    /* Get/check the sizes of the slope maps: */
    uint32_t NX_G = (uint32_t)G->sz[1];
    uint32_t NY_G = (uint32_t)G->sz[2];
    
    /* Get the size of the system: */
    uint32_t NX_Z = NX_G+1;
    uint32_t NY_Z = NY_G+1;
    uint32_t NXY_Z = NX_Z*NY_Z; /* Number of unknowns (heights) in full system. */
    
    /* Gather the equations {eq[0..NXY_Z-1]} of the system, using the full variable numbering: */
    bool_t fillHoles = full; /* TRUE to set isolated pixels to Poisson interpolant. */
    pst_imgsys_equation_t *eq = pst_slope_map_get_eqs(G, GW, fillHoles);
    
    /* Build the index {ix[0..N-1]} and inverse indices {col[0..N-1],row[0..N-1]}: */
    int32_t *ix = talloc(NXY_Z, int32_t);
    uint32_t *col = talloc(NXY_Z, uint32_t);
    uint32_t *row = talloc(NXY_Z, uint32_t);
    for (uint32_t xy = 0; xy < NXY_Z; xy++) { ix[xy] = (int32_t)xy; col[xy] = xy % NX_Z; row[xy] = xy / NX_Z; }
     
    uint32_t N; /* Number of valid equations. */
    if (full)
      { /* Keep all equations. */
        N = NXY_Z;
      }
    else
      { /* Compress the valid equations {eq[0..N-1]} of the system, and fill the table {ix[0..NXY_Z-1]}: */
        N = 0;
        for (uint32_t xy = 0; xy < NXY_Z; xy++)
          { if (eq[xy].wtot > 0.0)
              { /* Keep equation: */
                assert(N < NXY_Z);
                eq[N] = eq[xy]; col[N] = col[xy]; row[N] = row[xy];
                ix[xy] = (int32_t)N;
                N++;
              }
            else
              { /* Discard equation: */
                ix[xy] = -1;
              }
          }

        /* Replace the temporary indices in the equations by the correct ones, eliminating the excluded vars: */
        for (uint32_t k = 0; k < N; k++)
          { pst_imgsys_equation_t *eqk = &(eq[k]);
            assert(eqk->ix[0] == k);
            uint32_t nt = eqk->nt;
            uint32_t mt = 0;
            for (uint32_t i = 0; i < nt; i++)
              { /* Get the temporay index {xyi}: */
                uint32_t xyi = eqk->ix[i];
                /* Get the definitive index {ki}: */
                int32_t ki = ix[xyi];
                if (ki >= 0)
                  { /* Append the term to the equation: */
                    uint32_t j = mt;
                    eqk->ix[j] = (uint32_t)ki;
                    eqk->cf[j] = eqk->cf[i];
                    assert(!isnan(eqk->cf[j]));
                    mt++;
                  }
              }
            eqk->nt = mt;
          }
      }
     
    /* Now package the equations as a system: */
    pst_imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, ix, col, row);  
      
    return S;
  }

pst_imgsys_equation_t *pst_slope_map_get_eqs(float_image_t *G, float_image_t *GW, bool_t fillHoles)
  {
    /* Get/check the sizes of the slope maps: */
    uint32_t NC = (uint32_t)G->sz[0]; assert(NC == 2);
    uint32_t NX_G = (uint32_t)G->sz[1];
    uint32_t NY_G = (uint32_t)G->sz[2];
    
    if (GW != NULL)
      { demand(GW->sz[0] == 1, "wrong {GW} channels");
        demand(GW->sz[1] == NX_G, "wrong {GW} cols");
        demand(GW->sz[2] == NY_G, "wrong {GW} rows");
      }
    
    /* Check the size of the system: */
    uint32_t NX_Z = NX_G+1;
    uint32_t NY_Z = NY_G+1;
    uint32_t NXY_Z = NX_Z*NY_Z; /* Number of unknowns (heights). */
    
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
      
      For each unknown height {z[x,y]} we select from these eight mismatch
      terms a subset of at most four terms. If an axial term has non-zero 
      weight, it is included in the system; otherwise one of the adjacent
      skew terms is included.      
      
      We want to minimize the sum {Q(h)} of all selected terms,
      each weighted by the corresponding weight. Minimizing {Q}
      is equivalent to finding the {z} variables that
      satisfy the equilibrium equations.
      
      Note that if the slope map {G} is zero, we can anihilate all
      those terms (and therefore Q) by setting all heights equal to 1.
      Therefore, if the system is to be solved iteratively, one must
      take care to exclude this solution.
    */
    
    auto void compute_axial_edge_delta_and_weight(uint32_t axis, int32_t x, int32_t y, int32_t u, double *dP, double *wP);
      /* Computes the delta {*dP} and the weight {*wP} of the axial
        edge that goes from {(x,y)} to {(x+u,y)} if {axis=0} or to
        {(x,y+u)} if {axis=1}, where {u} is either {+1} or {-1}.
        Assumes that {x,y} is a vertex of the cell grid. */

    void compute_axial_edge_delta_and_weight(uint32_t axis, int32_t x, int32_t y, int32_t u, double *dP, double *wP)
      {
        assert((u == -1) || (u == +1));
        (*dP) = 0.0; (*wP) = 0.0; /* If edge does not exist. */
        if (axis == 0)
          { /* X derivative */
            int32_t xg; /* Indices of relevant samples in {G,GW}. */
            if (u < 0) 
              { xg = x - 1; if (xg < 0) { return; } }
            else
              { xg = x; if (xg >= NX_G) { return; } }
            pst_interpolate_four_samples(G,GW, 0, xg, y-1, xg, y, dP,wP);
            (*dP) = (*dP) * u;
          } 
        else if (axis == 1)
          { /* Y derivative */
            int32_t yg; /* Indices of relevant samples in {G,GW}. */
            if (u < 0) 
              { yg = y - 1; if (yg < 0) { return; } }
            else
              { yg = y; if (yg >= NY_G) { return; } }
            pst_interpolate_four_samples(G,GW, 1, x-1,yg, x,yg, dP,wP);
            (*dP) = (*dP) * u;
          } 
        else
          { demand(FALSE, "invalid axis"); }
        return;
      }

    auto void compute_skew_edge_delta_and_weight(int32_t x, int32_t y, int32_t ux, int32_t uy, double *dP,  double *wP);
      /* Computes the delta {*dP} and the weight {*wP} of the skew edge that goes
        from {(x,y)} to {(x+ux,y+uy)}, where {ux} and {uy} are either {+1} or {-1}.
        Assumes that {x,y} is a vertex of the cell grid. */
    
    void compute_skew_edge_delta_and_weight(int32_t x, int32_t y, int32_t ux, int32_t uy, double *dP,  double *wP)
      {
        assert((ux == -1) || (ux == +1));
        assert((uy == -1) || (uy == +1));
        int32_t xg, yg; /* Indices of relevant samples in {G,GW}. */
        (*dP) = 0.0; (*wP) = 0.0; /* If edge does not exist. */
        if (ux < 0) 
          { xg = x - 1; if (xg < 0) { return; } }
        else
          { xg = x; if (xg >= NX_G) { return; } }
        if (uy < 0) 
          { yg = y - 1; if (yg < 0) { return; } }
        else
          { yg = y; if (yg >= NY_G) { return; } }
        double dx = float_image_get_sample(G, 0, xg, yg);
        double dy = float_image_get_sample(G, 1, xg, yg);
        double d = dx*ux + dy*uy;
        double w = (GW == NULL ? 1.0 : float_image_get_sample(GW, 0, xg, yg));
        assert((! isnan(w)) && (w >= 0));
        (*dP) = d;
        (*wP) = w;
      }

    auto void append_equation_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy, double d, double w);
      /* Appends to equation {eqk} a term corresponding to an edge with delta {d} and weight {w}
        that starts at {(x,y)} and ends at {(x+ux,y+uy)}. */
       
    void append_equation_term(pst_imgsys_equation_t *eqk, int32_t x, int32_t y, int32_t ux, int32_t uy, double d, double w)
      { uint32_t nt = eqk->nt;
        int32_t k = (x+ux) + (y+uy)*(int32_t)NX_Z;
        assert(k >= 0);
        eqk->ix[nt] = (uint32_t)k; 
        eqk->cf[nt] = -w; 
        eqk->rhs += -w*d; 
        eqk->wtot += w; 
        nt++;
        eqk->nt = nt;
      }

    /* Make sure that we can build the system: */
    assert(MAXCOEFS >= 5);
    
    pst_imgsys_equation_t *eq = (pst_imgsys_equation_t*)notnull(malloc(NXY_Z*sizeof(pst_imgsys_equation_t)), "no mem");

    for (uint32_t xy = 0; xy < NXY_Z; xy++)
      { /* Get the indices of pixel number {k}: */
        int32_t x = (int32_t)(xy % NX_Z);
        int32_t y = (int32_t)(xy / NX_Z);
        
        /* Get the axial deltas and their weights: */
        double wmo, wpo, wom, wop;  /* Edge weights. */
        double dmo, dpo, dom, dop;  /* Edge deltas. */
        compute_axial_edge_delta_and_weight(0, x,y, -1, &dmo,&wmo);
        compute_axial_edge_delta_and_weight(0, x,y, +1, &dpo,&wpo);
        compute_axial_edge_delta_and_weight(1, x,y, -1, &dom,&wom);
        compute_axial_edge_delta_and_weight(1, x,y, +1, &dop,&wop);
        assert((! isnan(wmo)) && (wmo >= 0));
        assert((! isnan(wpo)) && (wpo >= 0));
        assert((! isnan(wom)) && (wom >= 0));
        assert((! isnan(wop)) && (wop >= 0));

        /* Get the skew deltas and weights: */
        double wmm, wmp, wpm, wpp;  /* Edge weights. */
        double dmm, dmp, dpm, dpp;  /* Edge deltas. */
        compute_skew_edge_delta_and_weight(x,y, -1,-1, &dmm,&wmm);
        compute_skew_edge_delta_and_weight(x,y, -1,+1, &dmp,&wmp);
        compute_skew_edge_delta_and_weight(x,y, +1,-1, &dpm,&wpm);
        compute_skew_edge_delta_and_weight(x,y, +1,+1, &dpp,&wpp);
        assert((! isnan(wmm)) && (wmm >= 0));
        assert((! isnan(wmp)) && (wmp >= 0));
        assert((! isnan(wpm)) && (wpm >= 0));
        assert((! isnan(wpp)) && (wpp >= 0));
        
        /* Select which mismatch terms to include in the {Q} function: */
        bool_t use_mo = FALSE, use_po = FALSE, use_om = FALSE, use_op = FALSE;
        bool_t use_mm = FALSE, use_mp = FALSE, use_pm = FALSE, use_pp = FALSE;

        if (wmo > 0) { use_mo = TRUE; } else if ( wmm > 0) { use_mm = TRUE; } else if (wmp > 0) { use_mp = TRUE; }
        if (wpo > 0) { use_po = TRUE; } else if ( wpm > 0) { use_pm = TRUE; } else if (wpp > 0) { use_pp = TRUE; }
        if (wom > 0) { use_om = TRUE; } else if ( wmm > 0) { use_mm = TRUE; } else if (wpm > 0) { use_pm = TRUE; }
        if (wop > 0) { use_op = TRUE; } else if ( wmp > 0) { use_mp = TRUE; } else if (wpp > 0) { use_pp = TRUE; }
        
        uint32_t  n_use = 
          (uint32_t)use_mo + (uint32_t)use_po + (uint32_t)use_om + (uint32_t)use_op + 
          (uint32_t)use_mm + (uint32_t)use_mp + (uint32_t)use_pm + (uint32_t)use_pp;
        if (n_use == 0)
          { /* Vertex is isolated. */
            if (fillHoles)
              { /* Provide axial equations for existing edges. */
                double wtiny = 1.0e-50;
                wmo = (x <= 0    ? 0.0 : wtiny);  dmo = 0.0; use_mo = (wmo > 0);
                wpo = (x >= NX_G ? 0.0 : wtiny);  dpo = 0.0; use_po = (wpo > 0);
                wom = (y <= 0    ? 0.0 : wtiny);  dom = 0.0; use_om = (wom > 0);
                wop = (y >= NY_G ? 0.0 : wtiny);  dop = 0.0; use_op = (wop > 0);
                n_use = use_mo + use_po + use_om + use_op;
                assert(n_use > 0);
              }
            else
              { /* Leave it isolated: */
              }
          }
        
        /* Combine the selected mismatch terms in one normalized equilibrium equation: */
        /* Also compute the total vertex weight weight {wtot}: */
        pst_imgsys_equation_t *eqk = &(eq[xy]);

        /* Initialize equation: */
        eqk->nt = 0;     /* Number of terms in equation. */
        eqk->rhs = 0.0;  /* Right-hand side. */  
        eqk->wtot = 0.0; /* Weight of equation. */
        
        /* Set the diagonal term of equation (already normalized): */
        eqk->ix[0] = xy;  eqk->cf[0] = 1.0; eqk->nt = 1;

        if (n_use > 0)
          {
            assert((1+n_use) <= MAXCOEFS);
            
            /* Collect non-diagonal equation terms of selected mismtch terms with unnormalized coeffs: */
            if ( use_mo ) { append_equation_term(eqk, x,y, -1,00, dmo,wmo); }
            if ( use_po ) { append_equation_term(eqk, x,y, +1,00, dpo,wpo); }
            if ( use_om ) { append_equation_term(eqk, x,y, 00,-1, dom,wom); }
            if ( use_op ) { append_equation_term(eqk, x,y, 00,+1, dop,wop); }
            
            if ( use_mm ) { append_equation_term(eqk, x,y, -1,-1, dmm,wmm); }
            if ( use_mp ) { append_equation_term(eqk, x,y, -1,+1, dmp,wmp); }
            if ( use_pm ) { append_equation_term(eqk, x,y, +1,-1, dpm,wpm); }
            if ( use_pp ) { append_equation_term(eqk, x,y, +1,+1, dpp,wpp); }
             
            assert(eqk->nt == 1 + n_use);

            /* Compute {wtot} and normalize coeficients: */
            assert(eqk->wtot > 0);
            for (uint32_t j = 1; j < eqk->nt; j++) { eqk->cf[j] /= eqk->wtot;  }
            eqk->rhs /= eqk->wtot;
          }
      }
    return eq;
  }

void pst_slope_map_extract_system_weight_image(pst_imgsys_t *S, float_image_t *OW)
  { 
    demand(OW->sz[0] == 1, "wrong {OW} channels");
    int32_t NX = (int32_t)OW->sz[1];
    int32_t NY = (int32_t)OW->sz[2];
    int32_t NXY = NX*NY; /* Number of unknowns (heights). */
    demand(S->N <= NXY, "wrong {OW} size");
    
    for (int32_t xy = 0; xy < NXY; xy++)
      { /* Get the indices of pixel number {k}: */
        int32_t x = xy % (int32_t)NX;
        int32_t y = xy / (int32_t)NX;
        /* Get total mismatch weight of equation for vertex {(x,y)}: */
        int32_t k = S->ix[xy];
        assert(k >= 0);
        double wtot = S->eq[k].wtot;
        float_image_set_sample(OW, 0, x,y, (float)wtot);
      }
  }
