/* See pst_fit_light.h */
/* Last edited on 2025-01-19 06:47:06 by stolfi */ 

#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <float_image.h>
#include <r3.h> 
#include <rn.h> 
#include <rmxn.h> 
#include <gausol_solve.h> 
#include <qmin_simplex.h> 
#include <affirm.h> 

#include <pst_normal_map.h>
#include <pst_basic.h>
#include <pst_geom.h>
#include <pst_lamp.h>
#include <pst_light.h>
#include <pst_shading.h>

#include <pst_fit_light.h>

/* INTERNAL PROTOTYPES - LEAST-SQUARES SYSTEM BUILDERS AND SOLVERS */

typedef double pst_flt_basis_func_t(uint32_t i, r3_t *nrm); 
  /* An element of the function basis to use for least-squares
    fitting. The integer {i} is the index of the element in the basis.
    The domain is the set of all normal directions; the value
    is the apparent intensity of a surface with normal {nrm} that
    results from this term alone, with unit power. */
  
typedef double pst_flt_target_func_t(r3_t *nrm, double smp); 
  /* The target function to use for least-squares fitting. The domain
    is the set of all normal directions; the value is the desired
    apparent intensity at any pixel where the normal is {nrm} and the
    image value in the relevant channel is {smp}. */

typedef double pst_flt_weight_func_t(r3_t *nrm, double smp); 
  /* The weight function to use for least-squares fitting. The domain
    is the set of all normal directions; the value is the weight to be
    given to any pixel where the normal is {nrm} and the image
    value in the relevant channel is {smp}. */

void pst_flt_build_lsq_system
  ( float_image_t *IMG, 
    int32_t c,
    float_image_t *NRM, 
    uint32_t n,
    pst_flt_basis_func_t *bas[],
    pst_flt_target_func_t *fun,
    pst_flt_weight_func_t *wht,
    double A[],
    double b[]
  );
  /* Computes the {n�n} matrix {A} and the right-hand-side {n}-vector {b} of
    the the least squares system {A u = b} for a linear shading model,
    for channel {c} of the image {IMG}.
    
    The procedure requires the scene's normal map {NRM}.

    For each pixel {p=(x,y)}, let {FUN[p]} be
    {fun(NRM[0..2,p],IMG[c,p])}. The approximation {APP[p]} to be fitted
    to it is assumed to be a linear combination of the images {BAS[i][p]
    = bas[i](NRM[0..2,p])}, where {i} ranges in {0..n-1}. The system's
    unknowns {u[0..n-1]} are the coefficients of this linear
    combination. The quadratic functional, to be minimized, is the sum
    of {W[p]*(FUN[p] - APP[p])^2} where {W =
    NRM[3,p]*wht(NRM[0..2,p],IMG[c,p])}.
    
    The output matrix {A} must have space for {n^2} elements, which
    are assumed stored by rows, without gaps between the rows. The
    right-hand-side vector {b} must have space for {n} elements. */

void pst_flt_solve_lsq_system(uint32_t n, double A[], double b[], double u[], bool_t nonNegative);
  /* Solves the least-squares system {A x = b}, where {A} is an {n �
    n} matrix stored by rows, {b} is a known {n}-vector, and {u} is
    the unknown {n}-vector. If {nonNegative} is TRUE, the elements of
    {u} are constrained to be non-negative. */

void pst_flt_print_lsq_system(FILE *wr, uint32_t n, double A[], double b[]);
  /* Writes to {wr} the {n � n} matrix {A}, stored by rows, and the
    {n}-vector {b}, side by side, in human-readable format.  */

void pst_flt_apply_ambient_dimming(pst_light_t *lht, int32_t c, pst_lamp_t *src, double adim);
  /* Multiplies by {adim} the power of all lights in {lht} except {src}. */

double pst_flt_compute_ambient(pst_light_t *lht, int32_t c, pst_lamp_t *src, r3_t *nrm);
  /* Computes the apparent intensity of a pixel with surface normal
    {nrm}, when illuminated by all lamps of {lht} except {src}. */

bool_t pst_flt_pixel_is_valid(double smp, r3_t *nrm, double minNormalZ);
  /* Returns TRUE iff the pixel with image value {smp} in the relevant
    channel and surface normal {nrm} is valid. That is, {nrm} and {smp}
    must be both finite and non-zero, and the Z coordinate of {nrm} must
    be at least {minNormalZ}. */

/* ITERATIVE LIGHT LOCATION */

#define pst_fit_light_single_max_dir_step 0.2
  /* Uncertainty of initial guess about the lamp's direction. */

void pst_fit_light_single_iterative
  ( float_image_t *IMG, /* Photo of object. */ 
    float_image_t *NRM, /* Normal map of object. */
    int32_t c,          /* Channel to consider. */
    pst_light_t *lht,   /* (IN/OUT) Light model. */
    pst_lamp_t *src,    /* Lamp of {lht} to adjust. */
    double dirStep,     /* Max change allowed in {src} direction. */
    bool_t pwrAdjust,   /* TRUE optimizes the power of {src}, if possible. */
    bool_t ambAdjust,   /* TRUE optimizes the overall power of the other lamps. */
    double weightBias,  /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,  /* Ignore pixels {p} where {NRM[p].z < minNormalZ}. */
    uint32_t iterations,     /* Max iterations to use in fitting. */
    double tolerance    /* Iteration stopping criterion. */
  )
  { bool_t debug = TRUE;
    
    /* Keep {dirStep} within reasonable bounds: */
    dirStep = fmin(fmax(dirStep, 0.0), pst_fit_light_single_max_dir_step);

    /* If no parameters are adjustable, do nothing: */
    if ((dirStep <= 0.0) && (!pwrAdjust) && (!ambAdjust)) { return; }
    
    /* Get the lamp array: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps, including {src}. */
    
    /* Save initial lamp powers, since dimming is relative to them: */
    double ipwr[NS];   /* Initial lamp powers */
    for (uint32_t i = 0; i < NS; i++) { ipwr[i] = lmpv->e[i]->pwr.c[c]; }

    /* Get pointers to the adjustable parameters of {src}: */
    r3_t *sdir = &(src->dir);           /* (ADJUST) Direction of target lamp. */
    float *spwr = &(src->pwr.c[c]);   /* (ADJUST) Intensity of target lamp. */
    
    auto void print_state(uint32_t iter, double change);
      /* Prints a line with the state of the iteration: {sdir}, {spwr}, and {rdim}. */
    
    uint32_t it = 0;
    if (debug) { print_state(it, +INF); }

    /* State before iteration: */
    r3_t odir;        /* Direction computed in the previous iteration: */
    double opwr[NS];  /* Lamp powers from previous iteration. */

    double change = +INF;
    while ((it < iterations) && (change > tolerance))
      { 
        /* Save state before iteration: */
        odir = *sdir;  /* Direction computed in the previous iteration: */
        for (uint32_t i = 0; i < NS; i++) { opwr[i] = lmpv->e[i]->pwr.c[c]; }

        if (it == 0)
          { /* Make an initial heuristic adjustment of the applicable parameters: */
            bool_t dirAdjust = (dirStep > 0.0);
            pst_fit_light_single_trivial
              ( IMG, c, NRM, lht, src, dirAdjust, pwrAdjust, ambAdjust, nonNegative, minNormalZ );
          }
        else
          { /* Adjust the direction of {src}, if appropriate: */
            if (dirStep > 0.0)
              { /* Use least squares on fully lit parts. */
                pst_fit_light_single_lsq
                  ( IMG, NRM, c, lht, src, dirStep, FALSE, FALSE, weightBias, FALSE, minNormalZ );
                if (debug) { print_state(it, +INF); }
              }
            
            /* Adjust the power of {src} and the rest, if appropriate: */
            if (pwrAdjust || ambAdjust)
              { /* Reset original power of all lamps: */
                for (uint32_t i = 0; i < NS; i++) { lmpv->e[i]->pwr.c[c] = (float)ipwr[i]; }
                /* Use least squares on all valid pixels: */
                pst_fit_light_single_lsq
                  ( IMG, NRM, c, lht, src, 0.0, pwrAdjust, ambAdjust, weightBias, nonNegative, minNormalZ );
              }
          }

        it++;
        
        /* Check change in {src.dir} iteration: */
        double dirChange = acos(r3_dot(sdir, &odir)); /* Change in the direction of {dir} (radians). */

        /* Re-estimate the uncertainty in {src.dir}: */
        dirStep = fmin(2.0*dirStep, dirChange);

        /* Compute max change {pwrChange} in power levels: */
        double pwrChange_pos = 0.0; /* Total change in lamps that have increased. */
        double pwrChange_neg = 0.0; /* Total change in lamps that have decreased. */
        double maxPower = 0; /* Max power among all lamps (including {src}). */
        for (uint32_t i = 0; i < NS; i++) 
          { pst_lamp_t *alt = lmpv->e[i]; /* Lamp number {i}. */
            double apwr = alt->pwr.c[c]; /* Its power. */
            double d = apwr - opwr[i]; /* Change since prev iteration (signed). */
            if (d > 0) 
              { pwrChange_pos += fabs(d); }
            else
              { pwrChange_neg += fabs(d); }
            apwr = fabs(apwr);
            if (apwr > maxPower) { maxPower = apwr; }
          }
        /* The intensity of any pixel is a sub-convex combination of the lamp powers. */
        /* Hence the change in the pixel intensity is between the + and - tot changes. */
        double pwrChange = fmax(pwrChange_pos, pwrChange_neg);
        
        /* Compute total change, assuming that {src} has the max power seen: */
        change = pwrChange + maxPower*sin(dirChange);

        if (debug) { print_state(it, change); }
      }
    
    void print_state(uint32_t iter, double chg)
      { fprintf(stderr, "    iter %2d", iter); 
        rn_gen_print(stderr, 3, sdir->c, "%8.5f", " sdir  = ( ", " ", " )");
        fprintf(stderr, "  spwr = %8.5f", (*spwr));

        /* Compute the average power of all lamps other than {src}. */
        double apwr = 0.0;
        for (uint32_t i = 0; i < NS; i++) 
          { pst_lamp_t *alt = lmpv->e[i]; /* Lamp numbr {i}. */
            if (alt != src) { apwr += alt->pwr.c[c]; }
          }
        apwr /= NS;
        fprintf(stderr, "  avg amb pwr = %8.5f", apwr);

        if ((iter > 0) && (fabs(chg) < INF))
          { fprintf(stderr, "  change = %8.5f", chg); }
        fprintf(stderr, "\n");
      }
  }

/* LIGHT LOCATION BY SINGLE STEP LEAST SQUARES */

#define pst_fit_light_single_lsq_min_radius  M_PI/6
  /* Least squares for direction needs at least this much of the sphere. */

void pst_fit_light_single_lsq
  ( float_image_t *IMG, 
    float_image_t *NRM, 
    int32_t c,          /* Channel to consider. */
    pst_light_t *lht,   /* Light model. */
    pst_lamp_t *src,    /* Lamp of {lht} to adjust. */
    double dirStep,     /* Max change allowed in {src} direction. */
    bool_t pwrAdjust,   /* TRUE tries to optimize the power of {src}, if possible. */
    bool_t ambAdjust,   /* TRUE tries to optimize the overall power of the other lamps. */
    double weightBias,  /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ   /* Ignore image points where the normal's Z is less than this. */
  )
  { bool_t debug = TRUE;
    
    /* Get the lamp array: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps, including {src}. */
    
    /* Treat negative {dirStep} as zero: */
    dirStep = fmax(0.0, dirStep);
    
    /* Trim {dirStep} to ensure a minimum of area for reliable least squares: */
    double dirStepMax = fmax(0.0, M_PI/2 - acos(src->crad) - pst_fit_light_single_lsq_min_radius);
    dirStep = fmin(dirStep, dirStepMax);
    
    /* Turn off {ambAdjust} if there are no other lamps except {dir}: */
    if (NS == 1) { ambAdjust = FALSE; }
    
    /* If no parameters are adjustable, do nothing: */
    if ((dirStep <= 0.0) && (!pwrAdjust) && (!ambAdjust)) { return; }
    
    /* The least squares problem: */
    pst_flt_target_func_t *fun; /* The target function to be approximated. */
    uint32_t NP = 0;
    pst_flt_basis_func_t *bas[4]; /* {bas[0..NP-1]} is the function basis. */
    
    /* Get pointers to the adjustable parameters of {src}: */
    r3_t *sdir = &(src->dir);           /* (ADJUST) Direction of target lamp. */
    float *spwr = &(src->pwr.c[c]);   /* (ADJUST) Intensity of target lamp. */
    
    /* Other auxiliary data for basis functions: */
    double cfit;  /* Exclude pixels where {dot(dir,nrm)} is less than this. */

    auto double bas_nrm_x(uint32_t i, r3_t *nrm); /* Value is {nrm.x} */
    auto double bas_nrm_y(uint32_t i, r3_t *nrm); /* Value is {nrm.y} */ 
    auto double bas_nrm_z(uint32_t i, r3_t *nrm); /* Value is {nrm.z} */ 
    auto double bas_shd(uint32_t i, r3_t *nrm);   /* Value is {shade(src.dir,nrm)}. */ 
    auto double bas_amb(uint32_t i, r3_t *nrm);   /* Value is the contributin of all other lamps. */  
      /* Basis functions for least-squares fitting: */
    
    auto double fun_img(r3_t *nrm, double smp); 
      /* Target function for full model fitting: the image value {smp} itself. */
    
    auto double fun_nam(r3_t *nrm, double smp); 
      /* Target function for {src}-only fitting: the image value {smp} 
        minus effect of other lamps. */
    
    auto double wht(r3_t *nrm, double smp); 
      /* Weight function: sperical cap around current {src.dir} for power and 
        direction fitting, all valid pixels for power-only fitting. */
    
    /* Interpreting the elements of the solution vector {u}: */
    int32_t idirx, idiry, idirz;  /* {u} indices corresp. to direction coords, or -1. */
    int32_t ipwr;                 /* {u} index corresp. to {src}'s power, or -1. */
    int32_t idim;                 /* {u} index corresp. to ambient dimming factor, or -1. */
    
    /* Assemble the basis, target, and weight functions: */
    {
      /* Terms for fitting directions: */
      if (dirStep > 0.0)
        { /* It doesn't matter whether {pwrAdjust} is true or not. */
          /* If we fixed {dirStep} correctly: */
          assert(src->crad > 0.0); 
          /* Make sure that {src.dir} is a proper direction: */
          if (r3_L_inf_norm(sdir) == 0.0) { (*sdir) = (r3_t){{ 0.0, 0.0, 1.0 }}; }
          /* Add the three basis functions {nrm.x}, {nrm.y}, {nrm.z}: */
          idirx = (int32_t)NP; bas[NP] = &bas_nrm_x; NP++;
          idiry = (int32_t)NP; bas[NP] = &bas_nrm_y; NP++;
          idirz = (int32_t)NP; bas[NP] = &bas_nrm_z; NP++;
          /* Compute radius of that cap: */
          double rfit = M_PI/2 - acos(src->crad) - dirStep;
          /* If we adjusted {dirStep} properly, we should have enough room for fitting: */
          assert(rfit >= pst_fit_light_single_lsq_min_radius - 1.0e-6);
          cfit = cos(rfit);
        }
      else
        { /* No {src.dir} in result: */
          idirx = idiry = idirz = -1;
          /* Just in case: */
          cfit = -1.0; 
        }

      /* Terms for fitting the power of {src}: */
      if ((dirStep <= 0.0) && pwrAdjust)
        { /* Add the basis function {shade(nrm,dir)}: */
          ipwr = (int32_t)NP; bas[NP] = &bas_shd; NP++;
        }
      else
        { /* No {src.pwr} in result: */
          ipwr = -1;
        }

      /* Terms for fitting the ambient dimming factor: */
      if (ambAdjust)
        { /* Add the contrib of all other lamps as a basis function: */
          idim = (int32_t)NP; bas[NP] = &bas_amb; NP++;
          /* The target function is the image itself: */
          fun = &fun_img;
        }
      else
        { /* No dimming factor in result: */
          idim = -1;
          /* Target function is the image value {smp} minus the contrib of other lamps: */
          fun = &fun_nam;
        }
    } 
    
    /* The linear system and its solution: */
    double A[NP*NP];  /* Basis rigidity matrix, {A[i,j] = <bas[i],bas[j]>}. */
    double b[NP];     /* Right-hand side vector, {b[i] = <bas[i],fun>}. */
    double u[NP];     /* Solution vector, {u[i] = coefficient of {bas[i]}. */
    
    /* Build system: */
    if (debug) { fprintf(stderr, "building the least squares system...\n"); }
    pst_flt_build_lsq_system(IMG, c, NRM, NP, bas, fun, &wht, A, b);
    
    /* Solve system: */
    if (debug) { fprintf(stderr, "solving the least squares system...\n"); }
    pst_flt_solve_lsq_system(NP, A, b, u, nonNegative);

    /* Extract data from solution vector: */
    if (dirStep > 0.0)
      { /* There should be three unknowns for {src.dir*src.pwr}: */
        assert((idirx >= 0) && (idiry >= 0) && (idirz >= 0));
        /* There should be no separate unknown for {src.pwr}: */
        assert(ipwr < 0);
        /* Extract best-fit value of {src.dir*src.pwr} from solution vector: */
        r3_t udir = (r3_t){{ u[idirx], u[idiry], u[idirz] }};
        /* Split into direction and power: */
        double upwr = r3_dir(&udir, &udir);
        /* Clip to the cap of radius {dirStep} around the original dir: */
        pst_geom_clip_dir(&udir, sdir, dirStep);
        /* Set {src.dir} to the new direction: */
        (*sdir) = udir;
        if (pwrAdjust) 
          { /* Set {src}'s power to the new power: */
            (*spwr) = (float)upwr;
          }
      }
    if((dirStep <= 0.0) && pwrAdjust)
      { /* There should be a separate unknown for {src.pwr}: */
        assert(ipwr >= 0);
        /* There should no unknowns for {src.dir*src.pwr}: */
        assert((idirx < 0) && (idiry < 0) && (idirz < 0));
        /* Extract power from solution vector: */
        double upwr = u[ipwr];
        /* Set {src}'s power to the new power: */
        (*spwr) = (float)upwr;
      }
    if (ambAdjust)
      { /* There should be an unnown for the ambient dimming factor: */
        assert(idim >= 0);
        /* Get the dimming factor: */
        double udim = u[idim];
        /* Multiply it into the power of all the other lamps: */
        pst_flt_apply_ambient_dimming(lht, c, src, udim);
      }
      
    /* We are done: */
    return; 
        
    /* Basis, target, and weight functions: */

    double bas_nrm_x(uint32_t i, r3_t *nrm)
      { return nrm->c[0]; }
      
    double bas_nrm_y(uint32_t i, r3_t *nrm)
      { return nrm->c[1]; }
      
    double bas_nrm_z(uint32_t i, r3_t *nrm)
      { return nrm->c[2]; }
      
    double bas_shd(uint32_t i, r3_t *nrm)
      { return pst_lamp_geom_factor(nrm, &(src->dir), src->crad); }
      
    double bas_amb(uint32_t i, r3_t *nrm)
      { /* Shade surface with all lamps in {lht} except {src}: */
        return pst_flt_compute_ambient(lht, c, src, nrm);
      }
    
    double fun_img(r3_t *nrm, double smp)
      { return smp; }
      
    double fun_nam(r3_t *nrm, double smp)
      { return smp - pst_flt_compute_ambient(lht, c, src, nrm); }
    
    /* The weight is 1 for normal fitting, {1/(smp + weightBias)} for dark-weighted fitting. */
    
    double wht(r3_t *nrm, double smp)
      { /* Lamp direction is relevant only when direction-fitting: */
        r3_t *dirP = (dirStep == 0.0 ? NULL : &(src->dir));
        return pst_fit_light_lsq_pixel_weight(smp, nrm, minNormalZ, dirP, cfit, weightBias);
      }
  }

bool_t pst_flt_pixel_is_valid(double smp, r3_t *nrm, double minNormalZ)
  { return
     isfinite(smp) &&
     (smp != 0) &&
     isfinite(nrm->c[0]) && isfinite(nrm->c[1]) && isfinite(nrm->c[2]) && 
     (r3_L_inf_norm(nrm) != 1.0e-16) &&
     (nrm->c[2] > minNormalZ);
  }

double pst_fit_light_lsq_pixel_weight
  ( double smp,
    r3_t *nrm, 
    double minNormalZ,
    r3_t *dir,
    double minCos,
    double weightBias
  )
  { /* Weight is 1 for valid pixels within {rfit} of {src.dir}, 0 elsewhere: */
    if (! pst_flt_pixel_is_valid(smp, nrm, minNormalZ))
      { return 0.0; }
    else if ((dir != NULL) && (r3_dot(nrm, dir) < minCos))
      { return 0.0; }
    else if (fabs(weightBias) == INF)
      { return 1.0; }
    else
      { smp += weightBias;
        demand(smp > 1.0e-5, "insufficient bias for dark-weighted fitting");
        return 1.0/smp;
      }
  }

void pst_flt_apply_ambient_dimming(pst_light_t *lht, int32_t c, pst_lamp_t *src, double adim)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps, including {src}. */
    for (uint32_t i = 0; i < NS; i++) 
      { pst_lamp_t *alt = lmpv->e[i]; /* Lamp number {i}. */
        if (alt != src) 
          { float *pwrc = &(alt->pwr.c[c]);
            (*pwrc) = (float)(adim * (*pwrc)); 
          }
      }
  }

double pst_flt_compute_ambient(pst_light_t *lht, int32_t c, pst_lamp_t *src, r3_t *nrm)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps, including {src}. */
    double sum = 0.0;
    for (uint32_t i = 0; i < NS; i++) 
      { pst_lamp_t *alt = lmpv->e[i]; /* Lamp number {i}. */
        if (alt != src)
          { double gf = pst_lamp_geom_factor(nrm, &(alt->dir), alt->crad);
            sum += gf*alt->pwr.c[c];
          }
      }
    return sum;
  }

void pst_flt_build_lsq_system
  ( float_image_t *IMG, 
    int32_t c,
    float_image_t *NRM, 
    uint32_t n,
    pst_flt_basis_func_t *bas[],
    pst_flt_target_func_t *fun,
    pst_flt_weight_func_t *wht,
    double A[],
    double b[]
  )
  { bool_t debug = FALSE;
  
    int32_t NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    demand((c >= 0) && (c < NC), "invalid photo channel");
    float_image_check_size(NRM, 4, NX, NY, "bad normal map");

    /* Initialize matrix and RHS vector with zeros: */
    for (uint32_t i = 0; i < n; i++)
      { for (uint32_t j = 0; j <= n; j++) { A[i*n + j] = 0.0; }
        b[i] = 0.0;
      }
      
    double basv[n]; /* Values of basis functions at current pixel. */
    /* Scan valid pixels and add terms to the cross products: */
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { double smp = (double)float_image_get_sample(IMG, c, x, y);
            r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
            double w = pst_normal_map_get_weight(NRM, x, y)*wht(&nrm, smp);
            if (w > 0.0)
              { /* Evaluate the basis elements and target function at this pixel: */
                if (debug)
                  { fprintf(stderr, "%5d %5d smp = %7.4f", x, y, smp);
                    rn_gen_print(stderr, 3, nrm.c, "%7.4f", " nrm  = ( ", " ", " )");
                    fprintf(stderr, " bas = "); 
                  }
                for (uint32_t i = 0; i < n; i++)
                  { basv[i] = (bas[i])(i, &nrm);
                    if (debug) { fprintf(stderr, " %7.4f", basv[i]); }
                  }
                double funv = fun(&nrm, smp);
                if (debug) { fprintf(stderr, " fun = %7.4f\n", funv); }
                /* Compute the cross product terms and accumulate: */
                for (uint32_t i = 0; i < n; i++)
                  { for (uint32_t j = i; j < n; j++)
                      { /* Accumulate scalar product of basis elements {i,j}: */
                        double prod = w*basv[i]*basv[j];
                        A[i*n + j] += prod;
                        if (i != j) { A[j*n + i] += prod; }
                      }
                    /* Accumulate scalar product of basis element {i} and image: */
                    b[i] += w*funv*basv[i];
                  }
              }
          }
      }

    /* Debugging printout: */
    if (debug) { pst_flt_print_lsq_system(stderr, n, A, b); }
  }

void pst_flt_print_lsq_system(FILE *wr, uint32_t n, double A[], double b[])
  { fprintf(wr, "system =\n");
    for (uint32_t i = 0; i < n; i++)
      { fprintf(wr, "  ");
        for (uint32_t j = 0; j < n; j++) { fprintf(wr, " %11.3f", A[i*n + j]); }
        fprintf(wr, "   %11.3f", b[i]);
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }

/* TRIVIAL LIGHT MODEL FITTING */

#define pst_fit_light_trivial_min_crad  (-0.5)
  /* Trivial heuristic for direction not likely to succeed for smaller {crad}. */

void pst_fit_light_single_trivial
  ( float_image_t *IMG, 
    int32_t c,          /* Color channel to consider */
    float_image_t *NRM, 
    pst_light_t *lht,   /* Light model. */
    pst_lamp_t *src,    /* Lamp of {lht} to adjust. */
    bool_t dirAdjust,   /* TRUE estimates the direction of {src}. */
    bool_t pwrAdjust,   /* TRUE estimates the power of {src}. */
    bool_t ambAdjust,   /* TRUE estimates a dimming factor for the other lamps. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ   /* Ignore image points where the normal's Z is less than this. */
  )
  { bool_t debug = TRUE;
  
    int32_t NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    demand((c >= 0) && (c < NC), "invalid photo channel");
    float_image_check_size(NRM, 4, NX, NY, "bad normal map");
    
    /* Get the lamp array: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps, including {src}. */
    
    /* Ignore {dirAdjust} if {src} is too big: */
    if (src->crad < pst_fit_light_trivial_min_crad) { dirAdjust = FALSE; }

    /* Turn off {ambAdjust} if there are no other lamps except {src}: */
    if (NS == 1) { ambAdjust = FALSE; }
    
    /* If nothing is adjustable, we are done: */
    if ((! dirAdjust) && (! pwrAdjust) && (! ambAdjust)) { return; }

    if (ambAdjust)
      { /* Estimate the dimming factor to apply to all other lamps. */
        /* If we fixed {ambAdjust} correctly: */
        assert(NS > 1);
        /* We assume that pixels where {smp/amb} is minimum are in {src}'s shadow: */
        double adim = +INF; /* Estimated dimming factor. */
        uint32_t NV = 0; /* Count pixels used in estimate. */
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { double smp = (double)float_image_get_sample(IMG, 0, x, y);
                r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
                if (pst_flt_pixel_is_valid(smp, &nrm, minNormalZ))
                  { /* Get computed contribution from all other lamps: */
                    double amb = pst_flt_compute_ambient(lht, c, src, &nrm);
                    /* Remember smallest {smp/amb} ratio: */
                    if (amb > 0)
                      { if (smp/amb < adim) { adim = smp/amb; }
                        NV++;
                      }
                  }
              }
          }
        /* Apply the ambient dimming factor, if OK: */
        if (adim == +INF)
          { fprintf(stderr, "%s: ", __FUNCTION__);
            fprintf(stderr, "could not estimate ambient dimming\n");
          }
        else
          { if (debug)
              { fprintf(stderr, "%s: ", __FUNCTION__);
                fprintf(stderr, "amb dimming = %10.6f based on %d pixels\n", adim, NV);
              }
            if (adim != 1.0) { pst_flt_apply_ambient_dimming(lht, c, src, adim); }
          }
      }
    
    if (dirAdjust)
      { /* It doesn't matter whether {pwrAdjust} is true or not. */
        /* If we fixed {dirAdjust} correctly: */
        assert(src->crad >= pst_fit_light_trivial_min_crad + 1.0e-6);
        /* Assume that the largest differnce {smp-amb} is where {nrm=src.dir}. */
        double difMax = -INF; /* Max difference between {smp} and {amb}. */
        r3_t nrmMax = (r3_t){{ 0.0, 0.0, 1.0 }};  /* Normal where {difMax} occurred. */
        uint32_t NV = 0; /* Count pixels used in estimate. */
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { double smp = (double)float_image_get_sample(IMG, 0, x, y);
                r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
                if (pst_flt_pixel_is_valid(smp, &nrm, minNormalZ))
                  { /* Get computed contribution from all other lamps: */
                    double amb = pst_flt_compute_ambient(lht, c, src, &nrm);
                    /* The difference should be due to {src}: */
                    double dif = smp - amb;
                    /* Remember max and min difference: */
                    if ( dif > difMax) { difMax = dif; nrmMax = nrm; }
                    NV++;
                  }
              }
          }
        /* Set {src}'s direction, if OK: */
        if (difMax == -INF)
          { fprintf(stderr, "%s: ", __FUNCTION__);
            fprintf(stderr, "could not estimate direction or power\n");
          }
        else
          { if (debug)
              { fprintf(stderr, "%s: ", __FUNCTION__); 
                rn_gen_print(stderr, 3, nrmMax.c, "%8.5f", "direction = ( ", " ", " )");
                fprintf(stderr, " based on %d pixels\n", NV);
              }
            /* Set direction: */
            src->dir = nrmMax;
            if (pwrAdjust)
              { /* Set power: */
                if (debug)
                  { fprintf(stderr, "%s: ", __FUNCTION__); 
                    fprintf(stderr, "power = %10.6f based on %d pixels\n", difMax, NV);
                  }
                src->pwr.c[c] = (float)difMax;
              }
          }
      }

    if (pwrAdjust && (! dirAdjust))
      { /* Assume that the power is the average of {smp-amb} divided by the average {geo}. */
        /* where the average is weighted by the geometric shading factor {geo}. */
        double sum_dif_geo = 0.0; /* Sum of {geo*(smp-amb)}. */
        double sum_geo_geo = 0.0; /* Sum of {geo*geo}. */
        uint32_t NV = 0; /* Count pixels used in estimate. */
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { double smp = (double)float_image_get_sample(IMG, 0, x, y);
                r3_t nrm = pst_normal_map_get_vector(NRM, x, y);
                if (pst_flt_pixel_is_valid(smp, &nrm, minNormalZ))
                  { /* Get the geometric factor for {src} at this pixel: */
                    double geo = pst_lamp_geom_factor(&nrm, &(src->dir), src->crad);
                    if (geo > 0.0)
                      { /* Get computed contribution from all other lamps: */
                        double amb = pst_flt_compute_ambient(lht, c, src, &nrm);
                        /* The difference should be due to {src}: */
                        double dif = smp - amb;
                        /* accumulate sums: */
                        sum_dif_geo += dif*geo;
                        sum_geo_geo += geo*geo;
                        NV++;
                      }
                  }
              }
          }
        /* Set {src}'s power, if OK: */
        if (sum_geo_geo == 0.0)
          { fprintf(stderr, "%s: ", __FUNCTION__);
            fprintf(stderr, "could not estimate power\n");
          }
        else
          { double pwr = sum_dif_geo/sum_geo_geo;
            if (debug)
              { fprintf(stderr, "%s: ", __FUNCTION__); 
                fprintf(stderr, "power = %10.6f based on %d pixels\n", pwr, NV);
              }
            src->pwr.c[c] = (float)pwr;
          }
      }
  }
  
void pst_fit_light_multi
  ( float_image_t *IMG, 
    float_image_t *NRM,
    int32_t c,          /* Color channel to consider. */
    pst_light_t *lht,
    double weightBias,  /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ
  )
  { bool_t debug = TRUE;
  
    auto double bas_gen(uint32_t i, r3_t *nrm); 
      /* Basis functions for lamp intensity fitting (geom factor of lamp {i}). */
    
    auto double fun(r3_t *nrm, double smp); 
      /* Target function for lamp intensity fitting (the image values). */
    
    auto double wht(r3_t *nrm, double smp); 
      /* Weight for intensity fitting (1 for valid pixels, 0 for invalid). */
    
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NS = lmpv->ne; /* Number of lamps in the light model. */
    
    /* Basis functions: */
    pst_flt_basis_func_t *bas[NS];
    for (uint32_t i = 0; i < NS; i++) { bas[i] = &bas_gen; }

    /* Full linear system: */
    double A[NS*NS];  /* Basis rigidity matrix, {A[i,j] = <bas[i],bas[j]>}. */
    double b[NS];     /* Right-hand side vector, {b[i] = <bas[i],fun>}. */
    if (debug) { fprintf(stderr, "building the least squares system...\n"); }
    pst_flt_build_lsq_system(IMG, c, NRM, NS, bas, fun, &wht, A, b);
    double u[NS]; /* Solution vector. */
    
    /* Solve the system: */
    if (debug) { fprintf(stderr, "solving the least squares system...\n"); }
    pst_flt_solve_lsq_system(NS, A, b, u, nonNegative);
      
    /* Set lamp powers from {u}: */
    for (uint32_t i = 0; i < NS; i++)
      { pst_lamp_t *alt = lmpv->e[i];
        alt->pwr.c[c] = (float)u[i];
      }
    
    return;
      
    /* Basis, target, and weight functions: */
    
    double bas_gen(uint32_t i, r3_t *nrm)
      { pst_lamp_t *src = lmpv->e[i];
        return pst_lamp_geom_factor(nrm, &(src->dir), src->crad);
      }

    double fun(r3_t *nrm, double smp)
      { return smp; }

    double wht(r3_t *nrm, double smp)
      { /* Weight is 1 for valid pixels, 0 elsewhere: */
        return pst_fit_light_lsq_pixel_weight(smp, nrm, minNormalZ, NULL, 0.0, weightBias);
      }
  }

void pst_flt_solve_lsq_system(uint32_t n, double A[], double b[], double u[], bool_t nonNegative)
  {
    if (nonNegative) 
      { qmin_simplex(n, A, b, u); }
    else
      { uint32_t rank; 
        gausol_solve(n, n, A, n, b, u, TRUE,TRUE, 0.0, NULL, &rank);
        demand(rank == n, "indeterminate system");
      }
  }

/* COMMAND LINE PARSING */

void pst_fit_light_parse_weightBias(argparser_t *pp, double *weightBias)
  { if ((weightBias != NULL) && argparser_keyword_present(pp, "-weightBias"))
      { *weightBias = argparser_get_next_double(pp, 0, +INF); }
  }

void pst_fit_light_parse_iterations(argparser_t *pp, uint32_t *iterations)
  { if ((iterations != NULL) && argparser_keyword_present(pp, "-iterations"))
      { *iterations = (uint32_t)argparser_get_next_int(pp, 0, MAXINT); }
  }

void pst_fit_light_parse_tolerance(argparser_t *pp, double *tolerance)
  { if ((tolerance != NULL) && argparser_keyword_present(pp, "-tolerance"))
      { *tolerance = argparser_get_next_double(pp, 0.0, +INF); }
  }

void pst_fit_light_parse_minNormalZ(argparser_t *pp, double *minNormalZ)
  { if ((minNormalZ != NULL) && argparser_keyword_present(pp, "-minNormalZ"))
      { *minNormalZ = argparser_get_next_double(pp, -1.0, +1.0); }
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
**
**   Copyright � 2004 by the Fluminense Federal University (UFF).
**
** Created on jul/2005 by Rafael Saracchini, IC-UFF.
** Modified by Jorge Stolfi, mar/2006.
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
