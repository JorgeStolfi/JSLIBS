/* See pst_fit_ellipse.h */
/* Last edited on 2016-03-16 16:09:09 by stolfilocal */ 

#define _GNU_SOURCE
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <float_image.h>
#include <float_image_gradient.h>
#include <float_image_mscale.h>
#include <argparser.h>
#include <r2.h>
#include <ellipse_crs.h> 
#include <ellipse_crs.h> 
#include <affirm.h> 
#include <sve_minn.h> 

#include <pst_fit_ellipse.h>
#include <pst_geom.h>
#include <pst_argparser.h>
#include <pst_basic.h>

/* !!! Avoid reallocations of rasterized images !!! */

/* !!! Avoid scanning the inside of the ellipse? !!! */

/* INTERNAL PROTOTYPES */

int pst_fit_ellipse_num_params(double ctrAdj, double radAdj, double strAdj);
  /* Number of adjustable parameters, depending on which adjustment ranges
    are nonzero: {ctrAdj} (2 params), {radAdj} (1 param), {strAdj}
    (2 params). */

#define Pr fprintf
#define Er stderr

/* IMPLEMENTATIONS */

double pst_fit_ellipse
  ( float_image_t *IGR, /* Gradient image of a spherical object. */  
    ellipse_crs_t *EP,  /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,      /* Maximum adjustment allowed in {EP.ctr} coordinates. */
    double radAdj,      /* Maximum adjustment allowed in {EP.rad}. */
    double strAdj,      /* Maximum adjustment allowed in {EP.str} coordinates. */
    int maxIts          /* Max iterations of the optimizer. */
  )
  { bool_t debug_pst = TRUE;
    bool_t debug_sve = FALSE;
  
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IGR, &NC, &NX, &NY);
    demand(NC == 1, "bad image depth");
    
    /* Save a local copy {E} of the initial guess: */
    ellipse_crs_t E_ini = (*EP);
    
    /* Paranoia: */
    demand(E_ini.rad > 0, "invalid radius");
    demand(ctrAdj >= 0, "invalid center adjustment");
    demand(radAdj >= 0, "invalid radius adjustment");
    demand(strAdj >= 0, "invalid stretch adjustment");
    
    /* Compute the number of parameters {NP} in optimization: */
    int NP = pst_fit_ellipse_num_params(ctrAdj, radAdj, strAdj);
    
    auto void gather_params(ellipse_crs_t *ET, double x[]);
      /* Packs the adjustable parameters of the ellipse into
        a parameter vector {x[0..NP-1]}. 
        
        The packed values are a subset of the fields of {*ET}. Each parameter is
        mapped to a number {x[k]} in such a way that its adjustment
        range maps approximately to the range {[-1_+1]}. The encodings
        are defined by these formulas:
        
          {ET.ctr.c[j] = E_ini.ctr.c[j] + x[k] * ctrAdj}
          
          {ET.rad = E_ini.rad * exp(x[k] * radNrm)}
          
          {ET.str.c[j] = E_ini.str.c[j] + x[k] * strAdj}
        
       Here {radNrm} is such that {ET.rad} is {E_ini.rad + radAdj}
       when {x[k]} is {+1}. */
    
    auto void scatter_params(double x[], ellipse_crs_t *ETP);
      /* The inverse of {gather_params}. Namely, unpacks the parameter
        vector {x[0..NP-1]} and stores its decoded elements into the
        adjustable fields of {ETP}. The fixed fields of {*ETP} are
        taken from from {E_ini}. */
        
    /* Global variables for {gather_params,scatter_params}: */
    assert(E_ini.rad != 0);
    double radNrm = log((E_ini.rad + radAdj)/E_ini.rad);
    
    /* Count iterations for debugging and budget control: */
    int nIts = 0;
    
    auto double sve_goal(int m, double x[]);
      /* Evaluates the goal function for {sve_minn_iterate} for the 
        parameters {x[0..m-1]}.   Increments {nEvals}. 
        
        Requires {m == NP}. Calls {pst_fit_ellipse_eval(IGR,ET)} with a
        geometry {ET} that derives from the fixed fields of
        {E_ini} with the variable parameters taken from
        {x[0..NP-1]} as per {scatter_params}. */
    
    auto bool_t sve_check(int m, double x[], double Fx);
      /* To be called by the minimizer before each major optimization.
        Currently stops when the number of iterations is exceeded. */
    
    /* Gather the adjustable fields of {E_ini} into the vector {x[0..NP-1]}: */
    double x[NP];  /* Packed parameters: */
    gather_params(&E_ini, x);
    double Q = sve_goal(NP, x);
    if (debug_pst)
      { fprintf(stderr, "%s: initial goal function = %+24.16e\n", __FUNCTION__, Q); }

    if (NP > 0)
      { /* Compute the radius of search {dMax} for {x}: */
        double dMax = sqrt(NP); /* Since each param ranges in {[-1_+1]}. */
        
       /* Call the nonlinear optimizer: */
        sve_minn_iterate
          ( /*n:*/        NP,
            /*F:*/        sve_goal,
            /*OK:*/       sve_check,
            /*x:*/        x,
            /*FxP:*/      &Q,
            /*dir:*/      -1,
            /*dMax:*/     dMax,
            /*rIni:*/     0.50*dMax,
            /*rMin:*/     0.05*dMax,
            /*rMax:*/     dMax,
            /*stop:*/     0.001*dMax,
            /*maxEvals:*/ maxIts,
            /*debug:*/    debug_sve
          );
      }
    if (debug_pst)
      { Pr(Er, "%s: %d iterations\n", __FUNCTION__, nIts); }
      
    /* Check the mismatch for the  final parameter vector: */
    double QN = sve_goal(NP, x);
    if (debug_pst)
      { fprintf(stderr, "%s: final goal function = %+24.16e %+24.16e\n", __FUNCTION__, Q, QN); }
    demand(Q == QN, "inconsistent function value");

    /* Unpack the vector {x} onto {E_fin}: */
    ellipse_crs_t E_fin;
    scatter_params(x, &E_fin);
      
    /* Return the adjusted params and the mismatch: */
    (*EP) = E_fin;
    return Q;
    
    /* IMPLEMENTATIONS OF INTERNAL PROCS */

    void gather_params(ellipse_crs_t *ET, double x[])
      { /* Store the variable fields of {ET} into {x[0..NP-1]}: */
        int k = 0;
        if (ctrAdj > 0) 
          { /* Store the center coords: */
            x[k] = (ET->ctr.c[0] - E_ini.ctr.c[0])/ctrAdj; k++;
            x[k] = (ET->ctr.c[1] - E_ini.ctr.c[1])/ctrAdj; k++;
          }
        if (radAdj > 0)
          { /* Store the radius: */
            x[k] = log(ET->rad/E_ini.rad)/radNrm; k++;
          }
        if (strAdj > 0)
          { /* Store the stretch vector: */
            x[k] = (ET->str.c[0] - E_ini.str.c[0])/strAdj; k++;
            x[k] = (ET->str.c[1] - E_ini.str.c[1])/strAdj; k++;
          }
        assert(k == NP);          
      }
    
    void scatter_params(double x[], ellipse_crs_t *ET)
      { /* Get the variable fields of {ET} from {x[0..NP-1]}: */
        int k = 0;
        /* Get the center coords from {x} or the saved guess: */
        if (ctrAdj > 0) 
          { ET->ctr.c[0] = E_ini.ctr.c[0] + x[k] * ctrAdj; k++;
            ET->ctr.c[1] = E_ini.ctr.c[1] + x[k] * ctrAdj; k++;
          }
        else
          {  ET->ctr = E_ini.ctr; }
          
        /* Get the radius from {x} or the saved guess: */
        if (radAdj > 0)
          { ET->rad = E_ini.rad * exp(x[k]*radNrm); k++; }
        else
          { ET->rad = E_ini.rad; }
          
        /* Get the stretch vector from {x} or the saved guess: */
        if (strAdj > 0)
          { ET->str.c[0] = E_ini.str.c[0] + x[k] * strAdj; k++;
            ET->str.c[1] = E_ini.str.c[1] + x[k] * strAdj; k++;
          }
        else
          { ET->str = E_ini.str; }
          
        assert(k == NP);     
      }
        
    double sve_goal(int m, double x[])
      { ellipse_crs_t ET;
        scatter_params(x, &ET);
        double Fx = pst_fit_ellipse_eval(IGR, &ET);
        if (debug_sve) 
          { Pr(Er, "      E = "); 
            ellipse_crs_print(stderr, &ET, "%7.2f");
            Pr(Er, " Q = %14.8f\n", Fx);
            Pr(Er, "\n");
          }
        return Fx;
      }

    bool_t sve_check(int m, double x[], double Fx)
      { 
        nIts++;
        if (debug_sve) 
          { ellipse_crs_t ET;
            scatter_params(x, &ET);
            Pr(Er, "[%03d] E = ", nIts);
            ellipse_crs_print(stderr, &ET, "%7.2f");
            Pr(Er, " Q = %14.8f\n", Fx);
            Pr(Er, "\n");
          }
        return (nIts > maxIts);
      }
  }
  
int pst_fit_ellipse_num_params(double ctrAdj, double radAdj, double strAdj)
  {
    int NP = (ctrAdj > 0 ? 2 : 0) + (radAdj > 0 ? 1 : 0) + (strAdj > 0 ? 2 : 0);
    return NP;
  }

double pst_fit_ellipse_multiscale
  ( float_image_t *IMG, /* Photo of object. */  
    double noise,       /* Std dev of noise in {IMG}, per channel. */
    ellipse_crs_t *EP,  /* (IN/OUT) Geometric parameters of sphere in {IMG}. */
    double ctrAdj,      /* Maximum adjustment allowed in {EP.ctr} coordinates. */
    double radAdj,      /* Maximum adjustment allowed in {EP.rad}. */
    double strAdj,      /* Maximum adjustment allowed in {EP.str} coordinates. */
    double minRadius,   /* Min acceptable radius for multiscale. */
    int maxIts          /* Max iterations of optimizer at initial scale. */
  )
  {
    bool_t debug = TRUE;
    
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    
    if (debug) { Pr(Er, "channels = %d  image size = %d × %d\n", NC, NX, NY); }
    if (debug) { Pr(Er, "[ini] E = "); ellipse_crs_print(Er, EP, "%7.2f"); Pr(Er, "\n"); }    
    
    /* Save a local copy {E} of the initial guess: */
    ellipse_crs_t E = (*EP);
    
    /* Decide whether to recurse: */
    int NP = pst_fit_ellipse_num_params(ctrAdj, radAdj, strAdj);
    if ((NP > 0) && (E.rad + radAdj > 2*minRadius))
      { /* Reduce problem to half-scale: */
        float_image_t *IMG_r = pst_fit_ellipse_image_shrink(IMG);
        ellipse_crs_t E_r = pst_fit_ellipse_geom_shrink(&E);
        double noise_r = noise; /* Assuming noise is fractalish. */
        /* Solve at a more coarse scale: */
        (void) pst_fit_ellipse_multiscale
          ( IMG_r, noise_r, &E_r, ctrAdj/2, radAdj/2, strAdj/2, minRadius, maxIts );
        /* Magnify solution to double scale, use as initial guess: */
        E = pst_fit_ellipse_geom_expand(&E_r);
        /* Adjust the adjustment parameters to finish off: */
        /* Assume min error of {±0.5} pixel left by recursive call, so: */
        ctrAdj = (ctrAdj == 0 ? 0 : 0.75);
        radAdj = (radAdj == 0 ? 0 : 0.75);
        strAdj = (strAdj == 0 ? 0 : 0.75);
        /* We don't need many iterations to finish off: */
        maxIts = (int)imin(maxIts, 5); 
      }

    /* Compute the relative gradient image {IGR} of original image {IMG}: */
    float_image_t *IGR = float_image_gradient_sqr_relative(IMG, noise, TRUE);

    /* Solve at present scale: */
    double Q = pst_fit_ellipse(IGR, &E, ctrAdj, radAdj, strAdj, maxIts);
    if (debug) 
      { Pr(Er, "[fin] E = \n"); 
        ellipse_crs_print(Er, &E, "%7.2f");
        Pr(Er, " Q = %+12.6f\n", Q);
        Pr(Er, "\n");
      }
      
    /* Return fitted parameters to client: */
    (*EP) = E;
    return Q;
  }
  
float_image_t *pst_fit_ellipse_image_shrink(float_image_t *IMG)
  {
    int NXR = (int)((IMG->sz[1]+1)/2);
    int NYR = (int)((IMG->sz[2]+1)/2);
    int dxy = (int)((pst_fit_ellipse_nw-1)/2);
    return float_image_mscale_shrink(IMG, NULL, NXR, NYR, dxy, dxy, pst_fit_ellipse_nw);
  }

ellipse_crs_t pst_fit_ellipse_geom_shrink(ellipse_crs_t *EP)
  {
    ellipse_crs_t E_r;
    int dxy = (pst_fit_ellipse_nw-1)/2;
    E_r.ctr = float_image_mscale_point_shrink(&(EP->ctr), dxy, dxy, pst_fit_ellipse_nw);
    E_r.rad      = 0.5 * EP->rad;
    E_r.str.c[0] = 0.5 * EP->str.c[0];
    E_r.str.c[1] = 0.5 * EP->str.c[1];
    return E_r;
  }

ellipse_crs_t pst_fit_ellipse_geom_expand(ellipse_crs_t *EP)
  {
    ellipse_crs_t E_x;
    int dxy = (pst_fit_ellipse_nw-1)/2;
    E_x.ctr = float_image_mscale_point_expand(&(EP->ctr), dxy, dxy, pst_fit_ellipse_nw);
    E_x.rad      = 2.0 * EP->rad;
    E_x.str.c[0] = 2.0 * EP->str.c[0];
    E_x.str.c[1] = 2.0 * EP->str.c[1];
    return E_x;
  }

double pst_fit_ellipse_eval(float_image_t *IGR, ellipse_crs_t *EP)
  {
    bool_t debug = FALSE;
    
    /* Get the dimensions of the gradient image: */
    demand(IGR->sz[0] == 1, "image must be monochromatic");
    int NX = (int)IGR->sz[1];
    int NY = (int)IGR->sz[2];
    
    /* Choose the half-thickness {he} of the ideal ellipse outline: */
    double he = 1.5;
    
    /* Choose half-width {hm} of the mask: */
    double hm = 2*he;
    
    /* Find a box around the ellipse, with at least {hm} pixels extra: */
    int xLo, xHi, yLo, yHi;
    ellipse_crs_int_bbox(EP, hm, &xLo, &xHi, &yLo, &yHi);
    
    /* Clip the box to the image's boundary: */
    xLo = (int)imax(0,    xLo);
    xHi = (int)imin(NX-1, xHi);
    yLo = (int)imax(0,    yLo);
    yHi = (int)imin(NY-1, yHi);
    
    if (debug) { Pr(Er, "bbox = [%d _ %d] × [%d _ %d]\n", xLo, xHi, yLo, yHi); }

    /* Convert the ellipse {EP} to a more efficient form {F}: */
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(EP, &F);
    
    /* Compute the correlation {corr} between {IGR} and the ideal gradient: */
    double sum_W_I2 = 0;
    double sum_W_E2 = 0;
    double sum_W_I_E = 0;
    int x, y;
    for (x = xLo; x < xHi; x++)
      { for (y = yLo; y < yHi; y++)
          { /* Compute the coords of the pixel center {p}: */
            double xp = x + 0.5;
            double yp = y + 0.5;
            /* Compute the coords {dp} of the pixel rel to the ellipse center: */
            r2_t dp = (r2_t){{ xp - EP->ctr.c[0], yp - EP->ctr.c[1] }};
            /* Compute the signed distance from {p} to the ellipse: */
            double drm = ellipse_ouv_border_position(&F, hm, &dp);
            if (debug) { Pr(Er, "pixel = [ %4d %4d ]  drm = %+7.4f\n", x, y, drm); }
            if (fabs(drm) < 1.0)
              { /* Compute the weight {Er} of this pixel rel to the ideal ellipse: */
                double dsm2 = (1 - drm*drm);
                double Er = dsm2*dsm2;
                
                /* Extract the value {I} from the gradient image: */
                double I = float_image_get_sample(IGR, 0, x, y);
                  
                /* Compute the value {EP} of the ideal outline at this pixel: */
                double E;
                double dre = drm*hm/he;
                if (fabs(dre) >= 1.0)
                  { E = 0.0; }
                else
                  { double dse2 = (1 - dre*dre);
                    E = dse2*dse2;
                  }
                  
                /* Accumulate the correlation sums: */
                sum_W_I2 += Er*I*I;
                sum_W_E2 += Er*E*E;
                sum_W_I_E += Er*I*E;
              }
          }
      }
    assert(sum_W_E2 > 0);
    sum_W_I2 += 1.0e-30; /* To prevent division by zero. */
    double corr = sum_W_I_E/sqrt(sum_W_E2*sum_W_I2); /* Correlation. */
    
    /* The result must be a quadratic mismatch metric, so: */
    double Q = 1 - corr;
    
    return Q;
  }
