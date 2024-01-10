/* See pst_fit_ellipse.h */
/* Last edited on 2009-02-23 11:12:56 by stolfi */ 

#define _GNU_SOURCE
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <float_image_gradient.h>
#include <argparser.h>
#include <r2.h>
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

float_image_t *pst_fit_ellipse_image_shrink(float_image_t *IMG);
  /* Produces a version {RMG} of image {IMG} reduced in size by 1/2 in
    each direction. The sample with indices {xR,yR} of {RMG}
    corresponds to the sample with indices {2*xR,2*yR} of {IMG}. Thus
    the point with coordinates {(xR,yR)} of {RMG} corresponds to the
    point with coordinates {(2*xR-0.5,2*yR-0.5) of {IMG}. */

ellipse_crs_t pst_fit_ellipse_geom_shrink(ellipse_crs_t *EP);
  /* Produces the specifications {grd} of an ellipse that results from
    the ellipse described by {EP} by the same domain reduction of
    {pst_fit_ellipse_image_shrink}. */

ellipse_crs_t pst_fit_ellipse_geom_expand(ellipse_crs_t *EP);
  /* Produces the specifications {gex} of an ellipse that results from
    the ellipse described by {EP} by the inverse of the domain reduction of
    {pst_fit_ellipse_image_shrink}. */

#define P fprintf
#define W stderr

/* IMPLEMENTATIONS */

double pst_fit_ellipse
  ( float_image_t *IMG,     /* Monochromatic photo of object. */
    ellipse_crs_t *EP,      /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,          /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,          /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,          /* Maximum ± adjustment allowed in {str} coordinates. */
    int maxIts              /* Max iterations of the optimizer. */
  )
  { bool_t debug_pst = TRUE;
    bool_t debug_sve = FALSE;
  
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    demand(NC == 1, "bad image depth");
    
    /* Compute the relative gradient image {IGR}: */
    float_image_t *IGR = float_image_gradient_sqr_relative(IMG, 0.05);
    
    /* Save a local copy {E} of the initial parameters: */
    ellipse_crs_t E = (*EP);
    
    /* Compute the number of parameters {NP} in optimization: */
    int NP = pst_fit_ellipse_num_params(ctrAdj, radAdj, strAdj);
    
    /* Global variables for {gather_params,scatter_params}: */
    double radNrm = log((E.rad + radAdj)/E.rad);
    
    /* Count iterations for debugging and budget control: */
    int nIts = 0;
    
    auto void gather_params(ellipse_crs_t *tryE, double x[]);
      /* Packs the adjustable parameters of the ellipse into a
        parameter vector {x[0..NP-1]}. The packed values are the
        adjustable fields of {*tryE}. The center and stretch are
        scaled to the range {[-1_+1]}. The radius is converted to
        relative log scale, also in {[-1_+1]}. */
    
    auto void scatter_params(double x[], ellipse_crs_t *tryE);
      /* The inverse of {gather_params}. Namely, unpacks the parameter
        vector {x[0..NP-1]} and stores its elements into the
        adjustable fields of {tryE}. The fixed fields of {*tryE} are
        taken from from {E}. */
        
    auto double sve_goal(int m, double x[]);
      /* Evaluates the goal function for {sve_minn_iterate} for the 
        parameters {x[0..m-1]}.   Increments {nEvals}. 
        
        Requires {m == NP}. Calls {pst_fit_ellipse_eval(IGR,tryE)} with a
        geometry {tryE) that combines the fixed parameters taken from
        {E} with the variable parameters taken from
        {x[0..NP-1]} as per {scatter_params}. */
    
    auto int sve_check(int m, double x[], double Fx);
      /* To be called by the minimizer before each major optimization.
        Currently stops when the number of iterations is exceeded. */
    
    if (NP > 0)
      { double x[NP];
        double dMax = sqrt(NP); /* Since each param ranges in {[-1_+1]}. */
        int maxEvals = maxIts*((NP+1)*(NP+2)/2 + 1); /* Should be enough. */

        /* Gather the adjustable params from {E} into the vector {x[0..NP-1]}: */
        gather_params(&E, x);

        /* Optimize the vector {x}: */
        sve_minn_iterate
          ( /*n:*/        NP,
            /*F:*/        sve_goal,
            /*OK:*/       sve_check,
            /*x:*/        x,
            /*rIni:*/     dMax,
            /*rMin:*/     0.05*dMax,
            /*rMax:*/     dMax,
            /*dMax:*/     dMax,  
            /*maxEvals:*/ maxEvals,
            /*debug:*/    debug_sve
          );

        /* Unpack the vector {x} onto {E}: */
        scatter_params(x, &E);
      }
      
    /* Compute the mismatch: */
    double Q = pst_fit_ellipse_eval(IGR, EP);
    if (debug_pst)
      { fprintf(stderr, "%s: %d iterations\n", __FUNCTION__, nIts); }
      
    /* Return the adjusted params and the mismatch: */
    (*EP) = E;
    return Q;
    
    /* IMPLEMENTATIONS OF INTERNAL PROCS */

    void gather_params(ellipse_crs_t *tryE, double x[])
      { /* Store the variable fields of {tryE} into {x[0..NP-1]}: */
        int k = 0;
        if (ctrAdj > 0) 
          { x[k] = tryE->ctr.c[0]/ctrAdj; k++;
            x[k] = tryE->ctr.c[1]/ctrAdj; k++;
          }
        if (radAdj > 0)
          { x[k] = log(tryE->rad/E.rad)/radNrm; k++; }
        if (strAdj > 0)
          { x[k] = tryE->str.c[0]/strAdj; k++;
            x[k] = tryE->str.c[1]/strAdj; k++;
          }
        assert(k == NP);          
      }
    
    void scatter_params(double x[], ellipse_crs_t *tryE)
      { /* Get the fixed fields of {tryE} from the initial ellipse {E}: */
        (*tryE) = E;
        /* Get the variable fields of {tryE} from {x[0..NP-1]}: */
        int k = 0;
        if (ctrAdj > 0) 
          { tryE->ctr.c[0] = x[k] * ctrAdj; k++;
            tryE->ctr.c[1] = x[k] * ctrAdj; k++;
          }
        if (radAdj > 0)
          { tryE->rad = E.rad * exp(x[k]*radNrm); k++; }
        if (strAdj > 0)
          { tryE->str.c[0] = x[k] * strAdj; k++;
            tryE->str.c[1] = x[k] * strAdj; k++;
          }
        assert(k == NP);     
      }
        
    double sve_goal(int m, double x[])
      { ellipse_crs_t tryE;
        scatter_params(x, &tryE);
        if (debug_sve) 
          { fprintf(stderr, "%s called\n", __FUNCTION__);
            pst_geom_sphere_show(stderr, &tryE, 0.0, 0.0, 0.0);
            fprintf(stderr, "\n");
          }
        return pst_fit_ellipse_eval(IGR, &tryE);
      }

    int sve_check(int m, double x[], double Fx)
      { 
        nIts++;
        if (debug_sve) 
          { ellipse_crs_t tryE;
            scatter_params(x, &tryE);
            fprintf(stderr, "-------------------------------------------------\n");
            fprintf(stderr, "%s: iteration %d\n", __FUNCTION__, nIts);
            fprintf(stderr, "current best ellipse:\n");
            pst_geom_sphere_show(stderr, &tryE, 0.0, 0.0, 0.0);
            fprintf(stderr, "\n");
            fprintf(stderr, "goal function value = %12.6f\n", Fx);
            fprintf(stderr, "-------------------------------------------------\n");
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
  ( float_image_t *IMG,  /* Monochromatic photo of object. */  
    ellipse_crs_t *EP,    /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,       /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,       /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,       /* Maximum ± adjustment allowed in {str} coordinates. */
    double minRadius,    /* Min acceptable radius for multiscale. */
    int maxIts              /* Max iterations of optimizer at initial scale. */
  )
  {
    bool_t debug = TRUE;
    
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    demand(NC == 1, "bad image depth");
    
    if (debug) { P(W, "image size = %d × %d\n", NX, NY); }
    if (debug) { P(W, "ellipse ini = "); ellipse_crs_print(W, EP); P(W, "\n"); }
    
    /* Save a local copy {E} of the initial parameters: */
    ellipse_crs_t E = (*EP);
    
    /* Decide whether to recurse: */
    int NP = pst_fit_ellipse_num_params(ctrAdj, radAdj, strAdj);
    if ((NP > 0) && (E.rad + radAdj > 2*minRadius))
      { /* Reduce problem to half-scale: */
        float_image_t *IMG_r = pst_fit_ellipse_image_shrink(IMG);
        ellipse_crs_t E_r = pst_fit_ellipse_geom_shrink(&E);
        /* Solve at a more coarse scale: */
        (void) pst_fit_ellipse_multiscale
          ( IMG_r, &E_r, ctrAdj/2, radAdj/2, strAdj/2, minRadius, maxIts );
        /* Magnify solution to double scale, use as initial guess: */
        E = pst_fit_ellipse_geom_expand(&E_r);
        /* Adjust the adjustment parameters to finish off: */
        maxIts = 5; 
        ctrAdj = fmin(1.5, ctrAdj);
        radAdj = fmin(1.5, radAdj);
        strAdj = fmin(1.5, strAdj);
      }
    
    /* Solve at present scale: */
    double Q = pst_fit_ellipse(IMG, &E, ctrAdj, radAdj, strAdj, maxIts);
    (*EP) = E;
    if (debug) { P(W, "ellipse fit = "); ellipse_crs_print(W, EP); P(W, "\n"); }
    if (debug) { P(W, "goal function = %+12.6f\n", Q); }
    return Q;
  }
  
float_image_t *pst_fit_ellipse_image_shrink(float_image_t *IMG)
  { int NC = IMG->sz[0];
    int NXI = IMG->sz[1]; int NXR = (NXI+1)/2;
    int NYI = IMG->sz[2]; int NYR = (NYI+1)/2;

    /* !!! Guessed weights --- correct values to be determined !!! */
    int hw = 1; /* Weight table has size {2*hw+1}, centered. */
    assert (2*hw+1 == 3);
    double wt[3] = { 3.0, 10.0, 3.0 }; /* Weights for resampling. */
    
    float_image_t *RMG = float_image_new(NC, NXR, NYR);
    int c;
    for (c = 0; c < NC; c++)
      { int xR, yR;
        for(yR = 0; yR < NYR; yR++)
          { for(xR = 0; xR < NXR; xR++)
              { double sum_w = 0;
                double sum_w_v = 0;
                int xD, yD;
                for(yD = -hw; yD <= +hw; yD++)
                  { int yI = 2*yR+yD;
                    double wy = ((yI < 0) || (yI >= NYI) ? 0.0 : wt[yD+hw]);
                    for(xD = -hw; xD <= +hw; xD++)
                      { int xI = 2*xR+xD;
                        double wx = ((xI < 0) || (xI >= NXI) ? 0.0 : wt[xD+hw]);
                        double w = wx*wy;
                        if (w > 0)
                          { double v = float_image_get_sample(IMG, c, xI, yI);
                            sum_w_v += w*v;
                            sum_w += w;
                          }
                      }
                  }
                double v = (sum_w == 0 ? 0.5 : sum_w_v / sum_w);
                float_image_set_sample(RMG, c, xR, yR, v);
              }
          }
      }
    return RMG;
  }

ellipse_crs_t pst_fit_ellipse_geom_shrink(ellipse_crs_t *EP)
  {
    ellipse_crs_t E_r;
    E_r.ctr.c[0] = 0.5 * (EP->ctr.c[0] + 0.5);
    E_r.ctr.c[1] = 0.5 * (EP->ctr.c[1] + 0.5);
    E_r.rad      = 0.5 * EP->rad;
    E_r.str.c[0] = 0.5 * EP->str.c[0];
    E_r.str.c[1] = 0.5 * EP->str.c[1];
    return E_r;
  }

ellipse_crs_t pst_fit_ellipse_geom_expand(ellipse_crs_t *EP)
  {
    ellipse_crs_t E_x;
    E_x.ctr.c[0] = 2.0 * EP->ctr.c[0] - 0.5;
    E_x.ctr.c[1] = 2.0 * EP->ctr.c[1] - 0.5;
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
    int NX = IGR->sz[1];
    int NY = IGR->sz[2];
    
    /* Choose the half-thickness {he} of the ideal ellipse outline: */
    double he = 1.5;
    
    /* Choose half-width {hm} of the mask: */
    double hm = 2*he;
    
    /* Find a box around the ellipse, with at least {hm} pixels extra: */
    int xLo, xHi, yLo, yHi;
    ellipse_crs_int_bbox(EP, hm, &xLo, &xHi, &yLo, &yHi);
    
    /* Clip the box to the image's boundary: */
    xLo = imax(0,    xLo);
    xHi = imin(NX-1, xHi);
    yLo = imax(0,    yLo);
    yHi = imin(NY-1, yHi);
    
    if (debug) { P(W, "bbox = [%d _ %d] × [%d _ %d]\n", xLo, xHi, yLo, yHi); }

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
            if (debug) { P(W, "pixel = [ %4d %4d ]  drm = %+7.4f\n", x, y, drm); }
            if (fabs(drm) < 1.0)
              { /* Compute the weight {W} of this pixel rel to the ideal ellipse: */
                double dsm2 = (1 - drm*drm);
                double W = dsm2*dsm2;
                
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
                sum_W_I2 += W*I*I;
                sum_W_E2 += W*E*E;
                sum_W_I_E += W*I*E;
              }
          }
      }
    assert(sum_W_E2 > 0);
    sum_W_I2 += 1.0e-30; /* To prevent division by zero. */
    double corr = sum_W_I_E/sqrt(sum_W_E2*sum_W_I2); /* Correlation. */
    
    /* The result must be a quadratic mismatch metric, so: */
    return 1 - corr;
  }
