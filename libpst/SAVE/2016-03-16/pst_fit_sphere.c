/* See pst_fit_sphere.h */
/* Last edited on 2013-05-24 04:29:18 by stolfilocal */ 

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
#include <affirm.h> 
#include <sve_minn.h> 

#include <pst_fit_sphere.h>
#include <pst_fit_ellipse.h>
#include <pst_camera.h>
#include <pst_geom.h>
#include <pst_sphere.h>
#include <pst_argparser.h>
#include <pst_basic.h>

/* !!! Avoid reallocations of rasterized images !!! */

/* !!! Avoid scanning the inside of the ellipse? !!! */

/* INTERNAL PROTOTYPES */

int pst_fit_sphere_num_params(double ctrAdj, double radAdj, double resAdj);
  /* Number of adjustable parameters, depending on which adjustment
    amounts are nonzero: {ctrAdj} (2 params), {radAdj} (1), {resAdj} (1). */

double pst_fit_sphere_good_spr_adjustment(r2_t *Q, double G, r2_t K, double R);
  /* Computes an increment for the spread of camera {C} that would
    cause a shift of 1 pixel in the outline {E} of a sphere. */

void pst_fit_sphere_debug(int k, char *tag, double G, r2_t K, double R, double H);
  /* Prints to stderr the index {k}, or {tag} if {k<0}, or blank if {tag} is NULL.
    Then prints the camera spread {G}, the sphere center {K}, the 
    sphere radius {R}, and the goal function value {H} (if not NAN). */

#define Pr fprintf
#define Er stderr

#define X c[0]
#define Y c[1]
#define Z c[2]
  /* Cartesian coordinates of {r2_t} or {r3_t}. */

#define HPI (0.5*M_PI)
  /* Constant {\pi/2}, used in {sigmoid} and {diomgis}. */

/* IMPLEMENTATIONS */
 
double pst_fit_sphere
  ( float_image_t *IGR, /* Gradient image of a spherical object. */  
    r2_t *Q,     /* Optical center of {IGR}. */
    double *G,    /* (IN/OUT) angular spread of camera. */
    r2_t *K,      /* (IN/OUT) Center of sphere in {IGR}. */
    double *R,    /* (IN/OUT) Radius of sphere. */
    double GAdj,  /* Maximum adjustment allowed in camera spread {G}. */
    double KAdj,  /* Maximum adjustment allowed in {K} coordinates. */
    double RAdj,  /* Maximum adjustment allowed in {R}. */
    int maxIts    /* Max iterations of the optimizer. */
  )
  { bool_t debug_pst = TRUE;
    bool_t debug_sve = FALSE;
  
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IGR, &NC, &NX, &NY);
    demand(NC == 1, "bad image depth");
    
    /* Save the initial guesses of the local parameters: */
    double GIni = (*G);
    r2_t KIni = (*K);
    double RIni = (*R);
    
    /* Paranoia: */
    demand(isfinite(Q->X), "invalid optical center X");
    demand(isfinite(Q->Y), "invalid optical center Y");
    demand(isfinite(GIni), "invalid camera spread");
    demand(isfinite(RIni) && (RIni > 0), "invalid radius");
    demand(isfinite(KAdj) && KAdj >= 0, "invalid center adjustment");
    demand(isfinite(RAdj) && RAdj >= 0, "invalid radius adjustment");
    demand(isfinite(GAdj) && GAdj >= 0, "invalid spread adjustment");
    
    /* Compute the number of parameters {NP} in optimization: */
    int NP = pst_fit_sphere_num_params(KAdj, RAdj, GAdj);
    
    auto void scatter_params(double x[], double *GT, r2_t *KT, double *RT);
      /* Maps the minimizer's argument vector {x[0..NP-1]} to the
        variable fields of {*CT,*ET}. The fixed fields of {*CT,*ET}
        are taken from from the saved values {CIni,EIni}. The decoding
        formulas are:
        
          {KT.c[j] = KMid.c[j] + x[k] * KNrm}
          
          {R = RMid * exp(x[k] * RNrm)}
          
          {G = fabs(GMid + sigmoid(x[k]) * GNrm)}
        
        The coefficients {KMid,KNrm,RMid,RNrm,GMid,GNrm}
        are such that the unpacked parameter will vary over the
        user-specified range when the corresponding {x[k]} varies over
        {[-1_+1]}. The {sigmoid} function is defined below.
        
        Note that the {fabs} in the formula for {G} does not harm
        the minimizer since the goal function is smooth and symmetric
        around zero on that parameter. */
        
    auto void gather_params(double GT, r2_t KT, double RT, double x[]);
      /* The approximate inverse of {scatter_params}. Namely, maps the
        pair {*CT,*ET} into a minimizer's argument vector
        {x[0..NP-1]}. This function is used only to encode the initial
        guess before calling the minimizer. The encoding formulas are
        
          {x[k] = (KT.c[j] - KMid.c[j])/KNrm}
          
          {x[k] = log(R/RMid)/RNrm}
          
          {x[k] = diomgis((G - GMid)/GNrm)}
        
        */
    
    auto double sigmoid(double x);
      /* A smooth sigmoid function that ranges over {[-1_+1]} as {x}
         ranges over {[-oo_+oo]}, with slope 1 at {x==0}. */
    
    auto double diomgis(double v);
      /* The inverse of {sigmoid}: ranges over {[-oo_+oo]} as {v} 
         ranges over {[-1_+1]}, with slope 1 at {v==0}. */

    /* Center coefficients for {gather_params,scatter_params}: */
    r2_t KMid = KIni;
    double KNrm = KAdj;
    
    /* Radius coefficients for {gather_params,scatter_params}: */
    double RMid = RIni;
    double RNrm = log((RMid + RAdj)/RMid);
    
    /* Spread coefficients for {gather_params,scatter_params}: */
    double GBig = 1/pst_camera_min_focal_length(Q, NX, NY);
    double GMax = fmin(GIni + GAdj, 0.8*GBig); 
    double GMin = GIni - GAdj; /* OK if negative. */
    assert(fabs(GMax) >= fabs(GMin));
    double GMid = 0.5*(GMin + GMax);
    double GNrm = 0.5*(GMin - GMax);
    if (GMax < GIni + GAdj)
      { Pr(Er, "** spread range curbed to [ %12.10f _ %12.10f ]\n", GMin, GMax); }
    if (GIni > GMid + 0.75*GNrm)
      { /* The minimizer may be too slow if started too close to {GMin,GMax}: */
        GIni = GMid + 0.75*GNrm;
        Pr(Er, "** initial spread curbed to %12.10f\n", GIni);
      }
    
    /* Count iterations for debugging and budget control: */
    int nIts = 0;
    
    auto double sve_goal(int m, double x[]);
      /* Evaluates the goal function for {sve_minn_iterate} for the 
        parameters {x[0..m-1]}.   Increments {nEvals}. 
        
        Requires {m == NP}. Calls {pst_fit_sphere_eval(IGR,ET)} with a
        geometry {ET} derived from the fixed fields of
        {GIni,KIni,RIni} with the variable parameters taken from
        {x[0..NP-1]} as per {scatter_params}. To that value, the
        procedure adds a small penalty term that favors placing the
        camera at infinity ({G == 0}). The purpose of this adjustment
        is to avoid a degenerate minimum (goal value independent of
        {G}) when the sphere is on the optical axis. */
    
    auto bool_t sve_check(int m, double x[], double Fx);
      /* To be called by the minimizer before each major optimization.
        Currently stops when the number of iterations is exceeded. */
    
    /* Gather the adjustable fields of {KIni,RIni,GIni} into the vector {x[0..NP-1]}: */
    double x[NP];  /* Packed parameters: */
    gather_params(GIni, KIni, RIni, x);
    double H = sve_goal(NP, x);
    if (debug_pst)
      { fprintf(stderr, "%s: initial goal function = %+24.16e\n", __FUNCTION__, H); }
        
    if (NP > 0)
      { /* Compute the radius of search {dMax} for {x}: */
        double dMax = sqrt(NP); /* Since each param ranges in {[-1_+1]}. */
        
        /* Call the nonlinear optimizer: */
        sve_minn_iterate
          ( /*n:*/        NP,
            /*F:*/        sve_goal,
            /*OK:*/       sve_check,
            /*x:*/        x,
            /*FxP:*/      &H,
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
    double HN = sve_goal(NP, x);
    if (debug_pst)
      { fprintf(stderr, "%s: final goal function = %+24.16e %+24.16e\n", __FUNCTION__, H, HN); }
    demand(H == HN, "inconsistent function value");

    /* Unpack the vector {x} to get the fitted parameters {ctrFin,radFin,sprFin}: */
    double GFin;
    r2_t KFin;
    double RFin;
    scatter_params(x, &GFin, &KFin, &RFin);
     
    /* Return the adjusted params and the mismatch: */
    (*G) = GFin;
    (*K) = KFin;
    (*R) = RFin;
    return H;
    
    /* IMPLEMENTATIONS OF INTERNAL PROCS */

    void scatter_params(double x[], double *GT, r2_t *KT, double *RT)
      { int k = 0;
        
        /* Get center {*KT}, from {x} or from the saved guess: */
        if (KAdj > 0) 
          { KT->c[0] = KMid.c[0] + x[k] * KNrm; k++;
            KT->c[1] = KMid.c[1] + x[k] * KNrm; k++;
          }
        else
          { (*KT) = KIni; }
          
        /* Get radius {(*RT)}, from {x} or from the saved guess: */
        if (RAdj > 0)
          { (*RT) = RMid * exp(x[k]*RNrm); k++; }
        else
          { (*RT) = RIni; }
          
        /* Get spread {(*GT)}, from {x} or from the saved guess: */
        if (GAdj > 0)
          { (*GT) = fabs(GMid + sigmoid(x[k]) * GNrm); k++; }
        else
          { (*GT) = GIni; }
        
        assert(k == NP);     
      }
        
    void gather_params(double GT, r2_t KT, double RT, double x[])
      { /* Store the variable fields of {GT,KT,RT} into {x[0..NP-1]}: */
        int k = 0;
        
        if (KAdj > 0) 
          { /* Store the center coords: */
            x[k] = (KT.c[0] - KMid.c[0])/KNrm; k++;
            x[k] = (KT.c[1] - KMid.c[1])/KNrm; k++;
          }

        if (RAdj > 0)
          { /* Store the radius: */
            x[k] = log(RT/RMid)/RNrm; k++;
          }

        if (GAdj > 0)
          { /* Store the camera spread : */
            double z = (GT - GMid)/GNrm;
            assert(fabs(z) < 0.9);
            x[k] = diomgis(z); k++;
          }

        assert(k == NP);          
      }
    
    double sigmoid(double x)
      { return atan(HPI*x)/HPI; }
    
    double diomgis(double v)
       { return tan(HPI*v)/HPI; }
    
    double sve_goal(int m, double x[])
      { /* Unpack {x} into local variables {GT,KT,RT}: */
        double GT;
        r2_t KT;
        double RT;
        scatter_params(x, &GT, &KT, &RT);
        /* Compute the sphere's projection: */
        hr3_point_t OT = pst_camera_viewpoint_from_center_spread(Q, GT);
        pst_sphere_t ST = (pst_sphere_t) { .K = KT, .R = RT };
        ellipse_crs_t ET = pst_sphere_to_ellipse(&ST, &OT);
        /* Compare with image: */
        double HT_base = pst_fit_ellipse_eval(IGR, &ET);
        /* Add a tiny penalty factor for cameras at finite distance: */
        double HT_bias = 1.0e-7 * GT*GT;
        double HT = HT_base + HT_bias;
        if (debug_sve) 
          { pst_fit_sphere_debug(-1, NULL, GT, KT, RT, HT);
            Pr(Er, "      E = "); 
            ellipse_crs_print(Er, &ET, "%7.2f");
            Pr(Er, "\n");
          }
        return HT;
      }

    bool_t sve_check(int m, double x[], double Fx)
      { nIts++;
        if (debug_sve) 
          { /* Unpack {x} into local variables {GT,KT,RT}: */
            double GT;
            r2_t KT;
            double RT;
            scatter_params(x, &GT, &KT, &RT);
            pst_fit_sphere_debug(nIts, NULL, GT, KT, RT, Fx);
          }
        return (nIts > maxIts);
      }
  }
  
int pst_fit_sphere_num_params(double KAdj, double RAdj, double GAdj)
  {
    int NP = (KAdj > 0 ? 2 : 0) + (RAdj > 0 ? 1 : 0) + (GAdj > 0 ? 1 : 0);
    return NP;
  }

double pst_fit_sphere_multiscale
  ( float_image_t *IMG, /* Photo of object. */  
    double noise,       /* Std dev of noise in {IMG}, per channel. */
    r2_t *Q,           /* Optical center of {IGR}. */
    double *G,          /* (IN/OUT) angular spread of camera. */
    r2_t *K,            /* (IN/OUT) Center of sphere in {IGR}. */
    double *R,          /* (IN/OUT) Radius of sphere. */
    double GAdj,        /* Maximum adjustment allowed in camera spread {G}. */
    double KAdj,        /* Maximum adjustment allowed in {K} coordinates. */
    double RAdj,        /* Maximum adjustment allowed in {R}. */
    double RMin,        /* Min acceptable radius for multiscale. */
    int maxIts          /* Max iterations of optimizer at initial scale. */
  )
  {
    bool_t debug = TRUE;
    
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(IMG, &NC, &NX, &NY);
    
    if (debug) { Pr(Er, "channels = %d  image size = %d × %d\n", NC, NX, NY); }
    if (debug) { pst_fit_sphere_debug(-1, "ini", *G, *K, *R, NAN); }    
    
    /* Save a local copy of the initial guesses: */
    double Gt = (*G);
    r2_t Kt = (*K);
    double Rt = (*R);
    
    /* Decide whether to recurse: */
    int NP = pst_fit_sphere_num_params(KAdj, RAdj, GAdj);
    if ((NP > 0) && (Rt + RAdj > 2*RMin))
      { /* Reduce problem to half-scale: */
        float_image_t *IMG_r = pst_fit_ellipse_image_shrink(IMG);
        int dxy = (pst_fit_ellipse_nw-1)/2;
        r2_t QCr = float_image_mscale_point_shrink(Q, dxy, dxy, pst_fit_ellipse_nw);
        double Gr;
        r2_t Kr;
        double Rr;
        pst_fit_sphere_data_shrink(Gt, Kt, Rt,  &Gr, &Kr, &Rr);
        /* Solve at a more coarse scale: */
        (void) pst_fit_sphere_multiscale
          ( IMG_r, noise, &QCr, &Gr, &Kr, &Rr,
            2*GAdj, KAdj/2, RAdj/2, 
            RMin, maxIts
          );
        /* Magnify solution to double scale, use as initial guess: */
        pst_fit_sphere_data_expand(Gr, Kr, Rr,  &Gt, &Kt, &Rt);
        /* Reset the adjustment parameters to finish off: */
        /* Assume min error of {±0.5} pixel left by recursive call, so: */
        KAdj = (KAdj == 0 ? 0 : 0.75);
        RAdj = (RAdj == 0 ? 0 : 0.75);
        double GBig = 1/pst_camera_min_focal_length(Q, NX, NY);
        double GTad = pst_fit_sphere_good_spr_adjustment(Q, Gt, Kt, Rt);
        GAdj = (GAdj == 0 ? 0 : 0.75 * fmin(GBig, GTad));
        /* We don't need many iterations to finish off: */
        maxIts = 5; 
      }

    /* Compute the relative gradient image {IGR} of original image {IMG}: */
    float_image_t *IGR = float_image_gradient_sqr_relative(IMG, noise, TRUE);

    /* Solve at present scale: */
    double H = pst_fit_sphere(IGR, Q, &Gt, &Kt, &Rt, GAdj, KAdj, RAdj, maxIts);
    if (debug) { pst_fit_sphere_debug(-1, "fin", Gt, Kt, Rt, H); }    
   
    /* Return fitted parameters to client: */
    (*G) = Gt;
    (*K) = Kt;
    (*R) = Rt;
    return H;
  }
  
void pst_fit_sphere_debug(int k, char *tag, double G, r2_t K, double R, double H)
  {
    if (k >= 0)
      { Pr(Er, "[%03d]", k); }
    else if (tag != NULL)
      { Pr(Er, "[%3s]", tag); }
    else
      { Pr(Er, " %3s ", ""); }
    Pr(Er, "  G = %14.8f", G); 
    Pr(Er, "  K = "); r2_print(Er, &K);  
    Pr(Er, "  R = %14.8f", R);
    Pr(Er, "  H = %14.8f", H);
    Pr(Er, "\n");
  }

void pst_fit_sphere_data_shrink
  ( double G,
    r2_t K,
    double R, 
    double *Gr,
    r2_t *Kr,
    double *Rr
  )
  {
    int dxy = (pst_fit_ellipse_nw-1)/2;
    (*Gr) = G / 2;
    (*Kr) = float_image_mscale_point_shrink(&K, dxy, dxy, pst_fit_ellipse_nw);
    (*Rr) = R / 2;
  }

void pst_fit_sphere_data_expand
  ( double G,
    r2_t K,
    double R, 
    double *Gx,
    r2_t *Kx,
    double *Rx
  )
  {
    (*Gx) = G * 2;
    int dxy = (pst_fit_ellipse_nw-1)/2;
    (*Kx) = float_image_mscale_point_expand(&K, dxy, dxy, pst_fit_ellipse_nw);
    (*Rx) = R * 2;
  }

double pst_fit_sphere_good_spr_adjustment(r2_t *Q, double G, r2_t K, double R)
  {
    bool_t debug = TRUE;
    /* Assemble the sphere {S} and viewpoint {O}: */
    double G0 = G;
    pst_sphere_t S = (pst_sphere_t) { .K = K, .R = R };
    hr3_point_t O = pst_camera_viewpoint_from_center_spread(Q, G0); 
    /* Compute the projection {E0} of {S}: */
    ellipse_crs_t E0 = pst_sphere_to_ellipse(&S, &O);
    /* Compute the major semiaxis {maj0} of the ellipse: */
    double maj0 = E0.rad + r2_norm(&(E0.str));
    /* Let {maj1} be the major semiaxis augmented by 1 pixel: */
    double maj1 = maj0 + 1;
    /* Get the distance {d} from the center of {S} to the optical center: */
    r2_t dsp; r2_sub(&K, Q, &dsp);
    double d = r2_norm(&dsp);
    /* Compute the camera distance {Z1} (in pixels) needed to get {maj1}: */
    double Z1 = R*(d + maj1)/sqrt(maj1*maj1 - R*R);
    /* Compute the spread corresponding to {z1}: */
    double G1 = 1/Z1;
    if (debug) { Pr(Er, "raw G0 = %12.10f  G1 = %12.10f\n", G0, G1); }
    /* assert(G1 > G0); */
    return fabs(G1 - G0);
  }
