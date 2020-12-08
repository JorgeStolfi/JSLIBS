/* See pst_gray_scale_fit.h */
/* Last edited on 2016-03-16 16:08:27 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <r2.h> 
#include <r3.h> 
#include <r3x3.h> 
#include <rn.h> 
#include <gauss_elim.h> 
#include <qmin_simplex.h> 
#include <affirm.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_gray_scale_fit.h>

/* INTERNAL DEFINITIONS */

typedef struct patch_data_t
  { double cX, cY;     /* Coordinates of center of patch (pixels). */
    double rX, rY;     /* Half-width and half-height of patch (pixels). */
    double V;          /* Mean sample value inside patch. */
    double dVdX, dVdY; /* Mean sample value gradient within patch. */
    double noise;      /* Root-mean-square residual. */
  } patch_data_t;
  /* A raw data point for photometric map fitting, extracted from
    a patch of the gray-scale or of a reference strip.
  
    The coordinates of the patch's center {(cx,cy)} are measured in
    pixels on the original chart, relative to the center of the
    first patch of the gray-scale.
    
    Field {V} is the average sample value in the patch. Fields
    {dVdX,dVdY} are the mean gradient {dV/dX,dV/dY} of {V} inside the
    patch, estimated by by least squares (linear regression). Field
    {noise} is the root-mean-square residual relative to the affine
    function with that mean and gradient. */

typedef struct diff_data_t
  { double sX, sY; /* Nominal sampling point. */
    double R_a;   /* Nominal reflectance of patch {a}. */
    double V_a;   /* Image sample value in patch {a} at {(sX,sY)}. */
    double R_b;   /* Nominal reflectance of patch {b. */
    double V_b;   /* Image sample value in patch {b} at {(sX,sY)}. */
  } diff_data_t;
  /* A differential data point for photometric map fitting.  It says that 
    two patches with reflectances {R_a} and {R_b}, under the same (unknown)
    illumination incident to point {(sX,sY)} have the `raw' sample values
    {V_a} and {V_b}, respectively, in the input image. */

vec_typedef(diff_data_vec_t,diff_data_vec,diff_data_t);
  /* A vector of differential data points. */

/* INTERNAL PROTOTYPES */

patch_data_t pgsf_get_patch_data
  ( float_image_t *img, 
    int c,
    double cX,  /* X coordinate of patch center in chart (pixels). */
    double cY,  /* Y coordinate of patch center in chart (pixels). */
    int loX,    /* Min col of patch in {img} (pixels). */
    int hiX,    /* Max col of patch in {img} (pixels). */
    int loY,    /* Min row of patch in {img} (pixels). */
    int hiY     /* Max row of patch in {img} (pixels). */
  );
  /* Extracts from channel {c} of {img} the raw data of a rectangular
    patch {{loX..hiX} × {loY..hiY}} of the chart. */
 
diff_data_t pgsf_compute_quad
  ( double R_a,
    patch_data_t *C_a,
    double R_b,
    patch_data_t *C_b
  );
  /* Computes a differential data point {ddi} from the reflectances
    {R_a,R_b} and the raw data records {C_a,C_b} of two nearby patches, by
    extrapolating their values to a point somewhere along the line
    connecting them, using the mean values and value gradients within
    the two patches as specified in {C_a,C_b}.
    
    The reflectances are merely stored into the result record.
    By convention, {-1} denotes `unknown reflectance'. */

void pgsf_get_patch_to_strip_data
  ( int c,                   /* Channel of {img} to consider. */
    int NS,                  /* Number of steps in gray scale. */
    double noise,            /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale, /* The extracted and rectified gray-scale patches. */
    double albScale[],       /* Nominal albedo of each patch. */
    double dX,               /* X displ between centers of successive patches (pixels). */
    float_image_t *imgStrip, /* Extracted image of a reference strip, or NULL. */
    double albStrip,         /* Albedo of {imgStrip0}. */
    double dY,               /* Y displ between centers of {imgScale} and {imgStrip} (pixels). */
    diff_data_vec_t *dd,     /* List of differential data points. */
    int *NG                  /* Number of differential data points. */
  );
  /* Extracts a set of differential constraints
    ({diff_data_t}s) from channel {c} of {imgScale}, by
    comparing each patch in it with the corresponding segment of the
    reference strip {imgStrip}.
    
    Assumes that patch {i} of {imgScale} has albedo {albScale[i]}, for {i} in
    {0..NS-1}. Also assumes that all patches in {imgScale} have
    the same (integer) width, and that their centers are spaced by {dX} in the
    chart. Ditto for the patches in {imgStrip}, whose albedo is assumed
    to be {albStrip}.

    The data records (one for each patch of the scale) are stored in
    the vector {dd}, starting at {dd.el[NG]}. The vector is expanded
    as needed, and {NG} is incremented by the number of data records
    collected. The procedure requires the homogeneous projective
    matrix {P} that maps chart coordianates to {img} pixel
    coordinates. */

void pgsf_get_patch_to_patch_data
  ( int c,                     /* Channel of {img} to consider. */
    int NS,                /* Number of steps in gray scale. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale,  /* The extracted and rectified gray-scale patches. */
    double albScale[],         /* Nominal albedo of each patch. */
    double dX,                 /* X displ between centers of successive patches (pixels). */
    diff_data_vec_t *dd,  /* List of differential data points. */
    int *NG                    /* Number of differential data points. */
  );
  /* Extracts a set of differential constraints ({diff_data_t}s)
    from channel {c} of {imgScale}, by comparing each patch of the gray
    scale with the next patch.
    
    Assumes that patch {i} of {imgScale} has albedo {albScale[i]},
    for {i} in {0..NS-1}. Also assumes that all patches in
    {imgScale} have the same (integer) width.

    The data records (one for each patch of the scale, except the last
    step) are stored in the vector {dd}, starting at {dd.el[NG]}. The
    vector is expanded as needed, and {NG} is incremented by the
    number of data records collected. The procedure requires the
    homogeneous projective matrix {P} that maps chart coordianates to
    {img} pixel coordinates. */

double pgsf_estimate_reflectance(diff_data_vec_t *dd, int NG0, int NG1);
  /* Estimates the reflectance {R} of a reference strip, from the
    differential data {dd[NG0..NG1-1]} collected by
    {pgsf_get_patch_to_strip_data}. Assumes that {dd[i].R_b} is the
    reflectance of patch {i} from the gray scale, {dd[i].V_b} is the
    corresponding pixel value, and {dd[i].V_a} is the pixel value on the
    reference strip adjacent to that patch. */  

void pgsf_add_lsq_term
  ( diff_data_t *ddi,
    double noise,
    pst_gray_scale_fit_basis_t B,
    int n, 
    double A[],
    double b[],
    double logVlo,          /* Log of max {V_a,V_b} in data. */
    double logVhi           /* Log of min {V_a,V_b} in data. */ 
  );
  /* Adds to the {n × n} least-squares system {A,b} a term that
    accounts for the differential data point {ddi}. */

void pgsf_solve_lsq_system(int n, double A[], double b[], bool_t nonneg, double z[]);
  /* Solves a least squares system with {n} basis elements, 
    given the {n × n} basis rigidity matrix {A} and the
    right-hand-side {n}-vector {b}, obtaining the {n}-vector {z}
    of coefficients for that basis. */

double pgsf_fudge_log(double V, double noise);
  /* Returns the natural log of {V} fudged by {noise}. 
    More precisely, computes {log(hypot(V,noise))}. */

double pgsf_nlog(double V, double noise, double logVlo, double logVhi);
  /* Returns the natural log of {V} fudged by {noise} and 
   affinely rescaled from {[logVlo, logVhi]} yo {[0_1]}. 
    More precisely, computes {log(hypot(V,noise))}. */

/* INTENSITY CORRECTION MAP MODEL
   
  An /intensity correction map/ is a function {H} that gives a
  linear-scale intensity {Y = H(V)} for each sample value {V} in some
  channel of an image {img}. The map is defined by three /rescaling
  parameteres/ {noise,logVlo,logVhi}, a /function basis/ {B}, and a
  /coefficient vector/ {z[0..B.N-1]}.
  
  Specifically, {H(V) = exp(h(nlog(V)))}; where {h} is a linear
  combination of {N} functions {bas[0..N-1]} with coefficients
  {z[0..B.N-1]}, and {nlog} is a normalized logarithm, {nlog(V) =
  (log(V)-logVlo)/(logVhi-logVlo)}. The
  specific basis {bas} and ist size {N} are determined by a
  {pst_gray_scale_fit_basis_t} record.
  
  The last basis element {bas[N-1]} is always the unit constant
  function, {bas[N-1]}; so {z[N-1]} determines an overall scaling
  factor, equal to {exp(z[N-1])}.
  
  The other elements {bas[0..N-2]} are monotonic; so if the
  coefficients {z[0..N-2]} are non-negative, the function {h} will
  be monotonic, too. For some bases, the condition is if-and-only-if,
  in which case one may want to constrain {z[0..N-2]} to be
  positive by setting {monotonic=TRUE}. */

double pst_gray_scale_fit_sigmoid_spline(double v, int p, int n);
  /* Evaluates the quadratic sigmoid spline number {p}, in a set of
    {n} such splines, at the sample value {v}. Requires {v} in {[0_1]}
    and {p} in {0..n-1}.
    
    For each {p}, the function is a C1 spline of degree 2, monotonic
    and sigmoid (a `smooth step' function). The underlying grid is a
    partition of the domain {[0_1]} into {n-1} equal intervals of
    width {s=1/(n-1)}. The function has value -1 for {v} in the
    interval {[0 _ (p-1)*s]}, value 0 at {v = p*s}, value 1 in the
    interval {[(p+1)*s _ 1]}.
    
    For {v < 0}, the function is affine (a first-degree polynomial)
    for {p=0}, and constant {-1} for other {p}. Symmetrically, for {v
    > 1}, the function is affine for {p=n-1}, and {+1} for other {p}.
    The C1 condition defines the spline from these data.
    
    The functions for {p} in {0..n-1} can reproduce exactly any linear
    or quadratic polynomial over the interval {[0_1]}, except for the
    constant term. It has the property that a linear combination is
    monotonic non-decreasing if and only if all coefficients are
    non-negative. */

double pst_gray_scale_fit_sigmoid_gaussian(double v, int p, int n);
  /* Evaluates the gaussian sigmoid spline number {p}, in a set of
    {n} such splines, at the sample value {v}. Requires {v} in {[0_1]}
    and {p} in {0..n-1}.
    
    Within the interval {[0_1]}, the function is a monotonic sigmoid
    (a `smooth step' function). It is a scaled and displaced copy of
    the error function {erf(x)}, which is the integral from 0 of a
    Gaussian distribution with zero mean and unit standard deviation.
    The mother element tends to +1 for large negative {x}, to -1 for
    large positive {x}, and has maximum slope at {x=0}, where it has
    value 0.
    
    Specifically, the function is {erf(S*((n-1)*v - p)}, where {S} is
    some fixed parameter (see implementation). The only exceptions
    are {p = 0} and {v < 0}, or {p = n-1} and {v > 1}, when the function is
    affine (a first-degree polynomial). Except fof these two cases
    (which are C1 at 0 and at 1, respectively), the function is C-oo
    at every {v} and {p}.
    
    The functions for {p} in {0..n-1} have the property that a linear
    combination {h} with non-negative coefficients is monotonic
    increasing. However, the converse is not true. Therefore,
    requiring the coefficients to be positive is a stronger constraint
    than requiring {h} to be monotonic. */

/* IMPLEMENTATIONS */

vec_typeimpl(diff_data_vec_t,diff_data_vec,diff_data_t);

void pst_gray_scale_fit_light_map
  ( int c,                     /* Channel of {img} to consider. */
    int NS,                /* Number of steps in gray scale. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale,  /* The extracted and rectified gray-scale patches. */
    double albScale[],         /* Nominal albedo of each patch. */
    double dX,                 /* X displ between centers of successive patches (pixels). */
    bool_t useSelf,            /* TRUE to use neighboring patches for lighting estimation. */
    float_image_t *imgStrip0,  /* Extracted image of a reference strip, or NULL. */
    double albStrip0,          /* Albedo of {imgStrip0}. */
    double dY0,                /* Y displ from center of {imgScale} to center of {imgStrip0} (pixels). */
    float_image_t *imgStrip1,  /* Extracted image of another reference strip, or NULL. */
    double albStrip1,          /* Albedo of {imgStrip1}. */
    double dY1,                /* Y displ from center of {imgScale} to center of {imgStrip1} (pixels). */
    pst_gray_scale_fit_basis_t B,  /* Function basis to use for fitting. */
    bool_t monotonic,     /* TRUE forces the map to be monotonic. */
    double z[],           /* (OUT) Coefficient vector. */
    double *logVlo,       /* (OUT) Low end of grid. */
    double *logVhi        /* (OUT) Hight end of grid. */ 
  )
  {
    bool_t debug = TRUE;
    
    if (debug) { fprintf(stderr, "- - - - channel %d - - - - - - - - - - - - - - - \n", c); }
    
    /* Get basis size {NB} (including unit function): */
    int NB = B.N;
    
    /* Get gray-scale dimensions: */
    int NCScale, NXScale, NYScale;
    float_image_get_size(imgScale, &NCScale, &NXScale, &NYScale);
    demand((c >= 0) && (c < NCScale), "invalid channel index");
    demand(NXScale % NS == 0, "width of {imgScale} is not a multiple of {NS}");
    
    /* Extract the differential data points from the chart: */
    int NG = 0; /* Total number of differential data points. */
    diff_data_vec_t dd = diff_data_vec_new(40);
    if (imgStrip0 != NULL)
      { pgsf_get_patch_to_strip_data
          ( c, NS, noise, imgScale, albScale, dX, imgStrip0, albStrip0, dY0, &dd, &NG );
      }
    if (imgStrip1 != NULL)
      { pgsf_get_patch_to_strip_data
          ( c, NS, noise, imgScale, albScale, dX, imgStrip1, albStrip1, dY1, &dd, &NG );
      }
    if (useSelf)
      { pgsf_get_patch_to_patch_data
          ( c, NS, noise, imgScale, albScale, dX, &dd, &NG );
      }
    diff_data_vec_trim(&dd, NG);
    
    /* Determine {logVlo,logVhi}: */
    int i;
    (*logVlo) = +INF;  (*logVhi) = -INF;
    for (i = 0; i < NG; i++)
      { diff_data_t *ddi = &(dd.e[i]);
        double logU = pgsf_fudge_log(ddi->V_a, noise);
        double logV = pgsf_fudge_log(ddi->V_b, noise);
        if (logU < (*logVlo)) { (*logVlo) = logU; }
        if (logU > (*logVhi)) { (*logVhi) = logU; }
        if (logV < (*logVlo)) { (*logVlo) = logV; }
        if (logV > (*logVhi)) { (*logVhi) = logV; }
      }

    /* The least squares system: */
    int n = NB-1;
    double A[n*n];  /* Basis rigidity matrix, {A[i,j] = <bas[i],bas[j]>}. */
    double b[n];     /* Right-hand-side vector, {b[i] = <bas[i],fun>}. */
    /* Clear {A} and {b}: */
    int p, q;
    for (p = 0; p < n; p++)
      { for (q = 0; q < n; q++) { A[p*n + q] = 0.0; }
        b[p] = 0.0; 
      }
    /* Accumulate the least squares gradient terms: */
    for (i = 0; i < NG; i++)
      { diff_data_t *ddi = &(dd.e[i]);
        pgsf_add_lsq_term(ddi, noise, B, n, A, b, *logVlo, *logVhi);
      }
      
    /* Solve the least squares system: */
    pgsf_solve_lsq_system(n, A, b, monotonic, z);
    z[NB-1] = 0.0; /* Just in case. */
    
    /* Compute the constant factor {z[NB]} so that {Vmax} maps to {Vmax}: */
    z[NB-1] = 0.0;
    double Vmax = 1.0; /* For now. */
    double val1 = pst_gray_scale_fit_eval_map(Vmax, noise, *logVlo, *logVhi, B, z);
    assert(val1 > 0.0);
    z[NB-1] = log(Vmax/val1);
  }

double pgsf_fudge_log(double V, double noise)
  { return log(hypot(V, noise)); }

double pgsf_nlog(double V, double noise, double logVlo, double logVhi)
  { double logV = pgsf_fudge_log(V, noise);
    return (logV - logVlo)/(logVhi - logVlo); 
  }

void pgsf_get_patch_to_strip_data
  ( int c,                     /* Channel of {img} to consider. */
    int NS,                /* Number of steps in gray scale. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale,  /* The extracted and rectified gray-scale patches. */
    double albScale[],         /* Nominal albedo of each patch. */
    double dX,                 /* X displacement between successive patches (pixels). */
    float_image_t *imgStrip,   /* Extracted image of a reference strip, or NULL. */
    double albStrip,           /* Albedo of {imgStrip0}. */
    double dY,                 /* Y displ between centers of {imgScale} and {imgStrip} (pixels). */
    diff_data_vec_t *dd,  /* List of differential data points. */
    int *NG                    /* Number of differential data points. */
  )
  { bool_t debug = TRUE;
  
    /* Get dimensions {DXScale,DYScale} of each patch in {imgScale}: */
    int NXScale = (int)(imgScale->sz[1]);
    int NYScale = (int)(imgScale->sz[2]);
    demand(NXScale % NS == 0, "{imgScale} width is not divisible by {NS}");
    int DXScale = NXScale/NS;
    int DYScale = NYScale;

    /* Get dimensions {DXStrip,DYStrip} of each patch in {imgStrip}: */
    int NXStrip = (int)(imgStrip->sz[1]);
    int NYStrip = (int)(imgStrip->sz[2]);
    demand(NXStrip % NS == 0, "{imgStrip} width is not divisible by {NS}");
    int DXStrip = NXStrip/NS;
    int DYStrip = NYStrip;

    if (debug) 
      { fprintf(stderr, "comparing scale patches with reference strip\n");
        fprintf(stderr, "  grayscale patch dimensions %3d × %3d\n", DXScale, DYScale);
        fprintf(stderr, "  ref strip patch dimensions %3d × %3d\n", DXStrip, DYStrip);
        fprintf(stderr, "  Y offset from grayscale to ref strip = %+6.1f pixels\n", dY);
      }

    /* Gather values and compute scale reflectances: */
    int NG0 = (*NG), NG1 = NG0; /* Data gathered here is {dd[NG0..NG1-1]}. */
    int i;
    for (i = 0; i < NS; i++)
      { /* Make sure that the next data entry exists, set {*ddi} to it: */
        diff_data_vec_expand(dd, NG1);
        diff_data_t *ddi = &(dd->e[NG1]);
        
        /* Grayscale patch domain X-range and center: */
        int loX_a = i*DXScale, hiX_a = loX_a + DXScale;
        double cX_a = i*dX, cY_a = 0.0;
        
        /* Ref strip patch domain X-range and center: */
        int loX_b = i*DXStrip, hiX_b = loX_b + DXStrip;
        double cX_b = i*dX, cY_b = dY;
        
        /* Get data from the gray-scale patch: */
        patch_data_t CScale = pgsf_get_patch_data(imgScale, c, cX_a, cY_a, loX_a, hiX_a, 0, DYScale);
        
        /* Get data from the ref strip patch: */
        patch_data_t CStrip = pgsf_get_patch_data(imgStrip, c, cX_b, cY_b, loX_b, hiX_b, 0, DYStrip);
          
        /* Combine into a differential datum: */
        (*ddi) = pgsf_compute_quad(albScale[i], &CScale, -1.0, &CStrip);
        NG1++;
      }
    (*NG) = NG1;

    /* Estimate the reflectance {R_r} of the reference strip: */
    double RStrip = pgsf_estimate_reflectance(dd, NG0, NG1);

    /* Set the reflectance  {r} of the gray background: */
    int ig;
    for (ig = NG0; ig < NG1; ig++)
      { diff_data_t *ddi = &(dd->e[ig]);
        ddi->R_a = RStrip;
        if (debug)
          { int ip = ig - NG0;
            fprintf(stderr, "  patch %02d", ip);
            fprintf(stderr, "  a: %8.5f -> %8.5f", ddi->R_a, ddi->V_a);
            fprintf(stderr, "  b: %8.5f -> %8.5f", ddi->R_b, ddi->V_b);
            fprintf(stderr, "  at ( %8.2f %8.2f )\n", ddi->sX, ddi->sY);
          }
      }
  }

void pgsf_get_patch_to_patch_data
  ( int c,                     /* Channel of {img} to consider. */
    int NS,                /* Number of steps in gray scale. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale,  /* The extracted and rectified gray-scale patches. */
    double albScale[],         /* Nominal albedo of each patch. */
    double dX,                 /* X displ between centers of successive patches (pixels). */
    diff_data_vec_t *dd,  /* List of differential data points. */
    int *NG                    /* Number of differential data points. */
  )
  { 
    bool_t debug = TRUE;
  
    /* Get dimensions {DXScale,DYScale} of each patch in {imgScale}: */
    int NXScale = (int)(imgScale->sz[1]);
    int NYScale = (int)(imgScale->sz[2]);
    demand(NXScale % NS == 0, "{imgScale} width is not divisible by {NS}");
    int DXScale = NXScale/NS;
    int DYScale = NYScale;

    /* Gather values and compute scale reflectances: */
    int NG0 = (*NG), NG1 = NG0; /* Data gathered here is {dd[NG0..NG1-1]}. */
    int i;
    for (i = 1; i < NS; i++)
      { /* Make sure that the next data entry exists, set {*ddi} to it: */
        diff_data_vec_expand(dd, NG1);
        diff_data_t *ddi = &(dd->e[NG1]);
        
        /* Get data from the previous patch: */
        double R_a = albScale[i-1]; /* Reflectance. */
        int loX_a = (i-1)*DXScale, hiX_a = loX_a + DXScale;
        double cX_a = (i-1)*dX; /* Central X. */
        patch_data_t C_a = pgsf_get_patch_data(imgScale, c, cX_a, 0.0, loX_a, hiX_a, 0, DYScale);
        /* Get data from main scale patch: */
        double R_b = albScale[i]; /* Reflectance. */
        int loX_b = (i-1)*DXScale, hiX_b = loX_b + DXScale;
        double cX_b = i*dX;   /* Central X. */
        patch_data_t C_b = pgsf_get_patch_data(imgScale, c, cX_b, 0.0, loX_b, hiX_b, 0, DYScale);
        /* Combine into a differential datum: */
        (*ddi) = pgsf_compute_quad(R_a, &C_a, R_b, &C_b);
        NG1++;
      }
    (*NG) = NG1;

    if (debug)
      { int ig;
        for (ig = NG0; ig < NG1; ig++)
          { int ip = ig - NG0 + 1;
            diff_data_t *ddi = &(dd->e[ig]);
            fprintf(stderr, "  patches %02d and %02d", ip-1, ip);
            fprintf(stderr, "  a: %8.5f -> %8.5f", ddi->R_a, ddi->V_a);
            fprintf(stderr, "  b: %8.5f -> %8.5f", ddi->R_b, ddi->V_b);
            fprintf(stderr, "  at ( %8.2f %8.2f )\n", ddi->sX, ddi->sY);
          }
      }
  }

patch_data_t pgsf_get_patch_data
  ( float_image_t *img, 
    int c,
    double cX,  /* X coordinate of patch center in chart (pixels). */
    double cY,  /* Y coordinate of patch center in chart (pixels). */
    int loX,    /* Min col of patch in {img} (pixels). */
    int hiX,    /* Max col of patch in {img} (pixels). */
    int loY,    /* Min row of patch in {img} (pixels). */
    int hiY     /* Max row of patch in {img} (pixels). */
  )
  { bool_t debug = TRUE;
    patch_data_t C;
    C.cX = cX; 
    C.cY = cY;
    C.rX = 0.5*(hiX - loX + 1); 
    C.rY = 0.5*(hiY - loY + 1); 
    
    if (debug) 
      { fprintf(stderr, "  domain [ %3d .. %3d ] × [ %3d .. %3d ]\n", loX, hiX, loY, hiY);
        fprintf(stderr, "  center ( %6.1f %6.1f ) radius ( %6.1f %6.1f )\n", C.cX, C.cY, C.rX, C.rY);
      }

    /* Sample the rectangle and compute the scalar products
        {<1|1>, <1|V>, <X|X>, <X|V>, <Y|Y>, <Y|V>, <V,V>}
      where {<f|g> = SUM{W(p)*f(p)*g(p)}} where the sum is taken over
      all the sampling grid points {p}. Assumes that {<1|X> = <1|Y> =
      <X|Y> = 0}, that is, the basis {1,X,Y} is orthogonal with
      respect to the product {<|>}.
    */
    double dot11 = 0.0, dot1V = 0.0, dotVV = 0.0;
    double dotXX = 0.0, dotXV = 0.0;
    double dotYY = 0.0, dotYV = 0.0;
    double dot1X = 0.0, dot1Y = 0.0, dotXY = 0.0; /* Jut for checking. */
    int ix, iy;
    for (iy = loY; iy <= hiY; iy++)
      { double Y = cY + iy - 0.5*(loY + hiY);
        for (ix = loX; ix <= hiX; ix++)
          { double X = cX + ix - 0.5*(loX + hiX);
            double V = float_image_get_sample(img, c, ix, iy);
            double W = 1.0; /* Weight of subsample. */
            dot11 += W;
            dot1V += W*V;
            dot1X += W*X;
            dot1Y += W*Y;
            dotXY += W*X*Y;
            dotVV += W*V*V;
            dotXX += W*X*X;
            dotXV += W*X*V;
            dotYY += W*Y*Y;
            dotYV += W*Y*V;
          }
      }
    if (debug)
      { fprintf(stderr, "  ");
        fprintf(stderr, "  <1|1> = %8.5f", dot11);
        fprintf(stderr, "  <X|X> = %8.5f", dotXX);
        fprintf(stderr, "  <Y|Y> = %8.5f", dotYY);
        fprintf(stderr, "\n");
        fprintf(stderr, "  ");
        fprintf(stderr, "  <1|X> = %8.5f", dot1X);
        fprintf(stderr, "  <1|Y> = %8.5f", dot1Y);
        fprintf(stderr, "  <X|Y> = %8.5f", dotXY);
        fprintf(stderr, "\n");
        fprintf(stderr, "  ");
        fprintf(stderr, "  <1|V> = %8.5f", dot1V);
        fprintf(stderr, "  <X|V> = %8.5f", dotXV);
        fprintf(stderr, "  <Y|V> = %8.5f", dotYV);
        fprintf(stderr, "  <V|V> = %8.5f", dotVV);
        fprintf(stderr, "\n");
      }

    /* Compute value, gradient, and variance of residual: */
    double V = dot1V/dot11;  /* Mean value. */
    double gX = dotXV/dotXX; /* Mean X derivative. */
    double gY = dotYV/dotYY; /* Mean Y derivative. */
    double var = (dotVV - (V*V*dot11 + gX*gX*dotXX + gY*gY*dotYY))/dot11;
    
    /* Save in data recrod: */
    C.V = V; C.dVdX = gX; C.dVdY = gY;
    C.noise = sqrt(fmax(0.0, var));

    if (debug)
      { fprintf(stderr, 
          "    V = %8.5f %+8.5f*dX %+8.5f*dY ± %10.7f\n", 
          C.V, C.dVdX, C.dVdY, C.noise
        );
      }

    return C;
  }

diff_data_t pgsf_compute_quad      
  ( double R_a,                  /* Reflectance of patch {a} (-1 if not known). */
    patch_data_t *C_a, /* Raw image data for patch {a}. */
    double R_b,                  /* Reflectance of patch {b} (-1 if not known). */
    patch_data_t *C_b  /* Raw image data for patch {b}. */
  )
  { bool_t debug = TRUE;
    
    diff_data_t dd;
    dd.R_a = R_a;
    dd.R_b = R_b;

    /* Find the optimum interpolation parameter {t}: */
    double dx = C_b->cX - C_a->cX;
    double dy = C_b->cY - C_a->cY;

    double wxa = dx/C_a->rX, wya = dy/C_a->rY;
    double wxb = dx/C_b->rX, wyb = dy/C_b->rY;
    
    double Ma = (wxa*wxa + wya*wya)/(C_a->rX*C_a->rY);
    double Mb = (wxb*wxb + wyb*wyb)/(C_b->rX*C_b->rY);
    
    double ta = Ma/(Ma+Mb), tb = Mb/(Ma+Mb);
    
    if (debug)
      { fprintf(stderr, "  ta = %8.5f tb = %8.5f\n",  ta, tb); }

    /* Extrapolate both patches to the point {ta*C_a->c + tb*C_b->c}. */
    dd.V_a = C_a->V + ta*(dx*C_a->dVdX + dy*C_a->dVdY);
    dd.V_b = C_b->V - tb*(dx*C_b->dVdX + dy*C_b->dVdY);

    dd.sX = ta*C_a->cX + tb*C_b->cX;
    dd.sY = ta*C_a->cY + tb*C_b->cY;

    return dd;
  }

double pgsf_estimate_reflectance(diff_data_vec_t *dd, int NG0, int NG1)
  {
    bool_t debug = TRUE;
    
    if (debug) 
      { fprintf(stderr, "  estimating the reflectance of the reference strip:\n"); }

    /* If the camera's light map was linear, the background's
      reflectance {Rr} would be {R_b[i]*V_a[i]/V_b[i]}, for any {i}. However,
      the non-linear map introduces errors in this formula, which get
      worse as the difference between {V_a[i]} and {V_b[i]} increases. On
      the other hand, if we had {V_a[i] = V_b[i]}, then {Rr} would
      obviously be {R_b[i]}. So, we find the two patches {ilo,ihi} 
      whose reflectance is closest to {Rr}, respectively from below
      and from above. We then estimate {Rr} by linear interpolating 
      between {R_b[ilo]} and {R_b[ihi]} in log scale. */
    
    assert(NG1 > NG0);
    if (NG1 - NG0 == 1) { /* The best we can do: */ return dd->e[NG0].R_b; }
    int i;
    
    /* Find indices
      {ilo} such that {V_b[ilo] <= V_a[ilo]} and {R_b[ilo]} is max,
      {ihi} such that {V_b[ihi] >= V_a[ihi]} and {R_b[ihi]} is min.
      Keep second-best in each class in {jlo,jhi}.
    */
    int ilo = -1, ihi = -1;
    int jlo = -1, jhi = -1;
    for (i = NG0; i < NG1; i++)
      { if (dd->e[i].V_b <= dd->e[i].V_a)
          { /* Main patch {i} is darker than reference strip. */
            if ((ilo < 0) || (dd->e[i].R_b > dd->e[ilo].R_b)) { jlo = ilo; ilo = i; }
            else if ((jlo < 0) || (dd->e[i].R_b > dd->e[jlo].R_b)) { jlo = i; }
          }
        if (dd->e[i].V_b >= dd->e[i].V_a)
          { /* Main patch {i} is lighter than reference strip. */
            if ((ihi < 0) || (dd->e[i].R_b < dd->e[ihi].R_b)) { jhi = ihi; ihi = i; }
            else if ((jhi < 0) || (dd->e[i].R_b < dd->e[jhi].R_b)) { jhi = i; }
          }
      }
    
    /* If background is too light or too dark, must use extrapolation: */
    if (ilo < 0) { ilo = ihi; ihi = jhi; }
    if (ihi < 0) { ihi = ilo; ilo = jlo; }
    
    /* Get the two data records: */
    diff_data_t *dlo = &(dd->e[ilo]);
    diff_data_t *dhi = &(dd->e[ihi]);
    
    if (debug) 
      { fprintf
          ( stderr, "   V_a[%d] = %8.5f   V_b[%d] = %8.5f R_b[%d] = %8.5f\n",
            ilo, dlo->V_a, ilo, dlo->V_b, ilo, dlo->R_b
          );
        fprintf
          ( stderr, "   V_a[%d] = %8.5f   V_b[%d] = %8.5f R_b[%d] = %8.5f\n",
            ihi, dhi->V_a, ihi, dhi->V_b, ihi, dhi->R_b
          );
      }

    /* Let {x[i]} be {log(V_a[i]/V_b[i])}, and {y[i]} be
      {log(R_b[i]*V_a[i]/V_b[i]) = log(R_b[i]) + x[i]}. As said above, ideally
      {y[i]} should be the constant {log(Rr)}, but in practice it is
      not, unless {x[i]} is zero. So we estimate the value that {y}
      would have if {x} was zero, by linear interpolation between
      {(x[ilo],y[ilo])} and {(x[ihi],y[ihi])}. Actually, since the `{+
      x[i]}' part interpolates to zero, we need only interpolate the
      {log(R_b[i])} */
      
    /* Interpolation abcissas: */
    double eps = 1.0e-4; /* A small perturbation. */
    double loX = log(dlo->V_a/dlo->V_b) - eps; 
    double hiX = log(dhi->V_a/dhi->V_b) + eps; 
    
    if (debug) 
      { fprintf(stderr, "    abscissas: loX = %8.5f  hiX = %8.5f\n", loX, hiX); }
    
    /* Interpolation ordinates (minus the `{+ x[i]}' part): */
    double logSlo = log(dlo->R_b);
    double logShi = log(dhi->R_b);
    
    /* Interpolate: */
    double logR = (hiX*logSlo - loX*logShi)/(hiX - loX);
    
    /* Undo the log scale: */
    double Rr = exp(logR);
      
    if (debug) 
      { fprintf(stderr, "    estimated ref strip reflectance = %8.5f\n", Rr); }
    
    return Rr;
  }

double pst_gray_scale_fit_eval_map
  ( double V, 
    double noise, 
    double logVlo, 
    double logVhi, 
    pst_gray_scale_fit_basis_t B, 
    double z[]
  )
  { 
    /* Convert {V} to normalized logscale: */
    double v = pgsf_nlog(V, noise, logVlo, logVhi);
    /* Evaluate the underlying spline: */
    double h = pst_gray_scale_fit_eval_raw_map(v, B, z);
    return exp(h);
  }

double pst_gray_scale_fit_eval_raw_map(double v, pst_gray_scale_fit_basis_t B, double z[])
  { /* Combine basis elements: */
    int p;
    double h = 0;
    for (p = 0; p < B.N; p++)
      { double baspv = pst_gray_scale_fit_eval_basis(v, p, B); 
        h += z[p]*baspv;
      }
    return h;
  }
  
void pst_gray_scale_fit_apply_map
  ( float_image_t *img,
    int c, 
    double noise, 
    double logVlo, 
    double logVhi,
    pst_gray_scale_fit_basis_t B,
    double z[]
  )
  { int NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    demand((c >= 0) && (c < NC), "invalid channel index");
    int x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { float *p = float_image_get_sample_address(img, c, x, y);
            double a = pst_gray_scale_fit_eval_map((double)(*p), noise, logVlo, logVhi, B, z);
            (*p) = (float)a;
          }
      }
  }
  
double pst_gray_scale_fit_eval_basis(double v, int p, pst_gray_scale_fit_basis_t B)
  { if (p == B.N-1)
      { /* The last element is always the unit constant function: */
        return 1.0;
      }
    else
      { /* Elements {0..B.N-2} are sigmoids: */
        switch (B.bt)
          { case pst_gray_scale_fit_btype_A: return pst_gray_scale_fit_sigmoid_spline(v, p, B.N-1);
            case pst_gray_scale_fit_btype_B: return pst_gray_scale_fit_sigmoid_gaussian(v, p, B.N-1);
            default: demand(FALSE, "bad basis type"); return 0.0;
          }
      }
  }

double pst_gray_scale_fit_sigmoid_spline(double v, int p, int n)
  { assert((p >= 0) && (p < n));
    int ni = n - 1; /* Number of intervals in {[0_1]} */
    double d0 = 2; /* Derivative of mother sigmoid at 0. */
    if (v <= 0.0)
      { /* Extrapolate with a straight line: */
        if (p == 0)
          { return d0*v*ni; }
        else
          { return -1.0; }
      }
    if (v >= 1.0)
      { /* Extrapolate with a straight line: */
        if (p == n-1)
          { return d0*(v-1)*ni; }
        else 
          { return 1.0; }
      }

    /* Now the argument {v} is in the {[0_1]} range: */
    assert((v >= 0.0) && (v <= 1.0)); 
    /* Map {v} in {[0 _ 1]} to {t} in {[0 _ ni]}: */
    double t = v*ni;
    /* Split {t} into an interval index {it} and a remainder {ft} in [0_1]: */
    int it = (int)floor(t);
    if (it >= ni) { it = ni - 1; }
    assert(it >= 0);
    assert(it < ni);
 
    double ft = t - it;
    assert(ft >= 0.0);
    assert(ft <= 1.0);
    
    /* Evaluate the mother sigmoid {sgm(it+ft - p)}: */
    if (it < p-1)
      { return -1.0; }
    else if (it+1 > p+1)
      { return +1.0; }
    else if (it < p)
      { return ft*ft - 1.0; }
    else
      { return 1.0 - (1 - ft)*(1 - ft); }
  }
  
double pst_gray_scale_fit_sigmoid_gaussian(double v, int p, int n)
  { assert((p >= 0) && (p < n));
    int ni = n - 1; /* Number of intervals in {[0_1]} */
    double S = 1.0; /* Smoothing factor */
    double d0 = S*M_2_SQRTPI;  /* Derivative of mother sigmoid at 0: S*2/sqrt(PI). */
    if (v <= 0.0)
      { /* Extrapolate with a straight line: */
        if (p == 0)
          { return d0*v*ni; }
        else
          { return -1.0; }
      }
    if (v >= 1.0)
      { /* Extrapolate with a straight line: */
        if (p == n-1)
          { return d0*(v-1)*ni; }
        else 
          { return 1.0; }
      }
    /* Now the argument {v} is in the {[0_1]} range: */
    assert((v >= 0.0) && (v <= 1.0)); 
    /* Evaluate the mother sigmoid {sgm(ni*v - p)}: */
    double t = v*ni - p;
    return erf(S*t);
  }
  
void pgsf_add_lsq_term
  ( diff_data_t *ddi,
    double noise,
    pst_gray_scale_fit_basis_t B,
    int n, 
    double A[],
    double b[],
    double logVlo,          /* Log of min {V_a,V_b} in data. */
    double logVhi           /* Log of max {V_a,V_b} in data. */ 
  )
  {
    bool_t debug = TRUE;
    
    /* Get the datum's weight (harmonic mean of reflectances): */
    double W = 2*ddi->R_a*ddi->R_b/(ddi->R_a + ddi->R_b);
    
    /* Map sample values to normalized log scale: */
    double u = pgsf_nlog(ddi->V_a, noise, logVlo, logVhi);
    double v = pgsf_nlog(ddi->V_b, noise, logVlo, logVhi);
    
    /* Compute the log {df} of the `true' reflectance ratio: */
    double df = log(ddi->R_a/ddi->R_b);

    if (debug) 
      { fprintf(stderr, "  data point: h(%+8.4f) - h(%+8.4f) = %+8.4f", u, v, df);
        fprintf(stderr, "  weight = %8.5f\n", W);
      }

    /* Add to the quadratic functional the term
    
        {Qi = W*0.5*(h(u) - h(v) - df)^2}
      
      where {W} is an importance weight for this datum.
      This means adding to each component {p} of the gradient {(A z - b)[p]}
      the term 
      
        {(SUM{q: z[q]*(bas[q](u) - bas[q](v))} - df)*(bas[p](u) - bas[p](v))}
      
      That is, we add to {A[p,q]} 
        {(bas[q](u) - bas[q](v))*(bas[p](u) - bas[p](v))}
      and add to {b[p]}
        {df*(bas[p](u) - bas[p](v))}
      both weighted by {W}. */
  
    int p, q;

    /* Evaluate {dbi[p] = bas[p](u) - bas[p](v)} for all {p}: */
    double db[n];
    for (p = 0; p < n; p++)
      { double basv = pst_gray_scale_fit_eval_basis(v, p, B);
        /* if (debug) { fprintf(stderr, "    bas[%02d](%+8.4f) = %+8.4f", p, vi, basv); } */
        double basu = pst_gray_scale_fit_eval_basis(u, p, B);
        /* if (debug) { fprintf(stderr, "    bas[%02d](%+8.4f) = %+8.4f\n", p, ui, basu); } */
        db[p] = basu - basv;
        /* if (debug) { fprintf(stderr, "    db[%02d] = %+8.4f\n", p, db[p]); } */
      }

    /* Add {GRAD(Qi)} to {A z - b}: */
    for (p = 0; p < n; p++)
      { for (q = 0; q < n; q++)
          { A[p*n + q] += W*db[p]*db[q];
            if (debug) 
              { if (fabs(db[p]*db[q]) > 0.1)
                  { /* fprintf(stderr, "    db[%02d]*db[%02d] = %+8.4f\n", p, q, db[p]*db[q]); */ }
              }
          }
        b[p] += W*df*db[p];
      }
  }

void pgsf_solve_lsq_system(int n, double A[], double b[], bool_t nonneg, double z[])
  {
    bool_t debug = TRUE;
    
    if (debug)
      { int p, q;
        fprintf(stderr, "  least squares system:\n");
        for (p = 0; p < n; p++)
          { for (q = 0; q < n; q++)
              { fprintf(stderr, " %7.3f", A[p*n + q]); }
            fprintf(stderr, " %7.3f\n", b[p]);
          }
        fprintf(stderr, "\n");
      }
    
    if (nonneg) 
      { qms_quadratic_min(n, A, b, z); }
    else
      { gsel_solve(n, n, A, 1, b, z, 0.0); }

    if (debug)
      { int p;
        fprintf(stderr, "  solution:\n");
        for (p = 0; p < n; p++)
          { fprintf(stderr, " %7.3f\n", z[p]); }
        fprintf(stderr, "\n");
      }
  }
