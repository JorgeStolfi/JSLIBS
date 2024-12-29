/* Last edited on 2024-12-21 14:05:44 by stolfi */
#define PROG_NAME "virtual_gauge"
#define PROG_DESC "fits a simple lighting model to a spherical light gauge image"

#define virtual_gauge_C_COPYRIGHT \
  "Copyright © 2004 by the Federal Fluminense University (UFF) and State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <jsfile.h>
#include <math.h>
#include <normais.h>
#include <r3.h>
#include <rn.h>
#include <rmxn.h>
#include <gausol_solve.h>
#include <assert.h>
#include <affirm.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <string.h>
#include <tabela.h>
#include <argparser.h>
#include <lighting_models.h>
#include <lighting_compact.h>
#include <lighting_radial.h>
#include <lighting_harmonic.h>
#include <lighting_stereopoly.h>
#include <lighting_nearest.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "      -gauge \\\n" \
  "          image {IMAGE_FNAME} \\\n" \
  "          center {IX} {IY} \\\n" \
  "          radius {IRADIUS} \\\n" \
  "          [ albedo {IDR} [ {IDG} {IDB} ] \\\n" \
  "          [ stretch {ISX} {ISY} ] \\\n" \
  "          [ view {IVX} {IVY} {IVZ} ] \\\n" \
  "          [ trim {TRIMVALUE} ] \\\n" \
  "          [ gamma {GAMMA} ] \\\n" \
  "      -model \\\n" \
  "      [harmonic degree {DEGREE} ] \\\n" \
  "      [stereopoly degree {DEGREE} ] \\\n" \
  "      [radial resolution {RESOLUTION} span {SPAN} \\\n" \
  "      [compact \\\n" \
  "           dirlight {Ux Uy Uz} [error {LIGHTERROR}] \\\n" \
  "           radlight {RADLIGHT} [error {RADERROR}] \\\n" \
  "           shine {SHINE} [error {SHINERROR}] \\\n" \
  "           [isotropic] \\\n" \
  "           [backlight {Xix Xiy Xiz} error {BACKLIGHTERROR}] \\\n" \
  "      ]\\\n" \
  "      [ -maskOut {MASK_FNAME} ] \\\n" \
  "      [ -paramsOut {PARAMS_FNAME} ] \\\n" \
  "      [ -gammaOut {GAMMAOUT} ] \\\n" \
  "      [ -ignorePosWeight ] \\\n" \
  "      [ -useCosPosWeight ] \\\n" \
  "      [ -virtual \\\n" \
  "          image {OIMAGE_FNAME} \\\n" \
  "          [ magnify {OMAG} ] \\\n" \
  "          [ albedo {ODR} [ {ODG} {ODB} ] \\\n" \
  "          [ view {OVX} {OVY} {OVZ} ] \\\n" \
  "      ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program fits a simple light flow model to an image" \
  " showing a spherical Lambertian light gauge.  Optionally, it" \
  " creates a synthetic image of a \"virtual\" spherical" \
  " gauge from that lighting model.\n" \
  "\n" \
  "  The light flow model is fitted independently for each" \
  " color channel.  It consists of a `direct' term" \
  " due to a distant point source and an `indirect' term" \
  " due to stray and scene-scattered light.  By default, the" \
  " indirect term is an isotropic ambient light field.  If" \
  " the \"-backPlane\" option is given, the indirect field" \
  " is that of an uniformly illuminated Lambertian plane" \
  " perpendicular to the Z axis and located behind the gauge.\n" \
  "\n" \
  "  The \"-gauge\" option specifies the name {IGAUGE_FNAME} of the" \
  " image file containing the gauge's photo, the geometry of the" \
  " gauge's projection in it, the gauge's albedo, and the viewing" \
  " direction.  The geometry is described by the" \
  " sub-parameters \"radius\", \"center\", and \"stretch\", as" \
  " in {ellipse_crs.h}.  The viewing direction is a vector that" \
  " points from the gauge towards the camera; only its" \
  " direction matters, and is interpreted as the surface" \
  " normal at the point of the gauge that is visible at" \
  " the exact center of the projection.   The gauge's albedo" \
  " {IDR} {IDG} {IDB} is not critical; it only affects the" \
  " estimated light intensities by a constant factor.\n" \
  "\n" \
  "  If \"-maskOut\" is specified, the" \
  " program will write to file {MASK_FNAME} a grayscale mask" \
  " showing which pixels of the input file were considered to" \
  " be inside the gauge's projection (value 1), outside it" \
  " (value 0), or straddling the projection's boundary" \
  " (fractions between 0 and 1). \n" \
  "\n" \
  "  If \"-paramsOut\" is specified, the" \
  " program will write to file {PARAMS_FNAME} containing one line" \
  " for each color channel.  Each line" \
  " gives the channel number {c} and the corresponding fitted" \
  " light flow parameters: the three components" \
  " {FX} {FY} {FZ} of the point source's direction, its" \
  " intensity {ED}, the intensity {EB} of the light" \
  " scattered by the backplane, and the intensity {EI} of" \
  " the isotropic `ambient' light. \n" \
  "\n" \
  "  If \"-virtual\" is specified, the program will write to the" \
  " file {OGAUGE_FNAME} a synthetic image of a gauge illuminated" \
  " with the fitted light model, against a background with" \
  " radiance {EB}.  The output image will have the same size" \
  " as the input, and the gauge's projection in it will have" \
  " the same geometry (radius, center and stretch) as in the" \
  " input image; except that the whole image will be magnified by" \
  " the integer factor {OMAG} (default 1) with an extra margin {OMARG} pixels" \
  " all around (default 0).   The virtual gauge will be" \
  " assumed to have albedo {ODR,ODG,ODB}.\n" \
  "\n" \
  "  The input and output images may have extension \".pgm\", \".ppm\", or" \
  " \".fni\".  The parameter {GAMMA} applies to \".pgm\" or \".ppm\" images" \
  " only.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  The Taj Mahal and the Great Pyramid of Giza.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in 2004 (?) by Rafael F. V Saracchini (IC-UFF).\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2004-2009  Extensive changes by J. Stolfi and R. Saracchini at IC-UNICAMP.\n" \
  "  2010       Ambient term added by R. Saracchini at UWE-Bristol.\n" \
  "  2010-03-17 Thorough rewrite by J. Stolfi at IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " virtual_gauge_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  




typedef struct light_field_t
  { r3_t f;         /* Direction of point source. */
    double Ed;      /* Intensity of point source term. */
    double Eb;      /* Intensity of backplane term. */
    double Ei;      /* Intensity of isotropic term. */
  } light_field_t;
  /* Describes a monocheomatic lighting field that consists of direct
    light from a point source in direction {f} with intensity {Ed},
    light scattered from a Z-perpendicular backplane with radiance
    {Eb}, and isotropic ambient light with intensity {Ei}. */

typedef struct options_t
  { gauge_data_t *gin;     /* Data about input gauge image. */
    double gammaOut;          /* Gamma for  output images. */
    
    char *maskOut;         /* Name of output mask image file, or {NULL}. */
    char *maskIn;         /* Name of input mask image file, or {NULL}. */
    char *paramsOut;       /* Name of output parameter file, or {NULL}. */
    bool_t virtual;        /* TRUE if synthetic image is required. */
    bool_t virtualUserDefined;        /* TRUE if user wants define the virtual output image params */
    gauge_data_t *got;     /* Data about the output gauge image. */
    bool_t estimateU;
    int nStepsU;
    bool_t useNonLinear;
    int nStepsNL;
    
    
    model_options_t lmodel;
    
    
    char* errorOut;
    int nSteps;
    bool_t ignorePosWeight;
    bool_t useCosPosWeight;
  } options_t;
  /* Options from the command line. */

void ensure_non_negative_light(double *Evalue, char *Ename);
  /* Makes sure that the light intensity {*Evalue} is nonnegative.
    Bombs out if too negative. */

double compute_shade(double D, r3_t *n, light_field_t *L);
  /* Shading function for a single channel of the lighting model.
    
    Computes the apparent radiance for a given Lambertian surface with
    albedo {D} and normal {n} under a monochromatic light flow {L}.
    Assumes that the backplane is perpendicular to the Z axis. */

double compute_shade(double D, r3_t *n, light_field_t *L)
  {
    bool_t debug = FALSE;
    /* Light flow terms: */
    double val_d = L->Ed*r3_dot(&(L->f),n);
    if(val_d < 0 ) val_d = 0;
    double val_b = L->Eb*(0.5*(1 - n->c[2]));
    double val_i = L->Ei;
    /* Surface radiance: */
    double rad = D*(val_d + val_b + val_i);
    if (debug)
      { r3_gen_print(stderr, n, "%+9.6f", "  n = (", " ", ")");
        fprintf(stderr, "  direct = %8.6f  backplane = %8.6f  isotropic = %8.6f", val_d, val_b, val_i);
        fprintf(stderr, "  radiance = %8.6f\n", rad);
      }
    return rad;
  }

typedef double basis_function(int i, r3_t *X) ;

void build_least_squares_system
  ( double A[],
    double B[],
    basis_function *phi,
    int NB,
    r3_t X[], 
    double F[],
    int NP,
    bool_t sel[]
  );
  /* Given a spherical function basis {phi(i,x} with {NB} elements, a
    list {X[0..NP-1]} of {NP} points of the sphere, and a list
    {F[0..NP-1]} of data values at those points, stores into {A} and
    {B} the matrix and the RHS vector, respectivey, of the least
    squares system for approximating {F} with the basis {phi}.
    
    Uses the function dot product {<F,G>} defined as the sum of
    {F(X[i])*G(X[i])} over every {i} such that {sel[i]} is true.
*/

void build_least_squares_system
  ( double A[],
    double B[],
    basis_function *phi,
    int NB,
    r3_t X[], 
    double F[],
    int NP,
    bool_t sel[]
  )
  { 
    int r, s, i;
    bool_t debug = FALSE;

    if (debug) { fprintf(stderr,"\n\n"); }
    for ( r = 0; r < NB; r++)
      { for ( s = 0; s < NB; s++)
          { double valueA = 0;
            for (i = 0; i < NP; i++)
              { if (sel[i])
                  { r3_t *Xi = &(X[i]);
                    double phiRi = phi(r,Xi);
                    double phiSi = phi(s,Xi);
                    valueA+= phiRi*phiSi;
                  }
              }
            A[(r*NB) + s] = valueA;
          }

        double valueB = 0;
        for (i = 0; i < NP; i++)
          { if (sel[i])
              { r3_t *Xi = &(X[i]);
                double phiRi = phi(r,Xi);
                double Fi = F[i];
                valueB += phiRi*Fi;
              }
          }
        B[r] = valueB;

        if (debug)
          { int done = r + 1;
            fprintf(stderr,"\033[1A");
            fprintf(stderr,"Processed [%04d of %04d] - %4.3f%%\n",done,NB,done*100.0/((double)NB));
          }
    }
        
  }

void compute_lighting_parameters
  ( double D,          /* Assumed albedo of gauge. */
    r3_t X[],          /* Normals at gauge sampling points. */
    double F[],        /* Observed radiance at those points. */
    int NP,            /* Number of sampling points. */
    bool_t BPL,        /* Indirect light term: TRUE for backplane, FALSE for isotropic. */
    light_field_t *LP  /* (OUT) Lighting parameters. */
  );
  /* Given a set of surface normals {X[0..NP-1]}, the corresponding
    monochromatic radiances {F[0..NP-1]}, and the surface albedo {D},
    computes the monochromatic light flow parameters for {compute_shade} 
    and stores them into {*LP}.  
    
    The {BPL} flag specifies the form of the indirect lighting terms,
    due to stray "ambient" light. If {BPL} is true, the procedure
    assumes that the indirect lighting comes entirely from a
    uniformly-colored plane perpendicular to the Z axis and located
    behind the gauges. In that case {LP->Ei} will be set to zero. If
    {BPL} is false, the procedure assumes that the indirect lighting
    is isotropic, and will set {LP->Eb} to zero. */

void compute_lighting_parameters
  ( double D,     /* Assumed albedo of gauge. */
    r3_t X[],     /* Normals at gauge sampling points. */
    double F[],   /* Observed radiance at those points. */
    int NP,       /* Number of sampling points. */
    bool_t BPL,   /* Indirect light term: TRUE for backplane, FALSE for isotropic. */
    light_field_t *LP  /* (OUT) Lighting parameters. */
  )
  {
    int NB = 4;  /* Number of elements in function basis. */
    
    r3_t f = (r3_t){{ 0,0,1 }}; /* Current guess for light direction.*/

    auto double phi(int r, r3_t *X);

    double phi(int r, r3_t *X)
      { if(r < 3)
          { /* Elements for direct lighting by point source: */
            double s = r3_dot(&f, X);
            return (s <= 0 ? 0 : X->c[r]);
          }
        else if (r == 3)
          { /* Element for indirect lighing: */
            return (BPL ? 0.5*(1 - X->c[2]) : 1);
          }
        else
          { /* No such element: */
            assert(FALSE);
          }
      }
    
    /* The least-squares system: */
    double *A = rmxn_alloc(NB,NB);  /* System's matrix. */
    double *B = rmxn_alloc(NB,1);   /* System's RHS vector. */
    double *C = rmxn_alloc(NB,1);   /* System's solution (basis coeffs). */

    /* Boolean vector that selects the points for the scalar product: */
    bool_t *sel = notnull(malloc(NP*sizeof(bool_t)), "no mem");

    /* Stop when the computed source direction changes by less than this amount: */
    double epsilon = 10e-6;

    /* Main loop for adjustment of {f}: */
    int MAX_ITERATIONS = 2000;
    int iter = 0;
    double Ed, Eb, Ei;  /* Computed intensities of lighting field components. */
    double diff;        /* Change in source direction in last iteration. */
    do
      { 
        fprintf(stderr, "iteration %d", iter); 
        r3_gen_print(stderr, &f, "%+9.6f", "  f = (", " ", ")");
        
        /* Mark in {sel[0..NP]} the gauge points illuminated by the point source: */
        int i;
        int nok = 0;
        for (i = 0; i < NP; i++) 
          { r3_t *Xi = &(X[i]);
            sel[i] = (r3_dot(&f, Xi) > 0);
            if (sel[i]) { nok++; }
          }
        fprintf(stderr, " %d out of %d illuminated points\n", nok, NP);
        
        /* Make sure that we have enough data: */
        if (nok < 3*NB)
          { fprintf(stderr, "** not enough illuminated points - aborted\n");
            exit(1);
          }

        /* Build and solve the least squares system: */
        build_least_squares_system(A, B, &phi, NB, X, F, NP, sel);
        gausol_solve(NB, NB, A, 1, B, C, TRUE,TRUE, 1.0e-14, NULL,NULL);
        
        /* Check for invalid values: */
        int r;
        for (r = 0; r < NB; r++) { assert(isfinite(C[r])); }
        
        /* Extract the lighting parameters from the basis coefficients: */
        r3_t f_new = (r3_t){{ C[0], C[1], C[2] }};
        double M = r3_dir(&f_new, &f_new);
        Ed = M/D;
        if (BPL)
          { Eb = C[3]/D; Ei = 0.0; }
        else
          { Eb = 0.0; Ei = C[3]/D; }
        r3_gen_print(stderr, &f_new, "%+9.6f", "  f(new) = (", " ", ")");
        fprintf(stderr,"  Ed = %+-14.7e  Eb = %+-14.7e  Ei = %+-14.7e", Ed, Eb, Ei);

        /* Compare with previous light direction: */
        diff = r3_dist(&f, &f_new);
        fprintf(stderr,"  diff = %10.7f\n", diff);

        /* Update {f} and proceed: */
        f = f_new;
        iter++;
        fprintf(stderr, "\n");
      }
    while((diff > epsilon) && (iter < MAX_ITERATIONS));
    fprintf(stderr, "did %d iterations\n", iter);
    
    if (diff <= epsilon)
      { fprintf(stderr, "converged!\n"); }
    else
      { fprintf(stderr, "** warning: failed to converge\n"); } 
      
    fprintf(stderr, "light field:"); 
    fprintf(stderr, "  direct = %8.6f", Ed);
    r3_gen_print(stderr, &f, "%+9.6f", "  from (", " ", ")");
    fprintf(stderr, "  backplane = %8.6f  isotropic = %8.6f\n", Eb, Ei);
      
    /* Basic sanitation and sanity check: */
    ensure_non_negative_light(&Ed, "Ed");
    ensure_non_negative_light(&Eb, "Eb");
    ensure_non_negative_light(&Ei, "Ei");

    /* Return the fitted parameters: */
    LP->f = f;
    LP->Ed = Ed;
    LP->Eb = Eb;
    LP->Ei = Ei;
  }
  
void ensure_non_negative_light(double *Evalue, char *Ename)
  {
    double light_fuzz = 1.0e-4;  /* Max abs of a negative indirect light coeff. */
    if ((*Evalue) < 0)
      { if ((*Evalue) < -light_fuzz) 
          { fprintf(stderr, "** warning: coeff %s = %+9.6e is negative\n", Ename, (*Evalue)); }
        (*Evalue) = 0;
      }
  }


  



  
void filterNormalsByLimit(r3_t*  normals, double* values,int N, double*** normals_fixed, double** values_fixed,int* N_fixed,double limite);

void filterNormalsByLimit(r3_t*  normals, double* values,int N, double*** normals_fixed, double** values_fixed,int* N_fixed,double limite){
  int new_N = 0;
  int i;
  for(i = 0;i < N; i++){
    if(values[i] < limite) new_N++;
  }
  
  double** new_normals =  (double**)malloc(sizeof(double*)*new_N);
  double* new_values = (double*)malloc(sizeof(double)*new_N);
  int count = 0;
  for(i = 0;i < N; i++){
    if(values[i] < limite){
      new_normals[count] = (double*)malloc(sizeof(double)*3);
      new_values[count] = values[i];
      //new_normals[count] = normals[i];
      new_normals[count][0] = normals[i].c[0];
      new_normals[count][1] = normals[i].c[1];
      new_normals[count][2] = normals[i].c[2];
      count++;
    }
  }
  
  fprintf(stderr,"%d values selected from %d\n",new_N,N);
  *values_fixed = new_values;
  *normals_fixed = new_normals;
  *N_fixed = new_N;
  
}



void write_light_field_parameters(char *filename, int NC, double D[], light_field_t L[]);
  /* Writes gauge albedo {D[c]}, lightsource direction {f[c]}, and intensities of direct, backplane and
    isotropic light terms {Ed[c],Eb[c],Ei[c]} to file called {filename}, for each channel {c} in {0..NC-1}. */

void write_light_field_parameters(char *filename, int NC, double D[], light_field_t L[])
  { FILE *wr = open_write(filename, TRUE);
    int c;
    for (c = 0; c < NC; c++)
      { light_field_t *Lc = &(L[c]);
        fprintf(wr, "%2d", c);
        fprintf(wr, "  %8.6lf", D[c]);
        r3_gen_print(wr, &(Lc->f), "%9.6lf", "  ", " ", "");
        fprintf(wr, "  %8.6lf %8.6lf %8.6lf", Lc->Ed, Lc->Eb, Lc->Ei);
        fprintf(wr, "\n");
      }
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
  }





options_t *parse_args(int argc, char** argv);
  /* Parses the command line arguments. */
  
options_t *parse_args(int argc, char** argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);      

    options_t *o = (options_t *)malloc(sizeof(options_t));
    
    argparser_get_keyword(pp, "-gauge");
    o->gin = parse_gauge_args(pp, TRUE);
    fprintf(stderr,"Gauge gamma is %lf\n",o->gin->gamma);
    //o->lightingModel = MODEL_LAMBERT;
    
    
    
    argparser_get_keyword(pp, "-model");
    o->lmodel = model_parse(pp,TRUE);
    o->nSteps = 3;
    o->estimateU = FALSE;
    o->useNonLinear = FALSE;
    if(o->lmodel.modelType  == COMPACT_MODEL){
      o->estimateU = argparser_keyword_present(pp, "-estimateU");
      if(o->estimateU && argparser_keyword_present_next(pp,"steps")){
	o->nStepsU = argparser_get_next_int(pp,1,100000);
      }
      else{
	o->nStepsU = 3;
      }
      o->useNonLinear = argparser_keyword_present(pp, "-useNonLinear");
      if(o->useNonLinear && argparser_keyword_present_next(pp,"steps")){
	o->nStepsNL = argparser_get_next_int(pp,1,100000);
      }
      else{
	o->nStepsNL = 3;
      }
    }
    
    if (argparser_keyword_present(pp, "-nSteps")){
      o->nSteps = argparser_get_next_int(pp, 1, 10000);
    }
    
    o->ignorePosWeight = argparser_keyword_present(pp,"-ignorePosWeight");
    o->useCosPosWeight = argparser_keyword_present(pp,"-useCosPosWeight");
    
    if (argparser_keyword_present(pp, "-gammaOut"))
      { o->gammaOut = argparser_get_next_double(pp, 0.0, 10.0); }
    else
      { o->gammaOut = 1.0; }
    
    if (argparser_keyword_present(pp, "-maskOut"))
      { o->maskOut = argparser_get_next(pp); }
    else
      { o->maskOut = NULL; }

    if (argparser_keyword_present(pp, "-maskIn"))
      { o->maskIn = argparser_get_next(pp); }
    else
      { o->maskIn = NULL; }

    if (argparser_keyword_present(pp, "-paramsOut"))
      { o->paramsOut = argparser_get_next(pp); }
    else
      { o->paramsOut = NULL; }

    o->virtual = argparser_keyword_present(pp, "-virtual");
    o->virtualUserDefined = argparser_keyword_present(pp, "-virtualUserDefined");
    if(o->virtual){
      if(o->virtualUserDefined){
	o->got = parse_gauge_args(pp, TRUE);
	o->got->margin = 12;
      }else{
	o->got = parse_gauge_args(pp, FALSE);
      }
    }

    o->errorOut  = NULL;
    if(argparser_keyword_present(pp, "-error")){
      o->errorOut = argparser_get_next(pp);
    }
    
    argparser_finish(pp);
    return o;
  }
    

       
void write_image(char *filename, float_image_t *img, double gamma);
  /* Writes an image to file {filename}, which must have
    extensions ".fni", ".ppm", or ".pgm".  The {gamma} is used
    for PNM images. */
    
void write_image(char *filename, float_image_t *img, double gamma)
  { char *ext = strrchr(filename, '.');
    demand(ext != NULL, "image name has no extension");
    if ((strcmp(ext, ".ppm") == 0) || (strcmp(ext, ".pgm") == 0))
      { float_pnm_image_write(filename, img,FALSE, gamma, 0.0, TRUE,TRUE,TRUE); }
    else if (strcmp(ext, ".fni") == 0)
      { FILE *wr = open_write(filename, TRUE);
        float_image_write(wr, img);
        if (wr == stdout) { fflush(wr); } else { fclose(wr); }
      }
    else
      { demand(FALSE, "unknown image name extension"); }
  }
  
void test_geometry_functions(ellipse_crs_t *E);
  /* Various tests. */
  
void test_geometry_functions(ellipse_crs_t *E)
  { 
    /* Test data: */
    fprintf(stderr, "Testing geometry functions...\n");
    fprintf(stderr, "E = "); ellipse_crs_print(stderr, E, "%8.6f"); fprintf(stderr, "\n");
    
    r2_t *ctr = &(E->ctr);
    double rad = E->rad;
    r2_t *str = &(E->str);
    
    /* Test {ellipse_crs_inside}: */
    r2_t u, v;
    double mstr = r2_dir(str, &u);
    if (mstr == 0) { u = (r2_t){{ 1, 0 }}; }
    v = (r2_t){{ -u.c[1], u.c[0] }};
    
    auto void test_inside(double cu, double cv, bool_t expected);
      /* Tests insideness  of {E.ctr + cu*u + cv*v}. */
      
    void test_inside(double cu, double cv, bool_t expected)
      { fprintf(stderr, "cu = %8.6f  cv = %8.6f  expected = %c\n", cu, cv, "FT"[expected]);
        r2_t p;
        r2_mix(cu*(mstr+rad), &u, cv*rad, &v, &p);
        r2_add(ctr, &p, &p);
        assert(ellipse_crs_inside(E, &p) == expected);
      }
    
    test_inside(0.000,0.000,TRUE);
    test_inside(0.999,0.000,TRUE);
    test_inside(0.000,0.999,TRUE);
    test_inside(0.707,0.707,TRUE);
    test_inside(1.001,0.000,FALSE);
    test_inside(0.000,1.001,FALSE);
    test_inside(0.708,0.708,FALSE);
    
    /* Test {pixel_position_in_gauge_image}: */
    int irad = (int)floor(rad);
    int x = (int)floor(ctr->c[0]);
    int y = (int)floor(ctr->c[1]);
    assert(irad > 5);
    assert(pixel_position_in_gauge_image(x, y, E) == +1);
    assert(pixel_position_in_gauge_image(x + irad - 2, y, E) == +1);
    assert(pixel_position_in_gauge_image(x, y + irad - 2, E) == +1);
  }

int main(int argc, char** argv)
  {
    options_t *o = parse_args(argc,argv) ;

    gauge_data_t *gin = o->gin;
    assert(! isnan(gin->E.rad));
    assert(! isnan(gin->E.ctr.c[0]));
    
    test_geometry_functions(&(gin->E));
    
    /* Get the input gauge image {img_in} and its extension {ext}: */
    float_image_t *img_in = read_gauge_image(gin->image, gin->gamma);
    
    int NC = img_in->sz[0];
    assert(NC <= MAX_NC);

    int NXin = img_in->sz[1];
    int NYin = img_in->sz[2];
    /*Bug fix for coordinates*/
   // gin->E.ctr.c[1] =   NYin - gin->E.ctr.c[1] -1;
    
    
    /* Allocate the pixel mask: */
    float_image_t *img_mk = float_image_new(1, NXin, NYin);
    float_image_t* img_msk_in = NULL;
    if(o->maskIn != NULL){
      FILE* arq = open_read(o->maskIn,TRUE);
      img_msk_in = float_image_read(arq);
      int mNX,mNY;
      mNX = img_msk_in->sz[1];
      mNY = img_msk_in->sz[2];
      assert( (mNX == NXin) && (mNY == NYin));
      int x,y;
      for(x = 0; x < NXin; x++){
	for(y = 0; y < NYin; y++){
	  double val = float_image_get_sample(img_msk_in,0,x,y);
	  if(val < 0.5) val = 0 ;
	  else val = 1;
	  float_image_set_sample(img_msk_in,0,x,y,val);
	}
      }
    }
    
    /* Extract the normal and intensity data from the input image: */
    r3_t *X;
    int NP;
    double **F; /* We have NC arrays of values. */ 
    extract_data_from_gauge_image(img_in, img_mk,img_msk_in, &(gin->E),gin->trim, &(gin->view), &X, &F, &NP);
    
    
    
    /* Write mask for input image: */
    if (o->maskOut != NULL)
      { write_image(o->maskOut, img_mk, o->gammaOut);
        //float_image_free(img_mk);
      }

    /* Allocate the lighting parameters for all channels: */

    

    /* Fit models to each channel: */
    int c;
    void* l_data[NC];
    
    fprintf(stderr,"VIEW DIR (%lf,%lf,%lf)\n",gin->view.c[0],gin->view.c[1],gin->view.c[2]);
    approx_model_t* lm = create_approx_lighting_model(o->lmodel.modelType);
   
    for (c = 0; c < NC; c++)
      {
        fprintf(stderr,"Processing Channel %d\n",c);
        /* Compute the field parameters {L[c]} for this channel: */
        //double Dc = gin->albedo[c];
        //light_field_t *Lc = &(L[c]);
	/*WE INIT THE LIGHTING DATA HERE*/
	l_data[c] = NULL;
	void* test_data = NULL;
	if(o->lmodel.modelType == COMPACT_MODEL){
	  lighting_compact_options_t* lop = o->lmodel.options;
	  test_data = lighting_compact_init_components(lop,gin->view);
	}else if(o->lmodel.modelType == HARMONIC_MODEL){
	  lighting_harmonic_options_t* lop = o->lmodel.options;
	  test_data = lighting_harmonic_init_components(lop,gin->view);
	}else if(o->lmodel.modelType == RADIALBASIS_MODEL){
	  lighting_radial_options_t* lop = o->lmodel.options;
// 	  test_data = radialbasis_init_components(o->radialResolution,o->radialSpan,gin->E.ctr,gin->E.rad,gin->E.str,gin->view);
	  double thetaMax = acos(gin->trim/gin->E.rad);
	  fprintf(stderr,"THETAMAX IS %lf\n",thetaMax);
	  test_data = lighting_radial_init_components(lop,gin->E.ctr,gin->E.rad,gin->E.str,gin->view,thetaMax);
	}else if(o->lmodel.modelType == STEREOPOLY_MODEL){
	  lighting_stereopoly_options_t* lop = o->lmodel.options;
	  test_data = lighting_stereopoly_init_components(lop,gin->view);
	}else if( o->lmodel.modelType == NEAREST_MODEL){
	  lighting_nearest_options_t* lop = o->lmodel.options;
	  test_data = lighting_nearest_init_components(lop,gin->view);
	}
	
//	lm->ls_data = l_data[c];
	
	assert(test_data != NULL);
	double* fixed_F;
	//r3_t* fixed_X;
	int fixed_N;
	double** fixed_X;
        filterNormalsByLimit(X,F[c],NP, &fixed_X,&fixed_F,&fixed_N,0.99);
       // compute_lighting_parameters(Dc, X, F[c], NP, o->backPlane, Lc);
       // compute_lighting_parameters(Dc, fixed_X, fixed_F, fixed_N, o->backPlane, Lc);
       /*HERE WE PUT OUR FITTING FUNCTION*/
       double* weights = NULL;
       double* wpos  = NULL;
       if(o->ignorePosWeight){
	 fprintf(stderr,"Ignoring precomputed Weighting for positions\n");
	 wpos = (double*)malloc(sizeof(double)*NP);
	 rn_all(NP,1,wpos);
       }else if(o->useCosPosWeight){
	 fprintf(stderr,"Using cosine Weighting for positions\n");
	 wpos = (double*)malloc(sizeof(double)*NP);
	 int ipeso;
	 for(ipeso = 0; ipeso < NP; ipeso++){
	   r3_t normal_peso = (r3_t){{fixed_X[ipeso][0],fixed_X[ipeso][1],fixed_X[ipeso][2]}};
	   double peso_cos = r3_dot(&(gin->view),&normal_peso);
 	   peso_cos = (peso_cos < 0 ? 0: peso_cos);
	    
	   wpos[ipeso] = peso_cos;
	   if(peso_cos <=0){
  	   fprintf(stderr,"(%lf,%lf,%lf) vs (%lf,%lf,%lf) - %lf\n",normal_peso.c[0],normal_peso.c[1],normal_peso.c[2],
 								    gin->view.c[0],gin->view.c[1],gin->view.c[2],
 								    peso_cos);
	   }
	 }
       }
       if(o->estimateU){
	 if(o->lmodel.modelType == COMPACT_MODEL){
	   fprintf(stderr,"COMPACT MODEL: Estimating Pointlike Light source direction");
	   light_compact_estimate_lightdir(fixed_X,fixed_F,wpos, fixed_N,test_data,o->nStepsU);
	   lm->write_parameters(stderr,test_data);
	 }else{
	   fprintf(stderr,"Light source direction estimation is not supported by your model !\n");
	   assert(FALSE);
	 }
       }
       if(o->useNonLinear){
	 double alpha = 1.5;
	 double beta = 0.7;
	 double gamma = 0.2;
	 fprintf(stderr,"Estimating Non-Linear Parameters\n");
	 plot_non_linear_goal_function(o->paramsOut, fixed_X, fixed_F, weights,wpos, fixed_N, lm, test_data);
	 approx_model_fit_non_linear_model_simplex(fixed_X, fixed_F, weights,wpos, fixed_N, lm, test_data, o->nStepsNL,alpha, beta,gamma);
       }
       fprintf(stderr,"Estimating Linear Parameters\n");
       approx_model_fit_linear_model(fixed_X, fixed_F, weights, wpos,fixed_N, lm, test_data);
       lm->write_parameters(stderr,test_data);
//        fitModelToFunction(fixed_X,fixed_F, fixed_N,lm,test_data,o->nSteps);
       l_data[c] = test_data;

      }
    if (o->paramsOut != NULL)
      { /* Write fitted parameters to disk: */
        //write_light_field_parameters(o->paramsOut, c, gin->albedo, L);
	/*HERE WE PUT THE WRITING FUNCTION*/
	FILE *wr = open_write(o->paramsOut, TRUE);
	for( c= 0 ; c < NC; c++){
	  lm->write_parameters(wr,l_data[c]);
	}
        if (wr == stdout) { fflush(wr); } else { fclose(wr); }
	
      }
     
    if(o->errorOut != NULL){
       /* Provide defaults for output gauge image: */
        gauge_data_t *ger = o->gin;
        
        
        /* Choose the image dimensions: */
        int NXot = NXin;
        int NYot = NYin;
        
       
	float_image_t *img_plt = float_image_new(NC, NXot, NYot);
        
        paint_synthetic_gauge_image(img_plt, &(ger->E), &(ger->view), l_data,lm->evaluate);
        
	/*DUMB MANNER TO DO ERROR AND FNI*/
	float_image_t* img_err = float_image_new(NC, NXot, NYot);
	int i,j,c;
	for(c = 0; c < img_plt->sz[0]; c++){
	  for(i = 0; i < img_plt->sz[1]; i++){
	    for(j = 0; j < img_plt->sz[2]; j++){
	     double val1 = float_image_get_sample(img_in,c,i,j);
	     double val2 = float_image_get_sample(img_plt,c,i,j);
	     double mskval = float_image_get_sample(img_mk,0,i,j);
	     float_image_set_sample(img_err,c,i,j,(val1-val2)*mskval);
	    }
	  }
	}
	
	FILE* file_fni_er = open_write(o->errorOut,TRUE);
	float_image_write(file_fni_er,img_err);
	fclose(file_fni_er);

	float_image_free(img_plt);
	float_image_free(img_err);
    }
     
    if (o->virtual)
      { 
        /* Provide defaults for output gauge image: */
        gauge_data_t *got = o->got;
        assert(got->image != NULL);

	int NXot;
        int NYot;
	
        if(! o->virtualUserDefined){
	  assert(got->magnify >= 1);
	  int mag = got->magnify;
	  
	  assert(got->margin >= 0);
	  int mrg = got->margin;
	  
	  
	  assert(isnan(got->E.rad));
	  got->E.rad = mag*gin->E.rad;
	  
	  assert(isnan(got->E.ctr.c[0]));
	  got->E.ctr.c[0] = mag*gin->E.ctr.c[0] + mrg;
	  got->E.ctr.c[1] = mag*gin->E.ctr.c[1] + mrg;
	  
	  assert(isnan(got->E.str.c[0]));
	  got->E.str.c[0] = mag*gin->E.str.c[0];
	  got->E.str.c[1] = mag*gin->E.str.c[1];
	  
	  assert(isnan(got->view.c[0]));
	  got->view = gin->view;
	  
	  /* Choose the image dimensions: */
	  NXot = mag * NXin + 2*mrg;
	  NYot = mag * NYin + 2*mrg;
        
	}else{
	  int mrg = got->margin;
	    
	  NXot = got->E.rad*2 + 2*(mrg);
	  NYot = got->E.rad*2 + 2*(mrg);
	
	  
	  
	    
	  got->E.ctr.c[0]+= mrg;
	  got->E.ctr.c[1]+= mrg;
  
	  
	}
	
        float_image_t *img_ot = float_image_new(NC, NXot, NYot);
        
        paint_synthetic_gauge_image(img_ot, &(got->E), &(got->view), l_data,lm->evaluate);
        
	char* filename_ppm_ot = NULL;
	char *filename_ppm_ot = jsprintf("%s.ppm",got->image);
        /* Write virtual gauge image: */
        write_image(filename_ppm_ot, img_ot, o->gammaOut);
	
	
	
	/*Write FNI too*/
	char* filename_fni_ot = NULL;
	char *filename_fni_ot = jsprintf("%s.fni",got->image);
	FILE* file_fni_ot = open_write(filename_fni_ot,TRUE);
	float_image_write(file_fni_ot,img_ot);
	fclose(file_fni_ot);
	
	
	
        float_image_free(img_ot);
      }
    return 0;
  }
