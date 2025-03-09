#define PROG_NAME "test_sve_img_match"
#define PROG_DESC "tests {sve_minn.h} on an image matching problem"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-16 20:23:55 by stolfi */

#define test_sve_img_match_C_COPYRIGHT "Copyright © 2009 by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <values.h>

#include <bool.h>
#include <sign.h>
#include <float_image.h>
#include <float_image_mask.h>
#include <float_image_transform.h>
#include <argparser.h>
#include <r3x3.h>
#include <r3.h>
#include <r2_extra.h>
#include <rn.h>
#include <rmxn.h>
#include <affirm.h>
#include <ix_reduce.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <jsfile.h>
#include <js.h>
#include <vec.h> 
#include <minn_plot.h>

#define PROG_HELP = \
  "  " PROG_NAME " \\\n" \
  "    -img {IMG_NAME}.fni \\\n" \
  "    -obj {OBJ_NAME}.fni \\\n" \
  "    [ -msk {MSK_NAME}.fni ] \\\n" \
  "    [ -shrink {SHRINK} ] \\\n" \
  "      \\\n" \
  "    [ -center {CTRX} {CTRY} ] \\\n" \
  "    [ -angle {ANGLE} ] \\\n" \
  "    [ -scale {SCALE} ] \\\n" \
  "      \\\n" \
  "    [ -rotate {ROTF} ] \\\n" \
  "    [ -translate {TRF} ] \\\n" \
  "    [ -resize {RSZF} ] \\\n" \
  "    [ -shear {SHRF} ] \\\n" \
  "    [ -perspect {PSPF} ] \\\n" \
  "    -outPrefix {OPREF} "

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests {sve_minn.h} on an image matching problem.  The" \
  " parameter is a 3x3 projective transformation matrix {M}.  The" \
  " goal function is the quadratic mismatch between a given" \
  " image {obj}, windowed by an optional mask {msk} and" \
  " transformed by {M}, and a second image {img}.\n" \
  "\n" \
  "  The three images {img}, {obj}, and {msk} are read" \
  " from the files \"{IMG_NAME}.fni\", \"{OBJ_NAME}.fni\"," \
  " and \"{MSK_NAME}.fni\", respectively.  If the mask" \
  " is not specified, or is \"NONE\", uses a circular" \
  " biquadratic window.  If the \"-shrink\" option is given, all three images" \
  " are reduced by {1/SHRINK} before running the" \
  " test.  The {CTRX,CTRY} parameters are adjusted accordingly.\n" \
  "\n" \
  "\n" \
  "  The program initially looks for {obj} rotated" \
  " by {ANGLE}, scaled by {SCALE}, and with its" \
  " center translated to {(CTRX,CTRY)} relative" \
  " to the lower left corner of {img}.\n" \
  "\n" \
  "  This initial mapping is" \
  " then subjected to minor adjustments, described my the boolean" \
  " arguments \"-rotate\", \"-translate\", \"-resize\", \"-shear\", and" \
  " \"-perspect\".\n" \
  "\n" \
  "  The program writes files \"{OPREF}-{ITER}-img.fni\" and" \
  " \"{OPREF}-{ITER}-obj.fni\", where {ITER} is \"ini\" or" \
  " \"opt\", showing the matched sub-image" \
  " of {img} and the transformed {obj}, in situ.\n" \
  "\n" \
  "OPTIONS\n" \
  "  See above.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  Hamlet(10), Venice(11), Luxor(11), Amazon(12).\n" \
  "\n" \
  "AUTHOR\n" \
  "  2009-01-05 Created by Jorge Stolfi, IC-UNICAMP.\n" \
  "  2012-01-25 J. Stolfi: Substantially modified the parameter\n" \
  "             encoding, added image output and goal function plots.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_sve_img_match_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS


#include <sve_minn.h>

/* GENERAL PARAMETERS */

typedef struct adjust_options_t
  { /* Modifications allowed in map: */
    bool_t rotate;    /* If TRUE, adjusts rotation angle. */
    bool_t translate; /* If TRUE, adjusts translation. */
    bool_t resize;    /* If TRUE, adjusts magnification. */
    bool_t shear;     /* If TRUE, adjusts area-preserving linear deformation. */
    bool_t perspect;  /* If TRUE, allows perspective deformation. */
  } adjust_options_t;
  /* Adjustments allowed in matrix. */

typedef struct options_t
  { char *img;        /* Target image (where to search). */
    char *obj;        /* Object image (to be sought). */
    char *msk;        /* Mask image for {obj}. */
    char *outPrefix;  /* Prefix for output file names.*/
    int32_t shrink;       /* Image reduction factor, or 1 if none. */
    /* Initial guess parameters of {obj} relative to {img}: */
    r2_t center;      /* Center position, or {(NAN,NAN)} to mean center of {img}. */
    double angle;     /* Rotation angle (degrees ccw). */
    double scale;     /* Scaling factor. */
    /* Modifications allowed in map: */
    adjust_options_t adop;
  } options_t;
  /* Command-line parameters. */

#define MAXNQ 100
  /* Max number of . */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

options_t *get_options(int32_t argc, char **argv);
  /* Parses the command-line options. */

float_image_t *read_input_image(char *img_name, int32_t shrink, int32_t *NCP, int32_t *NXP, int32_t *NYP);
  /* If {img_name} is NULL, returns NULL. 
  
    Otherwise reads an FNI image from file "{img_name}", and shrinks it
    by {1/shrink}. Then, if {*NCP} is negative, sets {*NCP} to the
    number of channels {NC} in the image; otherwise requires {*NCP ==
    NC}. Ditto for {*NXP} and {NYP}. */

float_image_t *shrink_image(float_image_t *A, int32_t shrink);
  /* Returns a copy of image {A} shrunk by the factor {1/shrink}. */

void initial_matrices
  ( int32_t NXO, 
    int32_t NYO, 
    double mag, 
    double ang, 
    r2_t ctr, 
    int32_t NXI, 
    int32_t NYI,
    r3x3_t *P,
    r3x3_t *Q
  );
  /* Creates the initial matrices {P} and {Q} from the size
    {NXO,NYO} of {obj}, the scale factor {mag}, the rotation angle {ang}, the
    desired center {ctr} in {img}, and the size {NXI,NYI} of {img}. If
    {ctr} is {(NAN,NAN)}, uses the center of {img}.
    
    Let {r} be the circumradius of {obj}'s domain. The matrix {P}
    maps the center of {obj} to the origin, then scales it by {mag/r}
    and rotates it by {ang}. The matrix {Q} translates the origin to
    the point {ctr}. The working map will be {M = P*A*Q}
    where {A} is a perturbation of the identity. */

double mismatch(float_image_t *img, float_image_t *obj, float_image_t *msk, r3x3_t *M);
  /* Relative mismatch between the image {obj}, weighted by the mask {msk}
    and transformed by {M}, and the image {img}.
    
    The mismatch is defined as {(1/A)\int\int msk(p)*(obj(p) -
    img(M(p)))^2 dp}, where the integration is over the domain of
    {obj} (which must be the domain of {msk}), and {A} is the area of
    that domain. */

int32_t num_parameters(adjust_options_t *adop);
  /* The number of parameters neeeded to specify a particular matrix in 
    the class of matrices allowed by {adop}. */

void param_to_matrix
  ( r3x3_t *P, 
    adjust_options_t *adop, 
    int32_t n, 
    double x[], 
    r3x3_t *Q,
    r3x3_t *M, 
    bool_t verbose
  );
  /* Sets {*M} to {P*A*Q} where {A} is the identity with the 
    perturbations allowed by {adop} and parametrized by {x[0..n-1]}.  
    If {verbose} is true, prints {A} and {M}. */

void find_best_matrix
  ( float_image_t *img, 
    float_image_t *obj, 
    float_image_t *msk, 
    r3x3_t *P,
    r3x3_t *Q,
    r3x3_t *M,
    adjust_options_t *adop,
    char *outPrefix
  );
  /* Stores into {M} the matrix that minimizes {mismatch(img,obj,msk,M)},
    The matrix will be {M = P*A*Q} where {A} is a perturbation
    of the identity of the class allowed by {adop}. */

void write_matched_object_images
  ( char *outPrefix,
    char *version,
    float_image_t *img, 
    float_image_t *obj, 
    float_image_t *msk, 
    r3x3_t *M,
    float bgrf[],
    bool_t object
  );
  /* Writes three FNI images that illustrates the matching of image {obj},
    masked by {msk} and transformed by {M}, against image {img}. 
    
    The images are written to files "{outPrefix}-{version}-{tag}.fni"
    where {tag} is either "img", "obj", or "msk".
    
    Each output image has the same size as {img}. The first two images
    (whose {tag}s are "img" and {obj}) consist of a uniform background
    image {bgr(x,y)} blended with a second image {val(x,y)}, which is
    either {img} itself, or {obj} mapped bt {M}, respectively. The
    images {bgr} and {val} are blended with weights {1-msk'(x,y)} and
    {msk'(x,y)} respectively, where {msk'} is the mask image {msk}
    transformed by {M}. The uniform background image is filled with the
    pixel {bgrf[0..NC-1]}.
    
    The third image ({tag = "msk"}) consists of the mapped mask {msk'} padded
    with a black background to the size of {img}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { options_t *o = get_options(argc, argv);
  
    /* Read the target image {img}: */
    int32_t NC = -1; /* Number of channels in {img,obj}. */
    int32_t NXI = -1, NYI = -1; /* Size of {img}. */
    float_image_t *img = read_input_image(o->img, o->shrink, &NC, &NXI, &NYI);
    
    /* Read the object image {obj}: */
    int32_t NXO = -1, NYO = -1; /* Size of {obj} and {msk}. */
    float_image_t *obj = read_input_image(o->obj, o->shrink, &NC, &NXO, &NYO);

    /* Read or manufacture the mask image {msk} for {obj}: */
    int32_t NCM = 1;  /* Number of channels in {msk}. */
    float_image_t *msk = read_input_image(o->msk, o->shrink, &NCM, &NXO, &NYO);
    if (msk == NULL) 
      { msk = float_image_new(1, NXO, NYO); 
        float_image_mask_window(msk, 0, 2, TRUE);
      }

    /* Adjust the parameters for the reduction factor: */
    if (o->shrink != 1) { r2_scale(1/((double)o->shrink), &(o->center), &(o->center)); }
    
    float bgrf[NC]; /* Background pixel value for output images. */
    for (uint32_t c = 0;  c < NC; c++) { bgrf[c] = 0.5; }
    
    /* Compute the left and right matrices {Q,P}: */
    r3x3_t P, Q;
    initial_matrices(NXO, NYO, o->scale, o->angle, o->center, NXI, NYI, &P, &Q);
    
    /* Compute the initial matrix and write its effect: */
    r3x3_t M;
    r3x3_mul(&P, &Q, &M);
    write_matched_object_images(o->outPrefix, "ini", img, obj, msk, &M, bgrf, TRUE);
    
    /* Find and print the optimum transform: */
    find_best_matrix(img, obj, msk, &P, &Q, &M, &(o->adop), o->outPrefix);
    
    /* Write the diagnostic images: */
    write_matched_object_images(o->outPrefix, "opt", img, obj, msk, &M, bgrf, TRUE);
    
    fclose(stdout);
    return (0);
  }

float_image_t *read_input_image(char *img_name, int32_t shrink, int32_t *NCP, int32_t *NXP, int32_t *NYP)
  {
    if ((img_name == NULL) || (strlen(img_name) == 0))
      { return NULL; }
    else
      { FILE *rd = open_read(img_name, TRUE);
        float_image_t *img = float_image_read(rd);
        fclose(rd);
        if (shrink != 1)
          { float_image_t *red = shrink_image(img, shrink);
            float_image_free(img);
            img = red;
          }
        /* Get/check dimensions: */
        int32_t NC = (int32_t)img->sz[0];
        int32_t NX = (int32_t)img->sz[1];
        int32_t NY = (int32_t)img->sz[2];
        if ((*NCP) == -1) { (*NCP) = NC; } else { demand((*NCP) == NC, "invalid channel count"); }
        if ((*NXP) == -1) { (*NXP) = NX; } else { demand((*NXP) == NX, "invalid column count"); }
        if ((*NYP) == -1) { (*NYP) = NY; } else { demand((*NYP) == NY, "invalid row count"); }
        return img;
      }
  }

float_image_t *shrink_image(float_image_t *A, int32_t shrink)
  { 
    demand(shrink >= 1, "invalid reduction factor");

    /* Generate the 1D weight mask: */
    int32_t nw = 3*shrink;
    double wt[nw];
    wt_table_hann_fill(nw, 0.0, wt, NULL);
    wt_table_normalize_sum(nw, wt);
    
    /* Get the image dimensions: */
    int32_t NC  = (int32_t)A->sz[0];
    int32_t NXA = (int32_t)A->sz[1];
    int32_t NYA = (int32_t)A->sz[2];

    /* Create the output image {R}: */
    int32_t NXR = (NXA + shrink - 1)/shrink;
    int32_t NYR = (NYA + shrink - 1)/shrink;
    float_image_t *R = float_image_new(NC, NXR, NYR);
    
    /* Fill the pixels of {R}: */
    for (uint32_t yR = 0;  yR < NYR; yR++)
      { for (uint32_t xR = 0;  xR < NXR; xR++)
          { /* Accumulate the weighted pixel sum over the window: */
            double sum_w = 0;    /* Sum of weights. */
            double sum_w_v[NC];  /* Sum of weighted pixels. */
            for (uint32_t c = 0;  c < NC; c++) { sum_w_v[c] = 0; }
            for (uint32_t yD = 0;  yD < nw; yD++)
              { int32_t yA = shrink*yR - shrink + yD;
                double wty = ((yA < 0) || (yA >= NYA) ? 0.0 : wt[yD]);
                for (uint32_t xD = 0;  xD < nw; xD++)
                  { int32_t xA = shrink*xR - shrink + xD;
                    double wtx = ((xA < 0) || (xA >= NXA) ? 0.0 : wt[xD]);
                    double w = wtx*wty;
                    if (w > 0)
                      { /* Multiply by the mask weight, if any: */
                        sum_w += w; 
                        for (uint32_t c = 0;  c < NC; c++)
                          { double v = float_image_get_sample(A, c, xA, yA);
                            sum_w_v[c] += w*v;
                          }
                      }
                  }
              }
            /* Store the weighted average (maybe NAN) in the {R} pixel: */
            if (sum_w == 0) { sum_w = 1; }
            for (uint32_t c = 0;  c < NC; c++)
              { double v = sum_w_v[c]/sum_w; 
                float_image_set_sample(R, c, xR, yR, (float)v);
              }
          }
      }
    return R;
  }

void initial_matrices
  ( int32_t NXO, 
    int32_t NYO, 
    double mag, 
    double ang, 
    r2_t ctr, 
    int32_t NXI, 
    int32_t NYI,
    r3x3_t *P,
    r3x3_t *Q
  )
  {
    /* Set {P} to translation of center of {obj} to origin: */
    r2_t obj_ctr = (r2_t){{ 0.5*NXO, 0.5*NYO }}; /* Center of {obj}. */
    double rO = r2_norm(&obj_ctr);
    r3x3_ident(P);
    P->c[0][1] = - obj_ctr.c[0];
    P->c[0][2] = - obj_ctr.c[1];
    
    /* Append a scaling that reduces {obj} to radius {mag}: */
    r3x3_t S; r3x3_ident(&S);
    S.c[1][1] = mag/rO;
    S.c[2][2] = mag/rO;
    r3x3_mul(P, &S, P);
    
    /* Append a rotation by {ang} degrees ccw: */
    double t = ang*M_PI/180; /* Angle ccw in radians. */
    double c = cos(t), s = sin(t);
    r3x3_t R; r3x3_ident(&R);
    R.c[1][1] = +c;
    R.c[1][2] = +s;
    R.c[2][1] = -s;
    R.c[2][2] = +c;
    r3x3_mul(P, &R, P);
    fprintf(stderr, "  P =\n");
    r3x3_gen_print(stderr, P, "%12.8f", "", "", "", "  ", " ", "\n");   
    
    /* Set {Q} to a matrix that rescales {obj} to radius {rO*mag}: */
    r3x3_ident(Q);
    Q->c[1][1] = rO;
    Q->c[2][2] = rO;

    /* Append to {Q} a translation of origin to {ctr} (or center of {img}): */
    if (isnan(ctr.c[0])) { ctr.c[0] = 0.5*NXI; }
    if (isnan(ctr.c[1])) { ctr.c[1] = 0.5*NYI; }
    Q->c[0][1] = + ctr.c[0];
    Q->c[0][2] = + ctr.c[1];
    fprintf(stderr, "  Q =\n");
    r3x3_gen_print(stderr, Q, "%12.8f", "", "", "", "  ", " ", "\n");   
  }    

void find_best_matrix
  ( float_image_t *img, 
    float_image_t *obj, 
    float_image_t *msk, 
    r3x3_t *P,
    r3x3_t *Q,
    r3x3_t *M,
    adjust_options_t *adop,
    char *outPrefix
  )
  { 
    int32_t NXO = (int32_t)obj->sz[1];
    int32_t NYO = (int32_t)obj->sz[2];
    
    int32_t nx = num_parameters(adop); /* Number of degrees of freedom in matrix. */
    
    bool_t plot_goal = TRUE;
    
    auto double sve_goal(int32_t n, double x[]); 
      /* The goal function for uptimization.  Also sets {M} to the corresp. matrix. */
      
    double xPrev[nx]; /* Parameter vector in previous call of {sve_OK} function. */
    r3x3_t MPrev;     /* Matrix in previous call of  {sve_OK} function. */
    int32_t nok = 0;      /* Counts iterations (actually, calls to {sve_OK}). */
    bool_t sve_debug = TRUE;
    bool_t sve_debug_probes = FALSE;
    
    auto bool_t sve_OK(int32_t iter, int32_t n, double x[], double Fx, double dist, double step, double radius); 
      /* Acceptance criterion function. */
      
    double xx[nx];     /* Initial guess and final solution. */

    /* Start with the null paramenter vector, which should be the identity matrix: */
    { int32_t ix; for (ix = 0; ix < nx; ix++) { xx[ix] = 0.0; } }

    fprintf(stderr, "initial matrix:\n");
    param_to_matrix(P, adop, nx, xx, Q, M, TRUE);
    
    double Fx = sve_goal(nx, xx);
    fprintf(stderr, "\n");
    fprintf(stderr, "initial goal function = %16.10f\n", Fx);

    if (nx > 0)
      { /* Optimize iteratively: */
        double *ctr = NULL;
        double dMax = 0.250;
        bool_t dBox = FALSE;
        double rMin = 1.0/hypot(0.5*NXO,0.5*NYO);
        double rMax = 0.5*dMax;
        double rIni = 0.25*dMax;
        double minStep = 0.0001*rMin;
        sign_t dir = -1;
        int32_t maxIters = 100;
        
        if (plot_goal)
          { /* Choose the plot directions: */
            int32_t nu = (nx < 6 ? nx : 6); /* Number of directions. */
            double U[nu*nx]; /* The rows are the directions. */
            rmxn_throw_directions(nu, nx, U);
            /* Plot goal function along those directions: */
            double Rad = rIni;
            double Step = Rad/30;
            char *fname = jsprintf("%s-f2-plot.txt", outPrefix);
            FILE *wr = open_write(fname, TRUE);
            for (uint32_t ku = 0;  ku < nu; ku++)
              { double *u = &(U[ku*nx]);
                minn_plot_1D_gnuplot(wr, nx, xx, u, Rad, Step, sve_goal);
              }
            fclose(wr);
            free(fname);
          }

        sve_minn_iterate
          ( nx, &sve_goal, &sve_OK, NULL,
            xx, &Fx, dir, 
            ctr, dMax, dBox, rIni, rMin, rMax, 
            minStep, maxIters, sve_debug, sve_debug_probes
          );
        fprintf(stderr, "\n");
      }

    /* Print final matrix: */
    fprintf(stderr, "final matrix:\n");
    param_to_matrix(P, adop, nx, xx, Q, M, TRUE);

    /* Ealuate (and print) the final mismatch: */
    double FxN = sve_goal(nx, xx); 
    fprintf(stderr, "\n");
    fprintf(stderr, "final goal function = %16.10f %16.10f\n", Fx, FxN);
    demand(Fx == FxN, "inconsistent function value on return");
    return;
    
    double sve_goal(int32_t n, double x[])
      { assert(n == nx);
        r3x3_t MM;
        param_to_matrix(P, adop, n, x, Q, &MM, FALSE);
        double Fx = mismatch(img, obj, msk, &MM);
        fprintf(stderr, "sve_goal");
        return Fx;
      }
      
    bool_t sve_OK(int32_t iter, int32_t n, double x[], double Fx, double dist, double step, double radius)
      { assert(n == nx);
        fprintf(stderr, "\n");
        fprintf(stderr, "iteration %d matrix:\n", nok);
        param_to_matrix(P, adop, n, x, Q, M, TRUE);
        if (nok > 0)
          { double d = rn_dist(n, xPrev, x);
            fprintf(stderr, "  change in x = %16.10f\n", d);
            r3x3_t D; r3x3_inv(&MPrev, &D); r3x3_mul(&D, M, &D);
            fprintf(stderr, "  change in M = \n");
            r3x3_gen_print(stderr, &D, "%12.8f", "", "", "", "  ", " ", "\n");
          }
        fprintf(stderr, "mismatch = %16.10f\n", Fx);
        fprintf(stderr, "\n");
        /* Save guess in {xPrev,MPrev} for next call: */
        rn_copy(nx, x, xPrev); MPrev = (*M);
        nok++;
        return FALSE;
      }
    
  }
 
void write_matched_object_images
  ( char *outPrefix,
    char *version,
    float_image_t *img, 
    float_image_t *obj, 
    float_image_t *msk, 
    r3x3_t *M,
    float bgrf[],
    bool_t object
  )
  {
    /* Get/check the input image dimensions: */
    int32_t NC  = (int32_t)img->sz[0];                    /* Number of channels of {img,obj}. */
    int32_t NXI = (int32_t)img->sz[1], NYI = (int32_t)img->sz[2]; /* Size of {img}. */
    int32_t NXO = (int32_t)obj->sz[1], NYO = (int32_t)obj->sz[2]; /* Size of {obj} and {msk}. */

    assert(msk != NULL);
    assert(obj->sz[0] == NC);
    assert(msk->sz[0] == 1);
    assert(msk->sz[1] == NXO); 
    assert(msk->sz[2] == NYO);
    
    /* Compute the inverse {R} of matrix {M}: */
    r3x3_t R; r3x3_inv(M, &R);

    auto void map_img_to_obj(r2_t *p, r2x2_t *J);
      /* Maps a point {p} of {img} domain to the corresponding point {q} in {obj} and {msk} domain. */ 
      
    /* Output canvases: */
    float_image_t *img_out = float_image_new(NC, NXI, NYI);
    float_image_t *obj_out = float_image_new(NC, NXI, NYI);
    float_image_t *msk_out = float_image_new(1,  NXI, NYI);
    
    int32_t order = 1;                             /* C1 bicubic interpolation. */
    ix_reduce_mode_t red = ix_reduce_mode_SINGLE;  /* Index reduction method. */
    float imgf[NC];                            /* Pixel from {img}, unmapped. */
    float objf[NC];                            /* Pixel from {obj} mapped by {M}. */
    float mskf;                                /* Sample from {msk} mapped by {M}. */
    int32_t ic, ix, iy;
    for (iy = 0; iy < NYI; iy++)
      for (ix = 0; ix < NXI; ix++)
        { /* Get the (mapped) mask value {mskf[0]} for the output pixel {ix,iy}: */
          bool_t debug_pix = FALSE;
          float_image_transform_get_pixel
            ( msk, red, ix, iy, &map_img_to_obj, 0.0, TRUE, order, &(mskf), debug_pix );
          if (isnan(mskf) || (mskf == 0))
            { float_image_set_pixel(img_out, ix, iy, bgrf); 
              float_image_set_pixel(obj_out, ix, iy, bgrf); 
              float_image_set_sample(msk_out, 0, ix, iy, 0.0); 
            }
          else
            { if ((mskf < 0.0) || (mskf > 1.0))
                { fprintf(stderr, "pixel [%d,%d] mskf = %24.15e\n", ix, iy, mskf);
                  mskf = (float)fmin(1.0, fmax(0.0, mskf));
                }
               /* Get the source pixel value {valf}: */
              ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
              float_image_transform_get_pixel
                ( obj, red, ix, iy, &map_img_to_obj, NAN, TRUE, order, objf, debug_pix );
              float_image_get_pixel(img, ix, iy, imgf);
              /* Blend {valf} with {bgvf} with weight {mskf}: */
              for (ic = 0; ic < NC; ic++) 
                { objf[ic] = mskf*objf[ic] + (1-mskf)*bgrf[ic];
                  imgf[ic] = mskf*imgf[ic] + (1-mskf)*bgrf[ic];
                }
              /* Store it into the output images: */
              float_image_set_pixel(img_out, ix, iy, imgf); 
              float_image_set_pixel(obj_out, ix, iy, objf); 
              float_image_set_sample(msk_out, 0, ix, iy, mskf); 
            }
        }
        
    /* Write the images: */
    for (uint32_t it = 0;  it < 3; it++)
      { /* Write image number {it}: */
        char *tag = ((char*[]){ "img", "obj", "msk" })[it];
        float_image_t *out = ((float_image_t *[]){ img_out, obj_out, msk_out })[it];
        
        char *fname = jsprintf("%s-%s-%s.fni", outPrefix, version, tag);
        FILE *wr = open_write(fname, TRUE);
        float_image_write(wr, out);
        fclose(wr);
        free(fname);
        float_image_free(out);
      }
    
    return;

    /* INTERNAL IMPLEMENTATIONS */
    
    void map_img_to_obj(r2_t *p, r2x2_t *J)
      { r2x2_ident(J);
        r2_map_projective(p, &R, p, J);
      }
  }

int32_t num_parameters(adjust_options_t *adop)
  { int32_t n = 0;

    if (adop->translate) 
      { n += 2;}

    if (adop->shear)
      { n += 3; }
    else if (adop->rotate)
      { n += 1; }

    if (adop->resize)
      { n += 1; }
      
    if (adop->perspect)
      { n += 2; }
      
    return n;
  }

void param_to_matrix
  ( r3x3_t *P, 
    adjust_options_t *adop, 
    int32_t n, 
    double x[], 
    r3x3_t *Q,
    r3x3_t *M, 
    bool_t verbose
  )
  { r3x3_t A; r3x3_ident(&A);
    int32_t k = 0;
    
    if (adop->translate) 
      { /* Translation terms: */
        A.c[0][1] = x[k]; k++;
        A.c[0][2] = x[k]; k++;
      }
      
    /* Shear and rotate are mutually exclusive: */
    if (adop->shear)
      { /* Set the linear submatrix to a general area-preserving linear map. */
        /* Pick one unit vector {(ux,uy)} with random direction given by angle parameter: */
        double t = x[k]; k++;
        double ux = cos(t), uy = sin(t);
        /* Get the orthogonal CCW vector {(vx,vy)}: */
        double vx = -uy, vy = ux;
        /* Apply shear according to parameter: */
        double s = x[k]; k++;
        vx += s*ux; vy += s*vy;
        /* Modify lengths according to parameter: */
        double f = x[k]; k++;
        double fe = exp(f);
        ux *= fe;
        uy *= fe;
        vx /= fe;
        vy /= fe;
        /* Pack in matrix: */
        A.c[1][1] = ux;
        A.c[1][2] = uy;
        A.c[2][1] = vx;
        A.c[2][2] = vy;
      } 
    else if (adop->shear)
      { /* Set the linear submatrix to a rotation matrix. */
        /* Get the angle parameter: */
        double t = x[k]; k++;
        /* Fill the submatrix: */
        double ct = cos(t), st = sin(t);
        A.c[1][1] = +ct;
        A.c[1][2] = +st;
        A.c[2][1] = -st;
        A.c[2][2] = +ct;
      }
        
    if (adop->resize)
      { /* Apply a uniform scaling before translation and rotation/shear: */
        double s = x[k]; k++;
        double mag = 1 + s;
        A.c[1][1] *= mag;
        A.c[1][2] *= mag;
        A.c[2][1] *= mag;
        A.c[2][2] *= mag;
      }
    
    if (adop->perspect) 
      { /* Perspective distortion coeffs: */
        A.c[1][0] = x[k]; k++;
        A.c[2][0] = x[k]; k++;
      }
      
    assert(k == n);
    
    if (verbose) 
      { fprintf(stderr, "  A = \n");
        r3x3_gen_print(stderr, &A, "%12.8f", "", "", "", "  ", " ", "\n");
      }

    /* Combine with {P,Q}: */
    r3x3_mul(P, &A, M);
    r3x3_mul(M, Q, M);
    if (verbose) 
      { fprintf(stderr, "  M = \n");
        r3x3_gen_print(stderr, M, "%12.8f", "", "", "", "  ", " ", "\n");
      }
  }

double mismatch
  ( float_image_t *img, 
    float_image_t *obj, 
    float_image_t *msk, 
    r3x3_t *M
  )
  {
    /* Get/check the input image dimensions: */
    int32_t NC  = (int32_t)img->sz[0];  /* Number of channels of {img,obj}. */
    int32_t NXO = (int32_t)obj->sz[1], NYO = (int32_t)obj->sz[2]; /* Size of {obj} and {msk}. */

    assert(msk != NULL);
    assert(obj->sz[0] == NC);
    assert(msk->sz[0] == 1);
    assert(msk->sz[1] == NXO); 
    assert(msk->sz[2] == NYO);

    auto void map_obj_to_img(r2_t *p, r2x2_t *J);
      /* Maps a point {p} of {obj} and {msk} domain to the corresponding point in {img} domain. */ 
      
    int32_t order = 1; /* C1 bicubic interpolation. */
    ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
    float undef = 0.5;
    bool_t avg = TRUE;

    float objf[NC];
    float imgf[NC];
    double sumw = 0;
    double sumwd2 = 0;
    for (uint32_t iy = 0;  iy < NYO; iy++)
      for (uint32_t ix = 0;  ix < NXO; ix++)
        { /* Get pixel {ix,iy} of {obj} in {objf[0..NC-1]} and of {msk} in {mskf[0]}: */
          float_image_get_pixel(obj, ix, iy, objf);
          float mskf = float_image_get_sample(msk, 0, ix, iy);
          if ((! isnan(mskf)) && (mskf > 0))
            { /* Get the corresponding pixel from {img}: */
              bool_t debug_pix = FALSE;
              float_image_transform_get_pixel
                ( img, red, ix, iy, &map_obj_to_img, undef, avg, order, imgf, debug_pix );
              /* Accumulate square diff: */
              for (uint32_t ic = 0;  ic < NC; ic++)
                { if (! isnan(imgf[ic]))
                    { double d = ((double)imgf[ic]) - ((double)objf[ic]);
                      sumw += mskf;
                      sumwd2 += mskf*d*d;
                    }
                }
            }
        }
    return sumwd2/sumw;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void map_obj_to_img(r2_t *p, r2x2_t *J)
      { r2x2_ident(J);
        r2_map_projective(p, M, p, J);
      }
  }    

options_t *get_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_get_keyword(pp, "-img");
    o->img = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-obj");
    o->obj = argparser_get_next_non_keyword(pp);
    
    if (argparser_keyword_present(pp, "-msk"))
      { o->msk = argparser_get_next_non_keyword(pp);
        if (strcmp(o->msk, "NONE") == 0) { o->msk = NULL; } }
    else
      { o->msk = NULL; }
    
    if (argparser_keyword_present(pp, "-center"))
      { o->center.c[0] = argparser_get_next_double(pp, -100000, 100000);
        o->center.c[1] = argparser_get_next_double(pp, -100000, 100000);
      }
    else
      { o->center = (r2_t){{ NAN, NAN }}; }
    
    if (argparser_keyword_present(pp, "-angle"))
      { o->angle = argparser_get_next_double(pp, -3600, +3600); }
    else
      { o->angle = 0.0; }
    
    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, 0.001, 1000.0); }
    else
      { o->scale = 1.0; }
    
    if (argparser_keyword_present(pp, "-shrink"))
      { o->shrink = (int32_t)argparser_get_next_int(pp, 1, 100000); }
    else
      { o->shrink = 1; }
    
    o->adop.rotate = argparser_keyword_present(pp, "-rotate");
    o->adop.translate = argparser_keyword_present(pp, "-translate");
    o->adop.resize = argparser_keyword_present(pp, "-resize");
    o->adop.shear = argparser_keyword_present(pp, "-shear");
    o->adop.perspect = argparser_keyword_present(pp, "-perspect");

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    argparser_finish(pp);
    return o;
  }
