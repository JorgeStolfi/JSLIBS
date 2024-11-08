#define PROG_NAME "test_mfok_diff_ops"
#define PROG_DESC "Analyzes relation between image sharpness and local quadratic operators."
#define PROG_VERS "1.0"

/* Last edited on 2024-10-26 09:23:46 by stolfi */ 
/* Created on 2023-01-24 by J. Stolfi, UNICAMP */

#define test_mfok_diff_ops_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -inDir {inDir} \\\n" \
  "    -imageSize {NX} {NY} \\\n" \
  "    { -frame {zFoc[ki]} {zDep[ki]} }.. \\\n" \
  "    -basisType {basisType} \\\n" \
  "    -windowType {windowType} \\\n" \
  "    -termSetFile {termSetFile} \\\n" \
  "    -sampleNoise {sampleNoise} \\\n" \
  "    -outDir {outDir}" \

#define PROG_INFO \
  "SYNOPSIS\n" \
  "  Tries to indentify an optimum local focus sharpness formula from" \
  " images where that information is known.\n" \
  "\n" \
  "  Reads a stack {istack} of {NI} synthetic input frames (images)" \
  " {istack,frame[0..NI-1]}  with focus blur.  For each frame, enumerates" \
  " the {NW} by {NW} windows centered at each pixel.  For each window" \
  " position, maps the {NS = NW*NW} image samples in the window to" \
  " coefficients {coeff[0..NB-1]} of a  specified local operator" \
  " basis with {NB} elements.  Writes {NI*NB} images showing those" \
  " coefficients.  Then computes {NT} specified quadratic" \
  " terms {term[0..NT-1]} from those coefficients. Writes {NI*NT} images showing those" \
  " terms.  Then writes this data as a text file for plotting and regression.\n" \
  "\n" \
  "INPUTS\n" \
  "\n" \
  "  INPUT FRAMES\n" \
  "    Each input frame {ifr = istack.frame[ki]}, for {ki} in {0..NI-1} is actually" \
  " a quadruple of images {ifr.sVal}, {ifr.hAvg}, {ifr.hDev}, {ifr.shrp}, where\n" \
  "\n" \
  "      {ifr.sVal} is sythetic view of a 3D scene  with simulated focus blurring.\n" \
  "      {ifr.hAvg} specifies the average {Z}-coord of the scene within each pixel.\n" \
  "      {ifr.hDev} specifies the deviation of that {Z} coordinate within each pixel.\n" \
  "      {ifr.shrp} specifies the actual focus sharpness indicator of {sVal_img[k]} at each pixel, from the ray tracing.\n" \
  "\n" \
  "    All input images in all frames should have the same pixel column and row" \
  " counts {NX,NY}.  The {sVal} images may be in color, but are converted internally to" \
  " grayscale before any further processing.  The other images are greyscale.  Each" \
  " frame has a nominal focus plane height {zFoc[ki]} and nominal depth of" \
  " focus {zDep[ki]} used to generate the virtual view.\n" \
  "\n" \
  "    The images are read from files \"{inDir}/{frameDir}/{tag}.png\" where" \
  " {tag} is \"sVal\", \"hAvg\", etc; {frameDir} is \"frame-zf{FFF.FFFF}-df{DDD.DDDD}\";" \
  " and {FFF.FFFF} and {DDD.DDDD} are the focus plane" \
  " height {zFoc[ki]} and depth of focus {zDep[ki]}, both, formatted as \"%08.4f\"." \
  "\n" \
  "  The samples are assumed to be encoded" \
  " in the file with linear scale (gamma = 1).  The {hAvg} and {hDev} values are" \
  " assumed to have been scaled from {[0 _ zMax]} to {[0 _ 1]} before being" \
  " writen and are unscaled to {[0 _ zMax]} as" \
  " they are read.\n" \
  "\n" \
  "  TERMS TABLE\n" \
  "    The formula terms {term[0..NT-1]} are" \
  " computed from the basis coeffs as described" \
  " under {multifok_focus_op_terms_from_basis_coeffs}, using" \
  " the terms descriptions read from the specified {termSetFile}. See" \
  " {multifok_term_indices_from_names} and {multifok_term_read_term_table} for the file format.\n" \
  "\n" \
  "OUTPUTS\n" \
  "\n" \
  "  OUTPUT FOLDERS\n" \
  "    All output files associated with each frame {ki} are written out to a" \
  " folder \"{outDir}/{frameDir}/\" where, as above, {KKKKK} is the frame index {ki}, formatted" \
  " as \"%05d\".   Some files, not associated with any specific frame, are" \
  " written to folder \"{outDir}/\" itself.\n" \
  "\n" \
  "  BASIS IMAGE FILES\n" \
  "    For each input frame {ki} and each basis element {kb} in {0..NB-1}, the program writes" \
  " to \"{outDir}/{frameDir}/be{BBB}-bVal.png\" a grayscale" \
  " image showing the coefficient of that basis element computed on normalized window" \
  " samples; where {BBB} is the element index {kb} zero-padded to 3 digits." \
  "\n" \
  "  TERM IMAGE FILES\n" \
  "    For each input frame {ki} and each quadratic term {kt} in {0..NT-1}, the program writes" \
  " to \"{outDir}/{frameDir}/tm{TTT}-tVal.png\" a grayscale" \
  " image showing the value of that quadratic term computed from the basis" \
  " coefficients above; where {TTT} is the term index {kt} zero-padded to 3 digits." \
  "\n" \R
  "  LOCAL AVERAGE, GRADIENT, AND DEVIATION IMAGE FILES\n" \
  "    For each input frame {ki} the program computes images {sAvg}, {sGrd}, and {sDev} with" \
  " the weighted window sample average, gradient, gradient modululs and deviation around each pixel of the scene view image {sVal[ki]}.  The {sGrd} image has three channels, respectively the {X} derivative {sGrx} of {sVal}, the {Y} derivative {sGry}, and the modulus {sGrm = hypot(sGrx,sGry)/sqrt(2)}.   The {sDev} image is the weighted RMS value of the residual window sample values after subtracting the affine ramp {sAvg + sGrx*x + sGry*y}, where {x} and {y} are sample indices relative to the window center.  The program writes" \
  " these images to \"{outDir}/{frameDir}/{tag}.png\" where {tag} is \"sAvg\", \"sGrd\", and \"sDev\", respectively.\n" \
  "\n" \
  "  LOCALLY NORMALIZED IMAGE FILES\n" \
  "    For each input frame {ki}, the program writes" \
  " to \"{outDir}/{frameDir}/sNrm.png\" a grayscale image {sNrm} that is the normalized  version of {ifr.sVal}, showing the central window pixel after" \
  " removing the local linear ramp component {sAvg + sGrx*x + sGry*y} and dividing the residual by {hypot(sDev,sampleNoise)}.  The image {sNrm} is scaled from {[-1_+1]} to {[0_1]} on writing.\n" \
  "\n" \
  "  REGRESSION DATA FILE\n" \
  "    Writes to \"{outDir}/odata.txt\" a file with columns\n" \
  "\n" \
  "    ?? P{ki}.{ix}.{iy} {sAvg} {sGrd} {sDev} {shrp} {hDif} {hDev} {coeff[0]} .. {coeff[NB-1]} {term[0]} .. {term[NT-1]} \n" \
  "\n" \
  " where \n" \
  "\n" \
  "    {ki} is the image index.\n" \
  "    {ix} and {iy} are the column and row of the pixel.\n" \
  "    {sAvg} is the average of the {sVal} samples in the window centered at that pixel.\n" \
  "    {sGrx} is the gradient modulus of {sVal} within that window.\n" \
  "    {sGry} is the gradient modulus of {sVal} within that window.\n" \
  "    {sGrm} is the gradient modulus {hypot(sGrx,sGry)}.\n" \
  "    {sDev} is the rms of window samples values after removing the average and gradient.\n" \
  "    {shrp} is the \"actual\" sharpness of focus indicator at that pixel (as read from the {ifr.shrp} image).\n" \
  "    {hAvg} is the average scene {Z} in the pixel.\n" \
  "    {hDev} is the deviation of the {Z}-coord of the scene within that pixel.\n" \
  "    {hDif} is the difference between {hAvg} and {zFoc}.\n" \
  "    {bVal[0..NB-1]} are the basis coefficients computed from the normalized window samples.\n" \
  "    {tVal[0..NT-1]} are the quadratic terms computed from those coefficients.\n" \
  "\n" \
  "  The values of {sAvg} and {sDev} are computed from the window of {sVal} centered" \
  " at the pixel, before they are normalized, taking the window sample weights" \
  " into account.\n" \
  "\n" \
  "  The lines are sorted by pixel indices {ix,iy}, then by frame index {ki}.  A blank" \
  " line is written after the {NI} lines of each pixel.\n" \
  "\n" \
  "  PLOT MASK IMAGE FILE\n" \
  "    Writes to \"{outDir}/pSel.png\" a binary image file that is 1 for pixels" \
  " that were written out to the regression data file, 0 otherwise.\n" \
  "\n" \
  "  BASIS NAMES, PRODUCTS, AND TERM FILE\n" \
  "    Also writes, for documentation, the following files:\n" \
  "\n" \
  "      \"{outDir}/bnames.txt\"  names of the {NB} basis elements, like \"FXY\";\n" \
  "      \"{outDir}/prix.txt\"    table mapping products of pairs of coefficients to terms;\n" \
  "      \"{outDir}/tnames.txt\"  list of quadratic term formulas, like \"FX*FX+FY*FY\"." \
  "\n" \
  "    See {multifok_term_read_index_table} for the format and semantics of the \"prix\" file." \
  "\n" \
  "OPTIONS\n" \
  "  ???\n" \

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <bool.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_map_channels.h>

#include <multifok_window.h>
#include <multifok_image.h>
#include <multifok_scene.h>
#include <multifok_frame.h>
#include <multifok_stack.h>
#include <multifok_basis.h>
#include <multifok_term.h>
#include <multifok_test.h>

typedef struct mfdo_options_t 
  { /* Input images: */
    char *inDir;                       /* Directory where input image files live. */
    char *sceneType;                   /* Scene type shown in the frame images. */
    char *pattern;                     /* Name of pattern used to texturize the scene. */
    int32_t imageSize_NX;              /* Image pixels per row.  */
    int32_t imageSize_NY;              /* Image pixels per column.  */
    int32_t NI;                        /* Number of frames in stack, */
    double *zFoc;                      /* Focus plane position of each frame. */
    double *zDep;                      /* Simulated camera depth of focus of each frame. */
    multifok_basis_type_t basisType;   /* Local operator basis type. */
    multifok_window_type_t windowType; /* Window weights distribution type. */
    char *termSetFile;                 /* Name of file with formulas of quadratic terms. */
    double sampleNoise;                /* Noise level to assume when normalizing window samples. */
    char *outDir;                      /* Prefix for output filenames. */
  } mfdo_options_t;
  /* Command line parameters. */
   
typedef struct mfdo_result_frame_t
  { 
    int32_t NX;              /* Image pixel colums. */
    int32_t NY;              /* Image pixel rows. */
    int32_t NB;              /* Count of basis elements. */
    int32_t NT;              /* Count of quadratic term elements. */
    int32_t HW;              /* Invalid margin around images. */
    
    double zFoc;             /* Focus plane {Z} coordinate. */
    double zDep;             /* Depth of focus. */
    
    float_image_t **bVal;      /* Basis coeffs, indexed {0..NB-1}. */
    float_image_t **tVal;   /* Quadratic terms, indexed {0..NT-1}. */
    
    float_image_t *sAvg;       /* Window average image. */
    float_image_t *sGrd;       /* Window gradient image ({X} and {Y} derivatives, and modulus). */
    float_image_t *sDev;       /* Weighted  RMS of window residual after removing averge and gradient. */
    
    float_image_t *sNrm;    /* Locally normalized image. */
  } mfdo_result_frame_t;
  /* A record {ofr} of type {mfdo_result_frame_t} 
    is a set of images derrived from an input frame {ifr}.
    
    All images in {ofr} have the same pixel size, {NX} columns and {NY}
    rows, and a single channel -- except where noted otherwise.

    Each pixel of each of the images {ofr.bVal[0..NB-1]} contain the
    coefficient of the corresponding element of a
    specified basis of {NB} linear window operators, fitted to the
    window of the {ifr.sVal} input image around that pixel.

    The images {ofr.tVal[0..NT-1]} contain {NT} quadratic operators
    derived from those coefficients.

    Each pixel in the images {ofr.sAvg}, {ofr.sGrd}, and {ofr.sDev} are
    the average {sAvg}, the window gradient {(sGrx,sGry,sGrm)}, and the
    RMS of the non-affine residual of {ifr.sVal}, respectively, computed
    over the window around that pixel. See {PROG_INFO} for details.
    
    Each pixel of {ofr.sNrm} is the locally normalized image {ifr.sVal}. See
    {PROG_INFO} for details. */
 
mfdo_result_frame_t *mfdo_result_frame_new
  ( int32_t NX,
    int32_t NY,
    int32_t NB,
    int32_t NT
  );
  /* Creates a {mfdo_result_frame_t} record and its images, with the given parameters.
    The images are allocated with size {NX} by {NY}, assuming {NB} elements
    in the linear window operator basis and {NT} quadratic terms.
    The {HW} field is set to {-1}. The fields {zFoc} and {zDep} are set to {NAN}*/

typedef struct mfdo_result_stack_t
  { int32_t NI;                  /* Number of frames. */
    
    int32_t NX;                  /* Image pixel colums. */
    int32_t NY;                  /* Image pixel rows. */
    int32_t HW;                  /* Invalid margin around images (window radius). */
    
    int32_t NB;                  /* Count of basis elements. */
    int32_t NT;                  /* Count of quadratic term elements. */
    
    mfdo_result_frame_t** frame; /* Result frames, indexed {0..NI-1}. */
    
    float_image_t *pSel;         /* Pixel selection mask. */
  } mfdo_result_stack_t;
  /* A stack of {NI} result frames. 
  
    The {pSel} image is a binary mask that tells which pixels have been
    written out to the correlation data file (1) or omitted from it (0). */
  
mfdo_result_stack_t *mfdo_result_stack_new(int32_t NI);
  /* Allocates the result stack {ostack} with space for {NI} frames. the
    table {ostack.frame} is allocated with {NI} elemets, but they are
    all set to {NULL}. The fields {.NX,.NY,,HW,.NB,.NT} of {ostack} are
    set to {-1}. */

int32_t main(int32_t argc, char **argv);
 
multifok_term_set_t *mfdo_get_term_set
  ( multifok_basis_type_t bType,
    char *termSetFile,
    int32_t NW,
    double ws[]
  );
  /* Creates a linear window operator basis with the window weights {ws[0..NW*NW-1]} 
    and reads from file {termSetFile} a set of quadratic operator terms.
    The file must NOT have the weights column. */
   
multifok_stack_t *mfdo_read_stack(mfdo_options_t *o);
  /* Reads the {NI=o.NI} frames, each with the four images (blurred scene view, {Z} average,
    {Z} deviation, actual sharpness).  Returns them as a 
    {multifok_stack_t} record.
    
    All images must have {NX=o.imageSize_NX} columns by {NY=o.imageSize_NY} rows of pixels. 
    Frame {ki} is assumed to have simulated focus
    plane height {o.zFoc[ki]} and depth of focus {o.zDep[ki]},
    for {ki} in {0..NI-1}.
    
    The frames are read from files
    "{o.inDir}/frame-zf{FFF.FFFF}-df{DDD.DDDD}/{tag}.png", where
    {FFF.FFFF} and {DDD.DDDD} are {o.zFoc[ki]} and {o.zDep[ki]}
    formatted as "%08.4f", and {tag} is "sVal" (simulated view), "hAvg"
    ({Z}-coordinates), "hDev" ({Z}-deviation), and "shrp" ("actual"
    sharpness). */
 
mfdo_result_stack_t *mfdo_process_stack
  ( multifok_stack_t *istack,
    int32_t NW,
    double ws[],
    mfdo_options_t *o
  );
  /* Processes a multifocus stack of synthetic limited-focus-depth
    frames, producing a stack of corresponding derived frames
    through {mfdo_process_frame}.
    
    ??Does NOT fill the pixel selection image {ostack.pSel}. */
  
mfdo_result_frame_t *mfdo_process_frame
  ( int32_t ki, 
    multifok_frame_t *ifr,
    int32_t NW,
    double ws[],
    double sampleNoise,
    multifok_term_set_t *tset
  );
  /* Processes the frame {ifr} whose index is {ki}. 
    Computes the derived images (basis coefficients, quadratic terms,
    average, gradient, and residue, etc.). */

void mfdo_write_correl_data
  ( char *outDir,
    multifok_stack_t *istack,
    mfdo_result_stack_t *ostack
  );
  /* Writes to {outDir}/odata.txt} a subset of the pixel data that can
    be used to determine the best focus estimator by regression. Writes
    one line for each pixel considered useful for plotting and analysis,
    as explained in {PROG_INFO}.
    
    Also, computes a mask image {ostack.pSel} that shows which pixels
    had their data written out. */
  
void mfdo_write_result_stack(mfdo_result_stack_t *ostack);
  /* Writes the resulting stack images. 
  
    For each frame {ofr} with index {ki} writes 
    a set of images whose names begin with {framePrefix = "{outDir}/{frameDir}/{tag}.png"}
    where {frameDir} and {tag} are as explained in {PROG_INFO}. */

void mfdo_fit_quadratic(int32_t NI, double u[], double v[], double ev[], double *u_max_P, double *v_max_P, double *A_P);
  /* Estimates the maximum of a quadratic function {v = A*u^2 + B*u + C} with {A < 0},
    given at least three points {u[i],v[i]} of its graph and the estimated uncertainty {ev[i]}
    of {v[i]}, for {i} in {0..NI}.  Returns in {*u_max_P} the estimated 
    value {u_max = -B/(2*A)} of {u} where {Q(u)} is maximum, in {*v_max_P} the 
    estimated maximum {Q(u_max) = C + B^2/(4*A)}, and in {*A_P} the coefficient {A}.  */

bool_t mfdo_pixel_is_useful
  ( double A,
    double B,
    double C,
    double hAvg_min, 
    double hAvg_max, 
    double hDev_max, 
    int32_t *NXY_R_sDev_bz_P,
    int32_t *NXY_R_hDev_hi_P,
    int32_t *NXY_R_hAvg_ex_P,
  );
  /* Decides if a pixel is useful for analysis or regression, based on statistics
    of that pixel over all frames {ki}.
    
    Considers the coefficients {A,B,C} of the quadratic {Q(Z)} fitted to {(zFoc[ki],sVar[ki])}, 
    the max and min {hAvg_min,hAvg_max} over all {ki} of of the average {istack.frame[ki].hAvg} of simulated scene {Z} in pixel,
    and the max {hDev_max} of the deviation {istack.frame[ki].hDev} of the simulated scene {Z} in pixel. 
    
    Increments the counters on these conditions:
      
      {*NXY_R_sDev_bz_P} if the quadratic {Q} is bizarre.
      {*NXY_R_hDev_hi_P} if {hDev_max} is too high.
      {*NXY_R_hAvg_ex_P} if {hAvg_min} or {hAvg_max} are too close to the scene's min and max {Z}.

    The quadratic {Q} is considered bizarre if {A} is positive or {|A|}
    is too small or the {Z} value that maximizes {Q} is not well within the
    scene's {Z} range.
  */

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    mfdo_options_t *o = mfdo_parse_options(argc, argv);
    
    int32_t NW = 3;
    double *ws = multifok_window_weights(NW, o->windowType);

    /* Read the images: */
    multifok_stack_t *istack = mfdo_read_stack(o);

    mfdo_result_stack_t *ostack = mfdo_process_stack(istack, NW, ws, o);
    mfdo_write_correl_data(ostak);
    mfdo_write_result_stack(ostack);
    
    return 0;
  }
    
multifok_stack_t *mfdo_read_stack(mfdo_options_t *o)
  {
    int32_t NI = o->NI;
    int32_t NX = o->imageSize_NX;
    int32_t NY = o->imageSize_NY;
    char *stackDir = NULL;
    asprintf(&stackDir, "%s/st%s-%04dx%04d-%s", o->inDir, o->sceneType, NX, NY, o->pattern);
    
    bool_t gray = TRUE;
    multifok_stack_t *istack = multifok_stack_read(stackDir, gray, NI, o->zFoc, o->zDep, o->zMax);
    
    assert(istack->NI == NI);
    demand(istack->NX == NX, "inconsistent frame {NX}");
    demand(istack->NY == NY, "inconsistent frame {NY}");
     
    return istack;    
  }

multifok_term_set_t *mfdo_get_term_set
  ( multifok_basis_type_t bType,
    char *termSetFile,
    int32_t NW,
    double ws[]
  )
  {
    fprintf(stderr, "generating the linear operator basis...\n");
    bool_t ortho = TRUE;
    multifok_basis_t *basis = multifok_basis_make(bType, NW, ws, ortho);
    multifok_basis_print(stderr, basis);
    multifok_basis_ortho_check(stderr, basis);
    
    fprintf(stderr, "reading the quadratic term set...\n");
    FILE *rd = open_read(termSetFile, TRUE);
    bool_t weights = FALSE;
    bool_t verbose = TRUE;
    multifok_term_set_t *tset = multifok_term_set_read(rd, weights, basis, NULL, verbose);
    return tset;
  }
  
mfdo_result_frame_t *mfdo_result_frame_new
  ( int32_t NX,
    int32_t NY,
    int32_t HW,
    int32_t NB,
    int32_t NT
  )
  {
    mfdo_result_frame_t *ofr = talloc(1, mfdo_result_frame_t);
    
    ofr->NX = NX;
    ofr->NY = NY;
    ofr->HW = HW;

    ofr->NB = NB;
    ofr->NT = NT;
    
    ofr->zFoc = NAN;
    ofr->zDep = NAN;
    
    fprintf(stderr, "allocating images for basis coefs squared...\n");
    ofr->bVal = talloc(NB, float_image_t*);
    for (int32_t kb = 0; kb < NB; kb++) 
      { ofr->bVal[kb] = float_image_new(1, NX, NY); }

    fprintf(stderr, "allocating images for quadratic term values...\n");
    ofr->tVal = talloc(NT, float_image_t*);
    for (int32_t kt = 0; kt < NT; kt++) 
      { ofr->tVal[kt] = float_image_new(1, NX, NY); }

    fprintf(stderr, "allocating window average, gradient, and deviation images...\n");
    ofr->sAvg = float_image_new(1, NX, NY);
    ofr->sGrd = float_image_new(3, NX, NY);
    ofr->sDev = float_image_new(1, NX, NY);
      
    fprintf(stderr, "allocating the locally normalized image...\n");
    ofr->sNrm = float_image_new(1, NX, NY);
    
    return ofr;
  }
  
mfdo_result_stack_t *mfdo_result_stack_new(int32_t NI)
  {
    mfdo_result_stack_t *ostack = talloc(1, mfdo_result_stack_t);
    
    ostack->NI = NI;

    ostack->NX = -1;
    ostack->NY = -1;
    ostack->HW = -1;
    ostack->NB = -1;
    ostack->NT = -1;
    
    ostack->frame = talloc(NI, mfdo_result_frame_t*);
    for (int32_t ki = 0; ki < NI; ki++) { ostack->frame[ki] = NULL; }

    fprintf(stderr, "allocating the pixel selection mask...\n");
    ostack->pSel = float_image_new(1, NX, NY);
      
    return ostack;
  }

mfdo_result_stack_t *mfdo_process_stack
  ( multifok_stack_t *istack,
    int32_t NW,
    double ws[],
    mfdo_options_t *o
  )
  {
    demand((NW%2) == 1, "window size must be odd");
    
    int32_t NI = istack->NI;
    int32_t NX = istack->NX;
    int32_t NY = istack->NY;

    multifok_term_set_t *tset = mfdo_get_term_set(o->basisType, o->termSetFile, NW, ws);
    multifok_basis_t *basis = tset->basis;
    int32_t NB = basis->NB;
    int32_t NT = tset->NT;
     
    /* File with basis element names: */
    multifok_basis_elem_names_write_named(o->outDir, basis);
    multifok_term_set_write_named(o->outDir, tset, FALSE, NULL);
    multifok_term_set_product_table_write_named(o->outDir, tset->NP, tset->prix, TRUE);
   
    mfdo_result_stack_t *ostack = mfdo_result_stack_new(NI);
    
    ostack->NX = NX;
    ostack->NY = NY;
    ostack->HW = (NW-1)/2;
    
    ostack->NB = NB;
    ostack->NT = NT;
 
    for (int32_t ki = 0; ki < NI; ki++)
      { multifok_frame_t *ifr = istack->frame[ki];
        assert(ifr->NX == NX);
        assert(ifr->NY == NY);
        
        mfdo_result_frame_t *ofr = mfdo_process_frame(ki, ifr, NW, ws, o->sampleNoise, tset);
        
        assert(ofr->NX == NX);
        assert(ofr->NY == NY);
        assert(ofr->NB == NB);
        assert(ofr->NT == NT);
        assert(ofr->HW == HW);

        assert(ostack->frame[ki] == NULL);
        ostack->frame[ki] = ofr;
      } 
      
    return ostack;
  }

mfdo_result_frame_t *mfdo_process_frame
  ( int32_t ki, 
    multifok_frame_t *ifr,
    int32_t NW,
    double ws[],
    double sampleNoise,
    multifok_term_set_t *tset
  )
  {
    int32_t NX = ifr->NX;
    int32_t NY = ifr->NX;
    
    multifok_basis_t *basis = tset->basis;
    
    int32_t NB = basis->NB;
    int32_t NT = tset->NT;
    
    int32_t NS = NW*NW;

    mfdo_result_frame_t *ofr = mfdo_result_frame_new(NX, NY, HW, NB, NT);

    ofr->zFoc = ifr->zFoc;
    ofr->zDep = ifr->zDep;
    
    float_image_check_size(ifr->sVal, 1, NX, NY);

    /* Enumerate and process all windows in the frame: */
    int32_t HW = NW/2;
    float fsmp[NS];
    double dsmp[NS];
    double bVal[NB];       /* Coefficients of normalized window in be basis. */
    double tVal[NT];        /* Quadratic terms. */
    int32_t NL_scan = 0;    /* Number of pixels processed. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Get the samples in the window and remove mean and deviation: */
            float_image_get_window_samples(ifr->sVal, 0,ix,iy, NW, NW, FALSE, fsmp);
            for (int32_t ks = 0; ks < NS; ks++) { dsmp[ks] = fsmp[ks]; }
            
            /* Compute basis coefficients of raw window samples and save them in the {bVal} images: */
            multifok_basis_compute_coeffs(dsmp, basis, bVal);
            for (int32_t kb = 0; kb < NB; kb++)
              { float_image_set_sample(ofr->bVal[kb], 0, ix, iy, (float)bVal[kb]); }
              
            /* Compute the quadratic terms and save in the {tVal} images: */
            multifok_term_values_from_basis_coeffs(bVal, tset, tVal);
            for (int32_t kt = 0; kt < NT; kt++)
              { float_image_set_sample(ofr->tVal[kt], 0, ix, iy, (float)(tVal[kt])); }

            /* Compute the average and deviation of scene view image: */
            double sAvg, sGrx, sGry, sDev; /* Weighted window average and deviation */
            multifok_window_compute_average_gradient_and_deviation(NW, dsmp, ws, &sAvg, &sGrx, &sGry, &sDev); 
            double sGrm = hypot(sGrx,sGry);
            float_image_set_sample(ofr->sAvg, 0, ix, iy, (float)sAvg);
            float_image_set_sample(ofr->sGrd, 0, ix, iy, (float)sGrx);
            float_image_set_sample(ofr->sGrd, 1, ix, iy, (float)sGry);
            float_image_set_sample(ofr->sGrd, 2, ix, iy, (float)(sGrm/M_SQRT2));
            float_image_set_sample(ofr->sDev, 0, ix, iy, (float)sDev);
            
            /* Compute the nrmalized image: */
            multifok_window_remove_average_and_gradient(NW, dsmp, ws, sAvg, sGrx, sGry);
            double sNrm = dsmp[HW*NW + HW]/sMag;
            double sMag = hypot(sDev, sampleNoise);
            float_image_set_sample(ofr->sNrm, 0, ix, iy, (float)sNrm);
            NL_scan++;
          }
      }
    return ofr;
  }
            
void mfdo_write_correl_data
  ( char *outDir,
    multifok_stack_t *istack,
    mfdo_result_stack_t *ostack
  )
  {
    int32_t NI = istack->NI; assert(ostack->NI == NI);
    int32_t NX = istack->NX; assert(ostack->NX == NX);
    int32_t NY = istack->NX; assert(ostack->NX == NX);

    int32_t NB = ostack->NB;
    int32_t NT = ostack->NT;
    int32_t HW = ostack->HW;

    /* File with pixel-wise data for regression: */
    FILE *wr_dat = multifok_test_open_text_file(outDir, "odata");
    
    double zFoc[NI];  /* Nominal focus plane {Z} for each frame. */
    double zDep[NI];  /* Nominal depth of focus for each frame. */
    for (int32_t ki = 0; ki < NI; ki++) 
      { multifok_frame_t *ifr = istack->frame[ki];
        zFoc[ki] = ifr->zFoc; zDep[ki] = ifr->zDep;
      }
    
    /* Enumerate and process all windows in the frame: */
    int32_t NXY_scan = 0;       /* Number of pixels examined. */
    int32_t NXY_R_sDev_bz = 0;  /* Number of pixels rejected because {sDev} is bizarre. */ 
    int32_t NXY_R_sDev_lo = 0;  /* Number of pixels rejected because {sDev} is too low. */ 
    int32_t NXY_R_hDev_hi = 0;  /* Number of pixels rejected for high {hDev}. */ 
    int32_t NXY_R_hAvg_ex = 0;  /* Number of pixels rejected because {hAvg} is too extreme. */ 
    int32_t NXY_dat = 0;        /* Number of pixels selected to be written to {wr_dat}. */
    for (int32_t ix = HW; ix < NX-HW; ix++)
      { for (int32_t iy = HW; iy < NY-HW; iy++) 
          { /* Extremal values of pixel badness criteria among all frames: */
            double hDev_max = -INF; /* Max simulated {Z} deviation at this pixel */
            double hAvg_min = +INF; /* Min simulated {Z} average at this pixel. */
            double hAvg_max = -INF; /* Max simulated {Z} average at this pixel. */
            double sVar[NI];  /* Value of {sDev^2} for all frames at this pixel. */
            for (int32_t ki = 0; ki < NI; ki++)
              { multifok_frame_t *ifr = istack->frame[ki];
                mfdo_result_frame_t *ofr = ostack->frame[ki];
                double sDev = float_image_get_sample(ofr->sDev, 0,ix,iy);
                sVar[ki] = sDev*sDev;
                double hAvg = float_image_get_sample(ifr->hAvg, 0,ix,iy);
                if (hAvg < hAvg_min) { hAvg_min = hAvg; }
                if (hAvg > hAvg_max) { hAvg_max = hAvg; }
                double hDev = float_image_get_sample(ifr->hDev, 0,ix,iy);
                if (hDev > hDev_max) { hDev_max = hDev; }
              }
            /* Fit a quadratic function {Q(z) = A*z^2 + B*z + C} to points {(sVar[i],zFoc[i])}: */
            double A, B, C; 
            mfdo_fit_quadratic(NI, zFoc, sVar, zDep, &A, &B, &C);

            /* Decide whether the pixel is useful for analysis: */
            bool_t pix_ok = mfdo_pixel_is_useful
              ( A, B, C, hAvg_min, hAvg_max, hDev_max, 
                &NXY_R_sDev_bz, &NXY_R_hDev_hi, &NXY_R_hAvg_ex 
              );

            float pSel;  /* Value to write in pixel selection mask image. */
            if (pix_ok)
              { pSel = 0.0; }
            else
              { /* Write the data for analysis: */
                for (int32_t ki = 0; ki < NI; ki++)
                  { multifok_frame_t *ifr = istack->frame[ki];
                    mfdo_result_frame_t *ofr = ostack->frame[ki];

                    fprintf(wr_dat, "P%d.%d.%d ", ki, ix, iy);

                    double sAvg = float_image_get_sample(ofr->sAvg, 0,ix,iy);
                    fprintf(wr_dat, " %12.6f ", sAvg);   /* Window sample average. */
                    double sGrx = float_image_get_sample(ofr->sGrd, 0,ix,iy);
                    double sGry = float_image_get_sample(ofr->sGrd, 1,ix,iy);
                    double sGrm = float_image_get_sample(ofr->sGrd, 2,ix,iy)*M_SQRT2;
                    fprintf(wr_dat, " %+12.6f %+12.6f %12.6f ", sGrx, sGry, sGrm);   /* Window sample gradient. */
                    double sDev = float_image_get_sample(ofr->sDev, 0,ix,iy);
                    fprintf(wr_dat, " %12.6f ", sDev);   /* Window sample deviation after removing average and gradient. */
                    double hAvg = float_image_get_sample(ifr->hAvg, 0,ix,iy);
                    fprintf(wr_dat, " %12.6f ", hAvg);  /* Avgrage scene {Z} in pixel. */
                    double hDev = float_image_get_sample(ifr->hDev, 0,ix,iy);
                    fprintf(wr_dat, " %12.6f  ", hDev);  /* Deviation of scene {Z} in pixel. */
                    double shrp = float_image_get_sample(ifr->shrp, 0,ix,iy);
                    fprintf(wr_dat, " %14.10f ", shrp); /* "Actual" sharpness indicator. */
                    fprintf(wr_dat, "   ");
                    for (int32_t kb = 0; kb < NB; kb++) 
                      { double bVal = float_image_get_sample(ofr->bVal[kb], 0,ix,iy);
                        fprintf(wr_dat, " %+16.12f", bVal);
                      }
                    fprintf(wr_dat, "   ");
                    for (int32_t kt = 0; kt < NT; kt++) 
                      { double tVal = float_image_get_sample(ofr->tVal[kt], 0,ix,iy);
                        fprintf(wr_dat, " %16.12f", tVal);
                      }
                    fprintf(wr_dat, "\n");

                NXY_dat++;
                pSel = 1.0;
              }
            float_image_set_sample(ostack->pSel, 0, ix, iy, pSel);
            NXY_scan++;
          }
      }
    fprintf(stderr, "scanned %d pixels\n", NXY_scan); 
    fprintf(stderr, "rejected %d pixels for high {sDev}\n", NXY_R_sDev_bz);  
    fprintf(stderr, "rejected %d pixels for high {hDev}\n", NXY_R_hDev_hi);  
    fprintf(stderr, "rejected %d pixels for low {hAvg}\n", NXY_R_hAvg_ex);  
    fprintf(stderr, "wrote %d pixels for analysis\n", NXY_dat); 

    demand(NXY_dat > 0, "ABORTED");
    fclose(wr_dat);
    
  }
  
void mfdo_write_result_stack(char *outDir, mfdo_result_stack_t *ostack)
  {
    fprintf(stderr, "writing the output stack to %s...\n", outDir);
    for (int32_t ki = 0; ki < NI; ki++)
      { mfdo_result_frame_t *ofr = ostack->frame[ki];
        char *frameDir = NULL;
        asprintf(&frameDir, "%s/frame-zf%08.4f-fd-%08.4f", outDir, ofr->zFoc, ofr->zDep);
        mfdo_write_result_frame(frameSir, ofr);
        free(frameDir);
      }

    /* Write the pixel selection mask image: */
    multifok_test_write_pixel_mask_image(outDir, pSel, NULL, stackPrefix);
  }

void mfdo_write_result_frame(int32_t ki, char *frameDir, mfdo_result_frame_t *ofr)
  {
    fprintf(stderr, "  writing output frame %d to %s...\n", ki, frameDir);

    /* Write the average, deviation, and normalized images: */
    multifok_test_write_window_average_image(ofr->sAvg, frameDir);
    multifok_test_write_window_gradient_image(ofr->sGrd, frameDir);
    multifok_test_write_window_deviation_image(ofr->sDev, frameDir);
    multifok_test_write_normalized_image(ofr->sNrm, frameDir);

    /* Write out the squared basis coefficient images: */
    for (int32_t kb = 0; kb < NB; kb++)
      { multifok_test_write_basis_coeff_squared_image(bVal[kb], kb, frameDir); }

    /* Write out the quadratic term images: */
    for (int32_t kt = 0; kt < NT; kt++)
      { multifok_test_write_quadratic_term_image(tVal[kt], kt, frameDir); }

  }

bool_t mfdo_pixel_is_useful
  ( double A,
    double B,
    double C,
    double hAvg_min, 
    double hAvg_max, 
    double hDev_max, 
    int32_t *NXY_R_sDev_bz_P,
    int32_t *NXY_R_sDev_lo_P,
    int32_t *NXY_R_hDev_hi_P,
    int32_t *NXY_R_hAvg_ex_P,
  )
  { 
    /* Reject pixels whose {hAvg} is too close to the scene extremal {Z}  values: */
    double hAvg_min_ok = 1.1;         /* Ignore pixels with actual scene {Z} below this level. */
    double hAvg_max_ok = zmax - 1.1;  /* Ignore pixels with actual scene {Z} above this level. */
    if ((hAvg_min < hAvg_min_ok) || (hAvg_max > hAvg_max_ok)) { (*NXY_R_hAvg_ex_P)++; return FALSE; }

    /* Reject pixels which have too hight {hDev}: */
    double hDev_max = 0.75; /* Ignore pixels which more than this variance. */
    if (hDev > hDev_max) { (*NXY_R_hDev_hi_P)++; return FALSE; }

    /* Reject pixels whose {sDev} varies in a bizarre way: */
    double A_min_ok = 0.001; /* Reject pixels with {|A|} below this value. */
    if ((A >= 0) || (fabs(A) < A_min_ok)) { *(NXY_R_sDev_bz_P) ++; return FALSE; }
    /* Estimate the max {sDev} from the quadratic: */
    sVar_max = C + B*B/(4*A); 
    if (sVar_max < 0) { *(NXY_R_sDev_bz_P) ++; return FALSE; }
    sVar_argmax = -B/(2*A);
    /* Check whether the position of the quadratic's maximum is within the valid {Z} range:
    if ((sVar_argmax < hAvg_min_ok) || (sVar_argmax > hAvg_max_ok)) { *(NXY_R_sDev_bz_P) ++; return FALSE; }
    
    /* Reject pixels which have low {sDev}: */
    sDev_max = sqrt(sVar_max);
    double sDev_max_ok = 0.02; /* Reject pixels with {sDev_max} below this value. */
    if (sDev_max < sDev_max_ok) { (*NXY_R_sDev_lo_P)++; return FALSE; } 

    return TRUE;
  }

mfdo_options_t *mfdo_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    mfdo_options_t *o = notnull(malloc(sizeof(mfdo_options_t)), "no mem");

    argparser_get_keyword(pp, "-inDir");
    o->inDir = argparser_get_next(pp);
    
    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_NX = (int32_t)argparser_get_next_int(pp,50,4096);
    o->imageSize_NY = (int32_t)argparser_get_next_int(pp,50,4096);
    
    double_vec_t zFoc = double_vec_new(10);
    double_vec_t zDep = double_vec_new(10);
    int32_t NI = 0; /* Number of input images. */
    while (argparser_keyword_present(pp, "-frame"))
      { double_vec_expand(&(zFoc), NI);
        imo.zFoc = argparser_get_next_double(pp, 0.0,2*o->zMax);
        double_vec_expand(&(zDep), NI);
        imo.zDep = argparser_get_next_double(pp, 0.01,99.99);
        NI++;
      }
    if (NI == 0) { argparser_error(pp, "must specify at least one \"-frame\""); }
    double_vec_trim(&(zFoc), NI); o->zFoc = zFoc.e
    double_vec_trim(&(zDep), NI); o->zDep = zDep.e

    argparser_get_keyword(pp, "-basisType");
    char *bType = argparser_get_next_non_keyword(pp);
    o->basisType = multifok_basis_type_from_text(bType, FALSE);
    if (o->basisType < 0) { argparser_error(pp, "invalid \"-basisType\""); }
      
    argparser_get_keyword(pp, "-windowType");
    char *wType = argparser_get_next_non_keyword(pp);
    o->windowType = multifok_window_type_from_text(wType, FALSE);
    if (o->windowType < 0) { argparser_error(pp, "invalid \"-windowType\""); }
      
    argparser_get_keyword(pp, "-termSetFile");
    o->termSetFile = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-sampleNoise");
    o->sampleNoise = argparser_get_next_double(pp, 0.0, 100.0);  

    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
