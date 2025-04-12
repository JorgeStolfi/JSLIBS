#define PROG_NAME "test_spots_image"
#define PROG_DESC "tests {sve_minn_single_step} by generating a image with random spots"
#define PROG_VERS "1.0"

/* Last edited on 2025-04-03 14:20:16 by stolfi */ 
/* Created on 2025-03-16 by J. Stolfi, UNICAMP */

#define test_spots_COPYRIGHT \
  "Copyright © 2024  by the State University of Campinas (UNICAMP)"
  
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -imageSize {NX} {NY} \\\n" \
  "    -spotSize {SPOT_SIZE_MIN} {SPOT_SIZE_VAR} ] \\\n" \
  "    -relSpotDist {REL_SPOT_DIST} \\\n" \
  "    [ -maxSpots {MAX_SPOTS} ] \\\n" \
  "    -numPasses {NUM_PASSES} \\\n" \
  "    " argparser_help_info_HELP ""
  
#define PROG_INFO \
  "The program creates an image of size {NX} by {NY} with many round" \
  " dots in random places (at most {MAX_SPOTS}), and writes it" \
  " to \"out/png-{NX}x{NY}/spots-{SPOT_SIZE_VAR}{SPOT_SIZE_MIN}.png\" with; where" \
  " {NX} and {NY} are formatted as \"%04d\" and {SPOT_SIZE_MIN} and {SPOT_SIZE_VAR} are" \
  " formatted as \"%01d\".  The actual spot radius {rad} is an exponential" \
  " function of {SPOT_SIZE_MIN} and {SPOT_SIZE_MIN+SPOT_SIZE_VAR}, where 0 is about one pixel.\n" \
  "\n" \
  "  The {REL_SPOT_DIST} parameter is" \
  " a multiplier for the sum of the radii of two spots that gives the ideal" \
  " distance between their centers.  The value for {REL_SPOT_DIST} affects" \
  " the number of spots that the program chose to create, which, in any" \
  " case, is limited to {MAX_SPOTS}.\n" \
  "\n" \
  "  The spot locations are chosen by minimizing a repulsion energy function" \
  " between spot pairs.   The energy is the sum of terms, one for each pair" \
  " of spots.  The energy of a pair of near spots is inversely" \
  " proportional to the distance between centers divided by the sum" \
  " of their radii.   For more distant pairs, there is an exponential decay" \
  " factor that effectively limits the influence of a spot to the few" \
  " nearest ones.  The image domain is implicitly assumed to have" \
  " toroidal topology+geometry, so that the energy of a pair actually" \
  " takes into account copies of the spot that are translated by" \
  " multiples of {NX} and {NY}.\n" \
  "\n" \
  "  The minimization algorithm is not particularly efficient.  The program" \
  " does {NUM_PASSES} iterations where each pass adjust the center of one spot" \
  " assuming all other spots are fixed.  It should be good enough since the initial" \
  " state has the spots fairly evenly distributed.\n" \
  "\n" \
  "  The program also writes a" \
  " file \"out/png-{NX}x{NY}/spots-{SPOT_SIZE_VAR}{SPOT_SIZE_MIN}-epair.txt\" with" \
  " a graph of the single-pair energy function as a function of" \
  " center-to-center distance assuming unit spot radius.  If {NX} and {NY} are zero, this" \
  " is the only output of the program.\n" \
  "\n" \
  "  The program also writes an" \
  " image \"out/png-{NX}x{NY}/spots-{SPOT_SIZE_MIN}-efunc.fni\" showing" \
  " the energy function for some arbitrarily chosen spot, determined" \
  " by all other spots, as a function for that spot's center."

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <argparser.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <bool.h>
#include <r2.h>
#include <rmxn.h>
#include <rmxn_throw.h>
#include <rmxn_spin.h>
#include <rmxn_shift.h>
#include <rmxn_regular_simplex.h>
#include <sve_minn.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_write_gen.h>
#include <sve_minn.h>

#include <float_image_paint.h>

typedef struct tspi_options_t
  { int32_t imageSize_NX;   /* Image cols count. */
    int32_t imageSize_NY;   /* Image rows count. */
    uint32_t spotSize_min;  /* Integer code for min spot radius. */
    uint32_t spotSize_var;  /* Integer code for variability of spot radii. */
    double relSpotDist;     /* Distance between nearest spots as fraction of sum of radii. */
    uint32_t maxSpots;      /* Max number of spots. */
    uint32_t numPasses;     /* Number of center optimization passes. */
  } tspi_options_t;
               
void tspi_compute_spot_radius_range(uint32_t spotSize_min, uint32_t spotSize_var, double *radMin_P, double *radMax_P);
  /* Returns in {*radMin_P} and {*radMax_P} the min and max radius (in
     pixels) of the spots for the user-specified
     {spotSize_min,spotSize_var}. */
    
uint32_t tspi_choose_number_of_spots
  ( int32_t NX, int32_t NY,
    double radMin, double radMax,
    double relSpotDist,
    uint32_t maxSpots
  );
  /* Chooses a proper number of spots for an image with {NX} cols and {NY} rows 
    assuming each spot has average area {avgArea} and that the typical distance
    between adjacent spots is {relSpotDist} times the sum of the radii.
    The number is limited to {maxSpots} or less. */

void tspi_throw_initial_centers_and_radii
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    double radMin, double radMax,
    r2_t sctr[], double srad[]
  );
  /* Throw {NS} initial spot centers {sctr[0..NS-1]} in the rectangle
    {[0 _ NX) x [0 _ NY)} and their radii {srad[0..NS-1]} in the range
    {[radMin _ radMax]}. The centers are mostly at jittered centers of
    cells of a grid. */

void tspi_relax_centers
  ( char *outPrefix,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    double stepMax,
    uint32_t numPasses
  );
  /* Adjusts the centers {sctr[0..NS-1]} so as to minimize the total
    energy as computed by {tspi_total_energy}. Assumes spot {ks} has
    center {sctr[ks]} and radius {srad[ks]}. The parameters
    {NX,NY,dRef,UX,UY} are passed to {tspi_total_energy} (q. v.). The
    {stepMax} is the max displacement allowed for each center in the
    minimization.
    
    Also writes images
    "{outPrefix}-{NNNNN}.png" showing the evolution of the spots; where
    {NNNNN} is the pass index from 0 to {numPasses-1}, formatted as
    '%05d'. */
  
void tspi_relax_center
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    uint32_t kc,
    double stepMax
  );
  /* Performs one step of the simplex-vertex-edge (SVE) minimizer on the
    center {sctr[kc]} trying to minimize the total energy of that spot,
    as given by {tspi_single_spot_energy}. The spot centers and radii
    are assumed to be {sctr[0..NS-1]} and {srad[0..NS-1]}. The
    parameters {NX,NY,dRef,UX,UY} are passed to {tspi_single_spot_energy}
    (q. v.).  The {stepMax} if the max displacement allowed for the center
    in the minimization. */
       
void tspi_reduce_center(int32_t NX, int32_t NY, r2_t *ck);
  /* Reduces the coordinates of {ck} to the image domain {[0 _ NX) × [0 _ NY)}
    by adding or subtracting {NX,NY}. */
 
void tspi_reduce_centers(int32_t NX, int32_t NY, uint32_t NS, r2_t sctr[]);
  /* Reduces the coordinates of the centers {sctr[0..NS-1]]}
    to the image domain {[0 _ NX) × [0 _ NY)} by adding or subtracting {NX,NY}. */

double tspi_total_energy
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  );
  /* Computes the total energy of the spots with centers {sctr[0..NS-1]}
    radii {srad[0..NS-1]}, and nominal relative distance {dRef}, by
    calling {tspi_single_spot_energy} for all spots and adding up the
    results. The parameters {NX,NY,dRef,UX,UY} are passed to
    {tspi_single_spot_energy} (q. v.). */
  
void tspi_single_spot_energy
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    uint32_t kc,
    uint32_t jc0, uint32_t jc1,
    double *ESpot_P,
    double *dRelMin_P
  );
  /* Returns in {*Espot_P} the sum of spot pair energies between the
    spot {kc} and every spot {jc} in {jc0..jc1}, except {kc} itself.
    Uses {tspi_single_pair_energy_toroidal}, therefore takes into
    account the toroidal geometry. The parameters
    {NX,NY,sctr,srad,dRef,UX,UY} are passed to
    {tspi_single_pair_energy_toroidal} (q. v.).
    
    If {dRelMin_P} is not null, also sets {*dRelmin_P} to the minimum of the
    distance between {sctr[kc]} and {sctr[jc]} for all {jc} above,
    including the replicated copies of {jc}, divided by the sum of the
    respective radii. */

void tspi_single_pair_energy_toroidal
  ( int32_t NX,  int32_t NY,
    r2_t *ck, double rk, 
    r2_t *cj, double rj,
    double dRef, 
    int32_t UX, int32_t UY,
    double *EToro_P,
    double *dRelMin_P
  );
  /* Returns in {*EToro_P} the energy of two spots {ks,js} with centers
    {ck,yk} and radii {rk,rj}, taking into account the toroidal
    geometry. Namely, calls {tspi_single_pair_energy} with centers {ck}
    and {cj+(ux*NX,uy*NY)} for all integer {ux} in {-UX..+UX} and {uy} in {-UY..+UY}.
    
    For each of those pairs of centers, computes the ratio
    {dRel=d/(rk+rj)}. Returns the minimum of those ratios in
    {*dRelMin_P}. */

int32_t tspi_find_toroidal_rep_max 
  ( int32_t NX, int32_t NY,
    double radMax, double dRef,
    uint32_t axis
  );
  /* Returns the largest non-negative {U} such that the single pair
    (non-toroidal) energy {E} between a spot at {(xk,yk)} and a
    {U}-replica of another spot {(xj,yj)} is significant. The
    {U}-replica is {(xj±U*NX,yj)} if {axis} is 0, or {(xj,yj±U*NY)} if
    {axis} is 1.  The spot radii are assumed to be both {radMax}. 
    
    The parameter {dRef} is passed to {tspi_single_pair_energy}. */
    
void tspi_single_pair_energy
  ( double dkj, double rk, double rj,
    double dRef,
    double *E_P,
    double *rTerm_P,
    double *eTerm_P
  );
  /* Returns in {*E_P} the energy of two spots with center-to-center
    distance {dkj} and radii {rk,rj}.
    
    The parameter {dRef} should be the typical distance (in pixels)
    between a spot and its efw nearest neighbors. The formula is such
    that the energy is proportional to {(rk+rj)/dkj} for small {dkj},
    but drops rapidly to zero as {dkj} becomes larger than a few times
    {dRef}.
    
    Does NOT take into account the toroidal geometry. However, the
    distance {dkj} may be between points outside the image domain rectangle.
    
    Also stores into {*rTerm_P} (if not {NULL}) the rational term of the
    energy. Also stores into {*eTerm_P} (if not {NULL}) the exponential
    decay term of the energy. */

void tspi_make_and_write_spots_image
  ( char *outPrefix,
    int32_t iter,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef
  );
  /* Creates a spots image with {tspi_make_spots_image} then writes 
    it out to "{outPrefix}-{NNNNN}.png" where {NNNNN} is {iter} formatted as '%05d'.
    However, if {iter} is 99999 or greater, it is printed as "{outPrefix}-final.png" instead.
    The parameter {dRef} is ultimately used by {tspi_single_pair_energy} (q. v.). */

float_image_t* tspi_make_spots_image
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef
  );
  /* Creates the spots image given the image size {NX,NY},
    the list of centers {sctr[0..NS-1]} and radii {srad[0..NS-1]}.
    The parameter {dRef} is ultimately used by {tspi_single_pair_energy} (q. v.). */
   
void tspi_write_image(char *outPrefix, char *tag, char *stage, float_image_t *img, int32_t indent);
  /* Writes {img} out to disk as file "{outPrefix}-{tag}-{stage}.png".  The "-{tag}"
    is omitted if {tag} is {NULL}.  Ditto for the "-{stage}".
    
    If {indent} is non-negative, writes a message with the file name to {stderr},
    indented by {indent} spaces.  If {indent} is negative, the message is suppressed. */

void tspi_write_2D_energy_plot
  ( char *outPrefix,
    char *stage,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  );
  /* Writes a gnuplot-compatible file "{outPrefix}-eplot-{stage}.txt" with 
    data of the total energy as a function of the spot centers,
    along some 2D affine subspace of the space {(\RR^2)^NS} of all centers.
    
    Specifically, chooses two lists {u0[0..NS-1]} and {u1[0..NS-1]} of
    unit 2D-vectors using {tspi_choose_2D_energy_plot_vectors}. Then
    perturns {sctr[kc]} by {d0*u0[kc] + d1*u1[kc]}, for various values
    of {d0} and {d1} in a regular grid centered at {(0,0)}, and computes
    the total energy {Etot} using {tspi_total_energy} for the perturned
    centers, with the given parameters {srad,dRef,UX,UX}. For each
    pair {d0,d1}, writes a line with "{d0} {d1} {Etot}", with blank
    lines between rows. */

void tspi_choose_2D_energy_plot_vectors(uint32_t NS, r2_t **u0_P, r2_t **u1_P);
  /* Chooses two lists of 2D vectors {u0[0..NS-1]} and {u1[0..NS-1]} such that 
    {u0[kc]} and {u1[kc]} are random mutually orthogonal unit vectors,
    for each {kc}. */

void tspi_make_and_write_single_spot_energy_image
  ( char *outPrefix,
    char *stage, 
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  );
  /* Creates an energy plot image with {tspi_make_single_spot_energy_image} then writes 
    it out to "{outPrefix}-espot-{stage}.png". The parameters {dRef,UX,UY} are passed to
    {tspi_single_spot_energy} (q. v.). */
   
float_image_t* tspi_make_single_spot_energy_image
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  );
  /* Creates an image that shows the variation of
    {tspi_single_spot_energy} for various positions of the center of a
    selected spot {kc} with radius {srad[kc]}, assuming that every spot
    {jc} other than {kc} is fixed at {sctr[jc]} and has radius
    {sctr[j]}. The parameters {dRef,UX,UY} are passed to
    {tspi_single_spot_energy} (q. v.). */

void tspi_write_single_pair_energy_graph(char *outPrefix, double radMax, double dRef);
  /* Writes to a file "{outPrefix}-epair.txt" with the energy {E(dkj)}
    contributed by a single pair of spots of radius {radMax} as a function of the
    Euclidean distance {dkj} between the centers.
    
    Uses {tspi_single_pair_energy} with reference distance {dRef} to compute {E(dkj)}. 
    Does not consider the toroidal topology of the image domain. */

tspi_options_t* tspi_parse_options(int32_t argc, char **argv); 
  /* Parses the command line options. */

int32_t main(int32_t argc, char **argv);

#define tspi_NS_MAX 10000
 /* Max number of spots in an image. */

int32_t main (int32_t argc, char **argv)
  {
    tspi_options_t *o = tspi_parse_options(argc, argv);
    
    int32_t NX = o->imageSize_NX, NY = o->imageSize_NY;

    mkdir("out", 0755);
    char *outDir = jsprintf("out/png-%04dx%04d", NX, NY);
    mkdir(outDir, 0755);
    char *outPrefix = jsprintf("%s/spots-%01d%01d", outDir, o->spotSize_var, o->spotSize_min);
    
    double radMin, radMax;
    tspi_compute_spot_radius_range(o->spotSize_min, o->spotSize_var, &radMin, &radMax);
    fprintf(stderr, "actual spot radius = [ %.2f _ %.2f ] pixels\n", radMin, radMax);

    uint32_t NS = tspi_choose_number_of_spots(NX, NY, radMin, radMax, o->relSpotDist, o->maxSpots);
    fprintf(stderr, "will generate image with {NS} = %d spots\n", NS);
    assert(NS >= 2);

    double avgSpotDist = sqrt(NX*NY/(double)NS/sqrt(3)); /* Typical dist between centers (px). */
    fprintf(stderr, "assumed average center spacing {avgSpotDist} = %12.8f\n", avgSpotDist);
    
    double dRef = avgSpotDist; /* Reference distance for energy cutoff. */
    fprintf(stderr, "energy cutoff distance {dRef} = %12.8f\n", dRef);

    fprintf(stderr, "writing single pair energy graph ...\n");
    tspi_write_single_pair_energy_graph(outPrefix, radMax, dRef);

    fprintf(stderr, "choosing %d initial centers ...\n", NS);
    r2_t *sctr = talloc(NS, r2_t);
    double *srad = talloc(NS, double);
    tspi_throw_initial_centers_and_radii(NX, NY, NS, radMin, radMax, sctr, srad);

    fprintf(stderr, "determining the replication limits {UX,UY} ...\n");
    int32_t UX = tspi_find_toroidal_rep_max(NX, NY, radMax, dRef, 0);
    int32_t UY = tspi_find_toroidal_rep_max(NX, NY, radMax, dRef, 1);
    fprintf(stderr, "  UX = %d  UY = %d ...\n", UX, UY);

    fprintf(stderr, "making initial single-spot energy image ...\n");
    tspi_make_and_write_single_spot_energy_image(outPrefix, "beg", NX, NY, NS, sctr, srad, dRef, UX, UY);

    fprintf(stderr, "writing initial 2D total energy plot ...\n");
    tspi_write_2D_energy_plot(outPrefix, "beg", NX, NY, NS, sctr, srad, dRef, UX, UY);

    double stepMax = avgSpotDist;
    fprintf(stderr, "relaxing the centers - max step = %12.8f ...\n", stepMax);
    tspi_relax_centers(outPrefix, NX, NY, NS, sctr, srad, dRef, UX, UY, stepMax, o->numPasses);

    fprintf(stderr, "making final energy image ...\n");
    tspi_make_and_write_single_spot_energy_image(outPrefix, "end", NX, NY, NS, sctr, srad, dRef, UX, UY);

    fprintf(stderr, "writing final 2D total energy plot ...\n");
    tspi_write_2D_energy_plot(outPrefix, "end", NX, NY, NS, sctr, srad, dRef, UX, UY);

    free(sctr);
    free(srad);

    free(outDir);
    free(outPrefix);
    free(o);
    
    fprintf(stderr, "done.\n");
    return 0;
  }

void tspi_compute_spot_radius_range(uint32_t spotSize_min, uint32_t spotSize_var, double *radMin_P, double *radMax_P)
  { 
    (*radMin_P) = pow(M_SQRT2, spotSize_min + 2);
    (*radMax_P) = pow(M_SQRT2, spotSize_min + spotSize_var + 2);
  }
 
uint32_t tspi_choose_number_of_spots
  ( int32_t NX, int32_t NY,
    double radMin, double radMax,
    double relSpotDist,
    uint32_t maxSpots
  )
  { uint32_t NS;
    if ((NX == 0) || (NY == 0))
      { NS = 0; }
    else
      { double areaImg = NX*NY; /* Area of image in pixels. */

        double areaMin = M_PI*radMin*radMin;     /* Area of smallest spot. */
        double areaMax = M_PI*radMax*radMax;     /* Area of largest spot. */

        double areaMed = (2*areaMax + areaMin)/3; /* Mean spot area. */

        double areaNom = areaMed*relSpotDist*relSpotDist; /* Mean area of spot plus spot spacing. */

        double gapFactor = sqrt(3)/(M_PI/2);   /* Accounts for gaps in hex packing. */

        double areaEff = areaNom*gapFactor;    /* Estimated avg area of spot with spacing and gaps. */
        NS = (uint32_t)floor(areaImg/areaEff); /* Number of spots. */
        
        /* Try to increment {NS} until it is rel prime to {NX} and {NY}: */
        fprintf(stderr, "  tentative spot count = %d\n", NS);
            
        while ((gcd(NS, (uint32_t)NX) != 1) || (gcd(NS, (uint32_t)NY) != 1)) { NS++; }
        fprintf(stderr, "  next relative prime spot count = %d\n", NS);

        if (NS >maxSpots) 
          { fprintf(stderr, "  clipping spot count from %d to %d\n", NS, maxSpots);
            NS = maxSpots;
          }
      }
    return NS;
  }

void tspi_throw_initial_centers_and_radii
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    double radMin, double radMax,
    r2_t sctr[], double srad[]
  )
  { 
    double aspect = ((double)NX)/((double)NY);
    int32_t KX = (int32_t)ceil(sqrt(NS*aspect)) + 2; /* Cand position cols. */
    int32_t KY = (int32_t)ceil(sqrt(NS/aspect)) + 2; /* cand position rows. */
    int32_t KP = KX*KY; /* Number of candidate positions. */
    
    uint32_t kp = 0; /* Candidate positions tried so far. */
    uint32_t ns = 0; /* Spot centers chosen so far. */
    for (int32_t ky = 0; ky < KY; ky++)
      { for (int32_t kx = 0; kx < KX; kx++)
          { /* Compute the probability {prob} of choosing the candidate: */
            double prob = ((double)NS - ns)/((double)KP - kp);
            assert(prob <= 1.0); 
            if (drandom() <= prob)
              { /* Compute jittered coords of candidate: */
                double y = ((double)NY*(ky + 0.5 + dabrandom(-0.3,+0.3)))/((double)KY);
                double x = ((double)NX*(kx + 0.5 + dabrandom(-0.3,+0.3)))/((double)KX);
                /* Store the candidate: */
                assert(ns < NS);
                sctr[ns] = (r2_t){{ x, y }};
                srad[ns] = dabrandom(radMin, radMax);
                ns++;
              }
            kp++;
          }
      }
    assert(ns == NS);
  }
    
void tspi_reduce_centers(int32_t NX, int32_t NY, uint32_t NS, r2_t sctr[])
  { 
    for (int32_t kc = 0; kc < NS; kc++)
      { tspi_reduce_center(NX, NY, &(sctr[kc])); }
  }

void tspi_reduce_center(int32_t NX, int32_t NY, r2_t *ck)
  { double *xk = &(ck->c[0]), *yk = &(ck->c[1]);
    /* Reduce to image domain modulo {NX,NY}: */
    /* Reduce to image domain modulo {NX,NY}: */
    while ((*xk) < 0) { (*xk) += NX; }
    while ((*xk) >= NX) { (*xk) -= NX; }
    while ((*yk) < 0) { (*yk) += NY; }
    while ((*yk) >= NY) { (*yk) -= NY; }
  }

void tspi_relax_centers
  ( char *outPrefix,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    double stepMax,
    uint32_t numPasses
  )
  { 
    bool_t debug = FALSE;
    
    /* Relax the spot centers: */
    fprintf(stderr, "performing %d passes ...\n", numPasses);

    for (int32_t kpass = 0; kpass < numPasses; kpass++)
      { /* One more pass over the whole dot collection: */
        fprintf(stderr, "  pass %d ...\n", kpass);
        tspi_make_and_write_spots_image(outPrefix, kpass, NX, NY, NS, sctr, srad, dRef);
        for (uint32_t kc = 0; kc < NS; kc++)
          { if (debug) { fprintf(stderr, "    relaxing spot center %d ...\n", kc); }
            tspi_relax_center(NX, NY, NS, sctr, srad, dRef, UX, UY, kc, stepMax);
            if (debug) { fprintf(stderr, "    ..............................\n"); }
          }

        fprintf(stderr, "\n");
      }

    fprintf(stderr, "creating final image ...\n");
    tspi_make_and_write_spots_image(outPrefix, INT32_MAX, NX, NY, NS, sctr, srad, dRef);
  }

void tspi_relax_center
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    uint32_t kc,
    double stepMax
  )
  { 
    bool_t debug = FALSE;
    uint32_t nz = 2;  /* Number of variables to optimize: */
    double z_ini[nz]; /* Main work vector for the optimization variables. */
               
    auto double goal(uint32_t nz_cur, const double z_cur[]);
      /* Stores the tentative values of variables {z_cur[0..nz_cur-1]} into
        {sctr[kc]} and evaluates the repulsion energy of the spot {kc}.
        Requires {nz_cur=nz=2}. */
    
    /* Copy coordinates of vertex {sctr[kc]} to {z[0..nz-1]}: */
    if (debug) { fprintf(stderr, "      extracting coordinates of spot center ...\n"); }
    z_ini[0] = sctr[kc].c[0];
    z_ini[1] = sctr[kc].c[1];
    
    /* Compute energy at this initial config: */
    if (debug) { fprintf(stderr, "      computing spot energy of current position ...\n"); }
    double Fz_ini = goal(nz, z_ini);
    if (debug) { fprintf(stderr, "      Fz_ini = %24.16e\n", Fz_ini); }
    
    if (debug) { fprintf(stderr, "      allocating work storage for simplex and probe values ...\n"); }
    uint32_t NV = nz+1; /* Number of simplex vertices. */
    double v[NV*nz];  /* Simplex vertices. */
    uint32_t NF = (nz+1)*(nz+2)/2;
    double Fv[NF];  /* Goal function at vertices and midpoints. */
    
    double z_fin[nz];
    rn_copy(nz, z_ini, z_fin);
    double Fz_fin = Fz_ini;
    
    double dMax = stepMax;               /* Max dist from starting position. */
    bool_t dBox = FALSE;                 /* Ball-shaped search domain. */
    double rSim = 1.0;  /* Initial simplex radius (pixels). */
    double rMin = 0.01; /* Min simplex radius (pixels). */
    uint32_t nEvals = 0;
    double dStep = NAN;
    bool_t debug_probes = FALSE;
    if (kc == NS/2) { fprintf(stderr, "        dMax = %12.8f  rSim = %12.8f  rMin = %12.8f\n", dMax, rSim, rMin); }
    uint32_t stepKind = sve_minn_single_step
      ( nz, z_fin, &Fz_fin, goal, -1,
        z_ini, dMax, dBox, &rSim, rMin, 
        v, Fv, &nEvals, &dStep, debug, debug_probes
      );
    if (debug) { fprintf(stderr, "        moved to %s point\n", (char*[3]){"stationary", "simplex", "previous"}[stepKind]); }
    
    if (Fz_fin > Fz_ini)
      { fprintf(stderr, "        failed, reverting ...\n");
        Fz_fin = goal(nz, z_ini);
      }
    else
      { if (debug) { fprintf(stderr, "        succeeded\n"); } }

    if (debug) { fprintf(stderr, "      Fz_fin = %24.16e\n", Fz_fin); }

    if (debug) { fprintf(stderr, "      reducing centers to domain ...\n"); }
    tspi_reduce_centers(NX, NY, NS, sctr);

    return;
    
    double goal(uint32_t nz_cur, const double z_cur[])
      { assert(nz == 2); 
        /* Insert the variables into the {sctr} vector: */
        if (debug) { fprintf(stderr, "        inserting vars into center vector ...\n"); }
        sctr[kc].c[0] = z_cur[0];
        sctr[kc].c[1] = z_cur[1];
        tspi_reduce_center(NX, NY, &(sctr[kc]));
        /* Compute the energy: */
        if (debug) { fprintf(stderr, "        computing the energy ...\n"); }
        double ESpot;
        tspi_single_spot_energy(NX, NY, NS, sctr, srad, dRef, UX, UY, kc, 0, NS-1, &ESpot, NULL);
        return ESpot;
      }
  }
 
float_image_t* tspi_make_spots_image
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef
  )
  {
    float_image_t *img = float_image_new(1, NX,NY);
    float_image_fill_channel(img, 0, 0.0);
    
    /* Reduce centers to domain, just in case: */
    tspi_reduce_centers(NX, NY, NS, sctr);
    
    for (int32_t kc = 0; kc < NS; kc++)
      { r2_t *ck = &(sctr[kc]);
        double rk = srad[kc];
        bool_t round = TRUE;
        bool_t diagonal = FALSE; /* Irrelevant actually. */
        float vfill = 1.0;
        float vdraw = NAN;
        uint32_t msamp = 3;
        for (int32_t ux = -2; ux <= +2; ux++)
          { for (int32_t uy = -2; uy <= +2; uy++)
              { float_image_paint_dot
                  ( img, 0, ck->c[0] + ux*NX, ck->c[1] + uy*NY, rk, 0.1*rk, 
                    round, diagonal, vfill, vdraw, msamp
                  );
              }
          }
      }
    return img;
  }

double tspi_total_energy
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  )
  {
    double ETot = 0;
    for (uint32_t kc = 1; kc < NS; kc++)
      { /* Consider only previous spots: */
        uint32_t jc0 = 0, jc1 = kc - 1; assert(jc1 < NS);
        double ESpot;
        tspi_single_spot_energy(NX, NY, NS, sctr, srad, dRef, UX, UY, kc, jc0, jc1, &ESpot, NULL);
        ETot += ESpot;
      }
    return ETot;
  }

void tspi_single_spot_energy
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY,
    uint32_t kc,
    uint32_t jc0, uint32_t jc1,
    double *ESpot_P,
    double *dRelMin_P
  )
  { double ESpot = 0;
    double dRelMin = +INF;
    for (uint32_t jc = jc0; jc <= jc1; jc++)
      { if (jc != kc)
          { double Ekj, dRelMinkj;
            tspi_single_pair_energy_toroidal
              ( NX, NY, 
                &(sctr[kc]), srad[kc], &(sctr[jc]), srad[jc], 
                dRef, UX, UY,
                &Ekj, &dRelMinkj
              );
            ESpot += Ekj;
            dRelMin = fmin(dRelMin, dRelMinkj);
         }
      }
    if (dRelMin_P != NULL) { (*dRelMin_P) = dRelMin; }
    (*ESpot_P) = ESpot;
  }
    
void tspi_single_pair_energy_toroidal
  ( int32_t NX,  int32_t NY,
    r2_t *ck, double rk, 
    r2_t *cj, double rj,
    double dRef,
    int32_t UX, int32_t UY,
    double *EToro_P,
    double *dRelMin_P
  )
  { demand((UX >= 1) && (UY >= 1), "invalid {UX,UY}");
    double xk = ck->c[0], yk = ck->c[1];
    assert((xk >= 0) && (xk < NX) && (yk >= 0) && (yk < NY));
    double xj = cj->c[0], yj = cj->c[1];
    assert((xj >= 0) && (xj < NX) && (yj >= 0) && (yj < NY));
    
    /* Initialize {EToro} with the pair energy with no toroidal replication: */
    double EToro = 0;
    double dRelMin = +INF;
    
    /* Add to {EToro} the energy from off-axis replications: */
    for (int32_t ux = -UX; ux <= +UX; ux++)
      { for (int32_t uy = -UY; uy <= +UY; uy++)
          { double xju = xj + ux*NX, yju = yj + uy*NY;
            double E;
            double dkju = hypot(xk - xju, yk - yju);
            double dRelkju = dkju/(rk + rj);
            tspi_single_pair_energy(dkju, rk, rj, dRef, &E, NULL, NULL);
            EToro += E;
            if (dRelkju < dRelMin) { dRelMin = dRelkju; }
          }
      } 
    (*EToro_P) = EToro;
    if (dRelMin_P != NULL) { (*dRelMin_P) = dRelMin; }
  }

int32_t tspi_find_toroidal_rep_max 
  ( int32_t NX, int32_t NY,
    double radMax, double dRef,
    uint32_t axis
  ) 
  {
    /* Explore replication in {X} or {Y} until insignificant: */
    double ETot = 0.0;
    int32_t U = 0; 
    bool_t signif;
    double xk = 0.5*NX, yk = 0.5*NY;
    do 
      { U++;
        signif = FALSE;
        for (int32_t sgn = -1; sgn <= +1; sgn += 2)
          { double xj = (axis == 0 ? sgn*U*NX : 0);
            double yj = (axis == 1 ? sgn*U*NY : 0);
            double dkj = hypot(xk - xj, yk - yj);
            double E;
            tspi_single_pair_energy(dkj, radMax, radMax, dRef, &E, NULL, NULL);
            if (fabs(E) > 1.0e-9*fabs(ETot)) { signif = TRUE; }
            ETot += E;
          }
      }
    while (signif);
    return U;
  }

void tspi_single_pair_energy
  ( double dkj, double rk, double rj,
    double dRef,
    double *E_P,
    double *rTerm_P,
    double *eTerm_P
  )
  { 
    double dRel = dkj/(rk + rj);
    
    /* Rational term: */
    double rTerm = 1/dRel;
    
    /* Exponential term: */
    double uRel = dkj/dRef;
    /* We want {zr} to be about 0 for {uRel <~ 1.0}, and about 8 for {uRel > 6}: */
    double C = 1.5;
    double zr = C*(hypot(1.0, uRel) - 1.0);
    double eTerm = exp(-zr*zr/2);

    double E = rTerm*eTerm;
    (*E_P) = E;
    if (rTerm_P != NULL) { (*rTerm_P) = rTerm; }
    if (eTerm_P != NULL) { (*eTerm_P) = eTerm; }
  }
  
void tspi_write_single_pair_energy_graph(char *outPrefix, double radMax, double dRef)
  { char *fname = jsprintf("%s-epair.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);

    /* Assumed parameters: */
    double rk = radMax, rj = radMax;    /* Assumed for this plot. */
    
    /* Compute an energy value {ERef} at typical distances: */
    double dkjRef = dRef;
    double ERef;
    tspi_single_pair_energy(dkjRef, rk, rj, dRef, &ERef, NULL, NULL);
    fprintf(stderr, "  reference distance {dkjRef} = %12.8f energy {ERef} = %+12.8f \n", dkjRef, ERef);
    
    /* Plot energy varying the distance until it is not significant: */
    double dMin = rk+rj;
    double dStep = 0.01*dMin;
    double E;
    double dkj = dMin;
    do
      { double rTerm, eTerm;
        tspi_single_pair_energy(dkj, rk, rj, dRef, &E, &rTerm, &eTerm);
        fprintf(wr, "%8.4f %+12.8f %+12.8f %+12.8f \n", dkj, E, rTerm, eTerm);
        dkj += dStep;
      }
    while (E > 1.0e-10*ERef);
    fclose(wr);
    free(fname);
  }
  
#define tspi_fname_iter_FMT "%05d"
#define tspi_fname_iter_MAX 99999
  
#define tspi_imageSize_MIN 16
#define tspi_imageSize_MAX 4096
    
void tspi_make_and_write_spots_image
  ( char *outPrefix,
    int32_t iter,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef
  )
  { float_image_t *img = tspi_make_spots_image(NX, NY, NS, sctr, srad, dRef);
    bool_t final = (iter > tspi_fname_iter_MAX);
    char *stage = (final ? jsprintf("final") : jsprintf(tspi_fname_iter_FMT, iter));
    tspi_write_image(outPrefix, NULL, stage, img, (final ? 0 : 2));
    float_image_free(img);
    free(stage);
  }

void tspi_make_and_write_single_spot_energy_image
  ( char *outPrefix,
    char *stage,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  )
  { float_image_t *erg = tspi_make_single_spot_energy_image(NX, NY, NS, sctr, srad, dRef, UX, UY);
    tspi_write_image(outPrefix, "espot", stage, erg, 0);
    float_image_free(erg);
  }
 
float_image_t* tspi_make_single_spot_energy_image
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  )
  {
    float_image_t *img = float_image_new(1, NX,NY);
    float_image_fill_channel(img, 0, 0.0);
    
    /* Reduce centers to domain, just in case: */
    tspi_reduce_centers(NX, NY, NS, sctr);
    
    /* Choose middling spot: */
    uint32_t kc = NS/2;
    
    /* Save its initial position: */
    r2_t ck_ini = sctr[kc];
    
    r2_t *ck = &(sctr[kc]);
    double rk = srad[kc];
    srad[kc] = 1.0e-6; /* Temporary. */
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { (*ck) = (r2_t){{ x + 0.5, y + 0.5 }};
            double ESpot, dRelMin;
            tspi_single_spot_energy(NX, NY, NS, sctr, srad, dRef, UX, UY, kc, 0, NS-1, &ESpot, &dRelMin);
            if (dRelMin >= 1.0) { float_image_set_sample(img, 0, x, y, (float)(ESpot)); }
          }
      }
    /* Restore {srad[kc]}: */
    srad[kc] = rk;
      
    /* Scale the samples to {0_1]: */
    float vMin = +INF, vMax = -INF;
    float_image_update_sample_range(img, 0, &vMin, &vMax);
    if (vMin < 0)
      { float vMag = fmaxf(fabsf(vMin), fabsf(vMax));
        vMin = -vMag;
        vMax = +vMag;
      }
    else
      { assert(vMax > 0.0);
        vMin= 0.0;
      }
      
    float_image_rescale_samples(img, 0, vMin, vMax, 0, 1);
      
    /* Restore spot {kc} to its original position: */
    sctr[kc] = ck_ini;
    return img;
  }

void tspi_write_2D_energy_plot
  ( char *outPrefix,
    char *stage,
    int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double dRef,
    int32_t UX, int32_t UY
  )
  {
    char *fname = jsprintf("%s-eplot-%s.txt", outPrefix, stage);
    FILE *wr = open_write(fname, TRUE);
    r2_t *u0, *u1;
    tspi_choose_2D_energy_plot_vectors(NS, &u0, &u1);
    double step = 0.125; /* Displacement of centers for each plot step (pixels). */
    int32_t np = 15;    /* Number of plot steps in each direction. */
    r2_t *pctr = talloc(NS, r2_t);
    for (int32_t k0 = -np; k0 <= +np; k0++)
      { for (int32_t k1 = -np; k1 <= +np; k1++)
          { for (int32_t kc = 0; kc < NS; kc++)
              { r2_t dctr;
                r2_mix(k0*step, &(u0[kc]), k1*step, &(u1[kc]), &dctr);
                r2_add(&dctr, &(sctr[kc]), &(pctr[kc]));
              }
            double ETot;
            if (hypot(k0, k1) <= np + 1.0e-6)
              { tspi_reduce_centers(NX, NY, NS, pctr);
                ETot = tspi_total_energy(NX, NY, NS, pctr, srad, dRef, UX, UY);
              }
            else
              { ETot = NAN; }
            fprintf(wr, "%+12.4f %+12.4f %20.10f\n", k0*step, k1*step, ETot);
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    free(fname);
    free(u0); free(u1); free(pctr);
  }
    
void tspi_choose_2D_energy_plot_vectors(uint32_t NS, r2_t **u0_P, r2_t **u1_P)
  { r2_t *u0 = talloc(NS, r2_t);
    r2_t *u1 = talloc(NS, r2_t);
    for (int32_t kc = 0; kc < NS; kc++)
      { /* Make {u0[kc],u1[kc]} two random orthogonal directions: */
        r2_throw_ortho_dirs(&(u0[kc]), &(u1[kc]));
      }
    (*u0_P) = u0;
    (*u1_P) = u1;
  }

void tspi_write_image(char *outPrefix, char *tag, char *stage, float_image_t *img, int32_t indent)
  {
    int32_t NC, NX, NY;
    float_image_get_size(img, &NC, &NX, &NY);
    
    if (tag == NULL) { tag = ""; }
    char *tagSep = (strlen(tag) == 0 ? "" : "-");
    if (stage == NULL) { stage = ""; }
    char *stageSep = (strlen(stage) == 0 ? "" : "-");
    char *fname = jsprintf("%s%s%s%s%s.png", outPrefix, tagSep, tag, stageSep, stage);
    if (indent >= 0) { fprintf(stderr, "%*swriting image %s\n", indent, "", fname); }

    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE;
    double expoEnc = 1.000;
    double bias = 0.000;
    bool_t verbose = FALSE;
    float vMin = 0.0;
    float vMax = 1.0;
    float_image_write_gen_named(fname, img, ffmt, yUp, vMin, vMax, expoEnc, bias, verbose);

    free(fname);
  }

tspi_options_t *tspi_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    tspi_options_t *o = talloc(1, tspi_options_t);

    /* PARSE POSITIONAL ARGUMENTS: */
    
    argparser_get_keyword(pp, "-imageSize");
    o->imageSize_NX = (int32_t)argparser_get_next_int(pp, tspi_imageSize_MIN, tspi_imageSize_MAX);
    o->imageSize_NY = (int32_t)argparser_get_next_int(pp, tspi_imageSize_MIN, tspi_imageSize_MAX);
    
    argparser_get_keyword(pp, "-spotSize");
    o->spotSize_min = (uint32_t)argparser_get_next_int(pp, 0, 9);
    o->spotSize_var = (uint32_t)argparser_get_next_int(pp, 0, 9 - o->spotSize_min);

    argparser_get_keyword(pp, "-relSpotDist");
    o->relSpotDist = argparser_get_next_double(pp, 1.1, 99.0);

    argparser_get_keyword(pp, "-numPasses");
    o->numPasses = (uint32_t)argparser_get_next_int(pp, 1, 1000);
    
    if (argparser_keyword_present(pp, "-maxSpots"))
      { o->maxSpots = (uint32_t)argparser_get_next_int(pp, 2, tspi_NS_MAX); }
    else
      { o->maxSpots = tspi_NS_MAX; }
    
    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
