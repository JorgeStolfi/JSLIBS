#ifndef pst_fit_light_H
#define pst_fit_light_H

/* pst_fit_light.h -- general fitting of light fields to images. */
/* Last edited on 2024-12-28 20:03:58 by stolfi */

#include <bool.h>
#include <r3.h>
#include <float_image.h>

#include <pst_shading.h>
#include <pst_light.h>
#include <pst_lamp.h>

/* GENERAL PARAMETERS */

/* 
  The following parameters have the same meaning for all procedures
  in this interface:
  
    {IMG}: a photo of the scene, the target of fitting.
    
    {c}: the channel of {IMG} to consider.
    
    {NRM}: the scene's normal map.
      
    {minNormalZ}: minimum Z coordinate for the surface normal.
      
    {lht}: the light model to assume for the scene.

*/
 
#define pst_fit_light_valid_pixels_INFO \
  "The light fitting procedures will ignore any" \
  " pixel {p} where the given surface normal {NRM[p]} is not finite or is the null" \
  " vector and/or the photo's intensity {IMG[p]} in any channel is exactly zero.  Thus, the user" \
  " may permanently exclude `bad' pixels from the light field fitting, by" \
  " setting their value to zero in the normal map and/or in the input photo."

#define pst_fit_light_steep_parts_INFO \
  "The procedures will also ignore any pixel where the normal {NRM[p]}" \
  " has a Z component less than {MINNORMALZ}.  This feature may be" \
  " helpful when the target of the fitting is placed over a white" \
  " background surface, which may contribute significant amounts" \
  " of scattered light to the steepest parts of the object's surface."

/* LIGHT FIELD MODEL */

#define pst_fit_light_model_INFO \
  pst_light_model_INFO "\n" \
  "  The light field is assumed be the same at every point on" \
  " the scene's surface, except that the flow is zero along any" \
  " direction that points out of the surface (which accounts" \
  " for the proper shadow effect)."
  
#define pst_fit_light_model_uniform_INFO \
  "  Except for proper shadows, the light field model does" \
  " not allow for any effect of the scene on the light"\
  " field.  In particular, it does not account for projected" \
  " shadows (that is, blockage of the light flow at one" \
  " point by parts of the scene located elsewhere), transparent" \
  " or translucent objects, radiosity, etc..\n"

/* SURFACE MODELS FOR LIGHT FITTING */

/* The procedures in this section compute some or all parameters of a
  light field model, so that the synthetic photo of a scene with given
  normal map {NRM}, under that light field, best matches a given photo
  {IMG} of the same scene. The light field model is described by
  {pst_fit_light_simple_model_INFO} and the scene's nature by
  {pst_fit_light_white_Lambertian_surface_INFO}.
  
  The fitting is done separately for each channel {c}. That produces a
  light model where all lamps have nonzero power only in channel {c}.
  After fitting all channels, fitted lamps with directions that are
  similar enough may be merged. If these lamps have separate
  singe-channel colors, the result will be a multichannel lamp.
  Alternatively, if lamps are known to be white, or of a specific color,
  the input image can be converted to grayscale first, fitted for
  channel 0 only, and the resulting monochomatic lamps converted to
  polychromatic ones. */

#define pst_fit_light_white_Lambertian_surface_INFO \
  "  The light field fitting procedures assume that"\
  " the visible surface of the scene is Lambertian"\
  " (purely diffusive, with no specular or glossy"\
  " reflection) and of uniform lightness (diffusion"\
  " coefficient).  The diffusion coefficient"\
  " is assumed to be 1; if that is not the case,"\
  " the fitted intensities for the bounded lamp"\
  " and ambient field will be multiplied by"\
  " the actual diffusion coefficient."

/* SINGLE-LIGHT FITTING */

void pst_fit_light_single_iterative
  ( float_image_t *IMG,  /* Photo of object. */ 
    float_image_t *NRM,  /* Normal map of object. */
    int32_t c,           /* Color channel to consider. */
    pst_light_t *lht,    /* (IN/OUT) Light model. */
    pst_lamp_t *src,     /* Lamp of {lht} to adjust. */
    double dirStep,      /* Max change allowed in {src} direction. */
    bool_t pwrAdjust,    /* TRUE optimizes the power of {src}, if possible. */
    bool_t ambAdjust,    /* TRUE optimizes the overall power of the other lamps. */
    double weightBias,   /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative,  /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,   /* Ignore pixels {p} where {NRM[p].z < minNormalZ}. */
    uint32_t iterations, /* Max iterations to use in fitting. */
    double tolerance     /* Iteration stopping criterion. */
  );
  /* Adjusts the direction and power of source {src}, as well as the
    overal power of the other lights in {lht}, so as to minimize the
    difference between the predicted image {F(lht,NRM)} and the given
    image {IMG}.

    If {dirStep} is positive, and {src} is not an ambient lamp, the
    procedure tries to optimize the direction of {src}, by at most that
    amount. Otherwise it treats that direction as fixed.
    
    If {pwrAdjust} is TRUE, the procedure tries to optimize the
    power of {src} in channel {c}. Otherwise it treats that
    power as fixed. 
    
    If {ambAdjust} is true, the procedure tries to mutiply the powers in
    channel {c} of the lamps of {lht} other than {src} by an optimal
    dimming factor. Otherwise it treats those light powers as fixed.
    
    If {weightBias} is finite, the least-squares functional uses weight
    {1/(IMG[p] + weightBias)} for each valid pixel {p}. If {weightBias} is
    infinite, all valid pixels have weight 1.
    
    If {nonNegative} is TRUE, lamp powers are constrained to be positive.
    
    In any case, the procedure may only adjust the power of lamps
    in channel {c}. The power of any lamp in other channels is not affected.
    
    When the procedure exits, the best-fitting parameters are stored
    into {lht}, including in {src}. 
    
    For more details, see {pst_fit_light_single_iterative_INFO},
    {pst_fit_light_simple_direction_by_least_squares_INFO} and
    {pst_fit_light_simple_intensities_by_least_squares_INFO} below.
    (In those doc strings, the parameter {iterations} is called
    {MAXITER}, {tolerance} is {TOL}, and {minNormalZ} is {MINZ}). */

#define pst_fit_light_single_iterative_INFO \
  "  The procedure adjusts the power of the selected lamp, and, if so" \
  " requested, also its direction. If there are two or more lamps, the" \
  " procedure may aslo be asked to adjust a single /dimming factor/ {dim}," \
  " that will be multiplied into the original power in channel {c} of all lamps other" \
  " than {src}.\n" \
  "\n" \
  "  The least squares criterion reduces to a linear system which is solved" \
  " directly.  However, when adjusting the direction of {src}, the fitting" \
  " can be applied only to those pixels where the lamp is fully" \
  " visible.  Those parts are found iteratively --- by guessing" \
  " the lamp's direction, excluding those parts of the image where the" \
  " normal points away from that direction, and recomputing the lamp's" \
  " direction from the remaining pixels.\n" \
  "\n" \
  "  This process usually requires several iterations to" \
  " converge (or may not converge at all).  The parameters {MAXITER} and" \
  " {TOLERANCE} provide the stopping criterion for this loop.  Namely, the" \
  " procedure stops after {MAXITER} iterations, or when the intensity of the" \
  " fitted light field between two consecutive iterations is {TOLERANCE}" \
  " or less."
  
#define pst_fit_light_single_iterative_dark_weight_INFO \
  "  When computing the least-squares error functional, valid" \
  " pixels in the image normally get the same weight (1), while" \
  " invalid pixels are excluded (0).  Optionally, the user may" \
  " request for each valid pixel {p} the weight {1/(IMG[p] + WTBIAS)}," \
  " where {IMG[p]} is the image value at {p}, and {WTBIAS} is a" \
  " specified positive constant.  With this option, the" \
  " darker parts of the image are considered more important" \
  " than the light parts."

/* NON-ITERATIVE, SINGLE-LIGHT FITTING BY LEAST SQUARES */

void pst_fit_light_single_lsq
  ( float_image_t *IMG, 
    float_image_t *NRM, 
    int32_t c,          /* Color channel to consider. */
    pst_light_t *lht,   /* Light model. */
    pst_lamp_t *src,    /* Lamp of {lht} to adjust. */
    double dirStep,     /* Max change allowed in {src} direction. */
    bool_t pwrAdjust,   /* TRUE tries to optimize the power of {src}, if possible. */
    bool_t ambAdjust,   /* TRUE tries to optimize the overall power of the other lamps. */
    double weightBias,  /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ   /* Ignore image points where the normal's Z is less than this. */
  );
  /* Similar to {pst_fit_light_single_iterative}, but instead
    iterative method uses a linear least-squares fitting method that 
    is appropriate for a Lambertian shading model {F(lht,NRM)}.
    
    When computing {src.dir}, the procedure considers only those parts
    of the scene which are fully illuminated by {src}, at its current
    position {src->dir}. The uncertainty {unc} is added to the lamp's
    angular radius for the purposes of deciding which pixels satisfy
    this condition.
    
    For more details, see {pst_fit_light_single_lsq_INFO} below. */

#define pst_fit_light_single_lsq_INFO \
  "  The procedure for adjusting a single light" \
  " uses the following fact: at any pixel {p} where surface" \
  " where surface is white Lambertian and the lamp is fully visible from" \
  " it, the apparent lightness predicted by the simple light" \
  " model is an affine (degree 1) function {F(NRM[p])} of the coordinates" \
  " of the normal direction vector {NRM[p]}.  The four" \
  " coefficients of the function {F} are the three coordinates" \
  " of the lamp's direction times its power, and the dimming factor" \
  " times the light of the remaining lamps.\n" \
  "\n" \
  "  Therefore, the procedure applies affine least squares fitting" \
  " of {F(NRM[p])} to {IMG[p]} (over the fully lighted parts" \
  " of the image); and then recovers the lamp's direction" \
  " and power, as well as the dimming factor, from the four" \
  " coefficients of the fitted function {F}."
  
#define pst_fit_light_single_lsq_caveats_INFO \
  "  The least-squares light fitting step may fail if the geometry" \
  " is not adequate --- in particular, if the lamp's angular radius is" \
  " greater than 60 degrees, or there are too few valid pixels" \
  " illuminated by the source in its current position, or their" \
  " normals are not sufficiently varied."

double pst_fit_light_lsq_pixel_weight
  ( double smp,
    r3_t *nrm, 
    double minNormalZ,
    r3_t *dir,
    double minCos,
    double weightBias
  );
  /* Computes the weight of a pixel with value {smp} in the 
    relevant channel and surface normal {nrm}, for least-squares fitting under 
    a single lamp model.  
    
    The weight is zero if {smp} is zero, if {nrm} is the null vector,
    or if the Z coordinate of {nrm} is less than {minNormalZ}.
    The weight is zero also if {dir} (the lamp's direction) is not NULL
    and the {dot(dir,nrm)} is less than {minCos}.
    
    Otherwise, the weight is 1 if {weightBias} is infinite.
    
    Otherwise it is {1/(smp+weightBias)}. */

/* NON-ITERATIVE, SINGLE LIGHT FITING BY TRIVIAL HEURISTIC */

void pst_fit_light_single_trivial
  ( float_image_t *IMG, 
    float_image_t *NRM, 
    int32_t c,          /* Color channel to consider. */
    pst_light_t *lht,   /* Light model. */
    pst_lamp_t *src,    /* Lamp of {lht} to adjust. */
    bool_t dirAdjust,   /* TRUE estimates the direction of {src}. */
    bool_t pwrAdjust,   /* TRUE estimates the power of {src}. */
    bool_t ambAdjust,   /* TRUE estimates a dimming factor for the other lamps. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ   /* Ignore image points where the normal's Z is less than this. */
  );
  /* Similar to {pst_fit_light_single_iterative}, but instead
    iterative method uses a rather crude heuristics. */

#define pst_fit_light_single_trivial_INFO \
  "  !!! WRONG DESCRIPTION, FIX IT !!! The parameters are estimated" \
  " by locating the brightest and" \
  " darkest pixels {pmax} and {pmin} in the image. The bounded lamp's" \
  " direction {bdir} is assumed to be the normal vector {NRM[pmax]}." \
  " The ambient lamp intensity {apwr} is assumed to be the color" \
  " {IMG[pmin]}, and the bounded lamp's intensity {bpwr} is assumed to" \
  " be {IMG[pmax] - IMG[pmin]}.\n" \
  "\n" \
  "  This heuristic assumes that some visible part of the surface is" \
  " perpendicular to the lamp's direction, while some other part is" \
  " turned completely away from the lamp. (These assumptions may not" \
  " be true for some scenes. In fact, they are never true if the lamp's" \
  " direction makes an obtuse angle with the camera's direction vector.)"

/* MULTI-LAMP LIGHT MODEL FITTING */

/* The procedures in this section compute some or all parameters of
  the the multi-lamp light field model so that the synthetic photo of a scene
  with given normal map {NRM}, under that light field, best matches a
  given photo {IMG} of the same scene. The light field model is
  described by {pst_fit_light_multi_model_INFO} and 
  the scene's nature by {pst_shading_Lambertian_INFO}. */

void pst_fit_light_multi
  ( float_image_t *IMG, 
    float_image_t *NRM,
    int32_t c, 
    pst_light_t *lht,
    double weightBias,  /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative, /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ
  );
  /* Updates the intensity parameter {pwr} of each lamp in {lht},
    without changing their directions and agular radii, so as to
    produce the best possible fit to {IMG} and {NRM}.
    
    If {weightBias} is finite, the least-squares functional uses weight
    {1/(IMG[p] + weightBias)} for each valid pixel {p}. If {weightBias} is
    infinite, all valid pixels have weight 1.
    
    The procedure considers only valid pixels (where {IMG} and {NRM}
    are non-zero). For best results, pixels in projected shadows and
    highlights should be excluded if possible. */

#define pst_fit_light_multi_INFO \
  "  The fitting procedure adjusts only the intensities of" \
  " the lamps inchannel {c}, without changing their directions or angular" \
  " radii.  It tries to minimize the sum of the squared" \
  " differences between {IMG[p]} and the expected apparent color" \
  " of a Lambertian white surface with normal {NRM[p]} under the" \
  " light field {lht}." \
  
/* COMMAND LINE ARGUMENT PARSING */

void pst_fit_light_parse_weightBias(argparser_t *pp, double *weightBias);
  /* Parses from the command line the {weightBias} parameter.
    The syntax is described in {pst_fit_light_parse_weightBias_HELP_INFO}.

    If {weightBias} is not NULL, the option "-weightBias {WTBIAS}" is
    allowed. In that case, the procedure sets {*weightBias=WTBIAS} if
    the option is present, and leaves it alone otherwise. */

#define pst_fit_light_parse_weightBias_HELP \
  "-weightBias {WTBIAS}"

#define pst_fit_light_parse_weightBias_INFO \
  "With this option, each valid pixel {p} in the image" \
  " gets weight {1/(IMG[p] + WTBIAS)}," \
  " in the least-squares quadratic functional, where" \
  " {IMG[p]} is the image value at {p}."

#define pst_fit_light_parse_weightBias_HELP_INFO \
  "  " pst_fit_light_parse_weightBias_HELP " \n" \
  "    " pst_fit_light_parse_weightBias_INFO


void pst_fit_light_parse_minNormalZ(argparser_t *pp, double *minNormalZ);
  /* Parses from the command line the {minNormalZ} parameter.
    The syntax is described in {pst_fit_light_parse_minNormalZ_HELP_INFO}.

    If {minNormalZ} is not NULL, the option "-minNormalZ {ZMIN}" is
    allowed. In that case, the procedure sets {*minNormalZ=ZMIN} if
    the option is present, and leaves it alone otherwise. */

#define pst_fit_light_parse_minNormalZ_HELP \
  "-minNormalZ {ZMIN}"

#define pst_fit_light_parse_minNormalZ_INFO \
  "If specified, this option excludes from the fitting" \
  " any parts of the photo where the Z coordinate of the" \
  " surface normal is less than {ZMIN}.  This option is" \
  " useful when analyzing the photo of an object against a bright" \
  " background, which may contribute significant secondary" \
  " lighting to near-vertical surfaces."

#define pst_fit_light_parse_minNormalZ_HELP_INFO \
  "  " pst_fit_light_parse_minNormalZ_HELP " \n" \
  "    " pst_fit_light_parse_minNormalZ_INFO

void pst_fit_light_parse_iterations(argparser_t *pp, uint32_t *iterations);
  /* Parses from the command line the {iterations} parameter.
    The syntax is described in {pst_fit_light_parse_iterations_HELP_INFO}.
    
    If {iterations} is not NULL, the option "-iterations {MAXITER}" is
    allowed. In that case, the procedure sets {*iterations} to {MAXITER} if
    the option is present, and leaves it alone otherwise. */

#define pst_fit_light_parse_iterations_HELP \
  "-iterations {MAXITER}"

#define pst_fit_light_parse_iterations_INFO \
  "Specifies the maximum number of iterations to use when" \
  " computing the bounded lamp's direction."

#define pst_fit_light_parse_iterations_HELP_INFO \
  "  " pst_fit_light_parse_iterations_HELP " \n" \
  "    " pst_fit_light_parse_iterations_INFO

void pst_fit_light_parse_tolerance(argparser_t *pp, double *tolerance);
  /* Parses from the command line the {tolerance} parameter.
    The syntax is described in {pst_fit_light_parse_tolerance_HELP_INFO}.
    
    If {tolerance} is not NULL, the option "-tolerance {TOL}" is
    allowed. In that case, the procedure sets {*tolerance} to {TOL} if
    the option is present, and leaves it alone otherwise. */

#define pst_fit_light_parse_tolerance_HELP \
  "-tolerance {TOL}"
  
#define pst_fit_light_parse_tolerance_INFO \
  "Stops the light field fitting algorithm when the lamp" \
  " directions computed in two cosecutive iterations differ" \
  " by at most {tol}."

#define pst_fit_light_parse_tolerance_HELP_INFO \
  "  " pst_fit_light_parse_tolerance_HELP " \n" \
  "    " pst_fit_light_parse_tolerance_INFO

#endif
