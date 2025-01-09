#ifndef pst_lamp_H
#define pst_lamp_H

/* pst_lamp.h -- modeling lamps (light sources). */
/* Last edited on 2025-01-03 12:32:40 by stolfi */

#include <bool.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <argparser.h>

#include <pst_basic.h>

typedef struct pst_lamp_t
  { r3_t dir;       /* Unit vector pointing towards the lamp's center. */
    frgb_t pwr;     /* {pwr[c]} is the lamps's intensity in channel {c}. */
    double crad;    /* Co-sine of the lamp's angular radius. */
  } pst_lamp_t;
  /* A canonical model for a /distant round lamp/. See
    {pst_lamp_model_INFO} below. The angular radius parameter {angrad}
    is actually specified as {crad=cos(angrad)} for greater
    efficiency. */

#define pst_lamp_model_INFO \
  "A lamp is specified by three parameters:" \
  " its /direction/ {DIR}, its /power/ {PWR}, and its /angular" \
  " radius/ {ANGRAD}. \n" \
  "\n" \
  "  The direction {DIR} is a unit-length vector pointing from the" \
  " scene towards the lamp's center. The power {PWR} a triplet of reals giving the intensity" \
  " of the lamp's light in each color channel; it is defined as" \
  " the apparent brightness in channel {c} of a white Lambertian surface" \
  " perpendicular to {DIR}, fully illuminated by the lamp and viewed" \
  " from a point anywhere in the illuminated" \
  " side." \
  "\n" \
  "  The lamp proper is modeled by a circular cap on the sky" \
  " sphere (a sphere of infinite radius surrounding the" \
  " scene), centered on the point at infinity with" \
  " direction {DIR}, with angular radius {ANGRAD}" \
  " (in radians, between 0 and {PI}).  For small values" \
  " of {ANGRAD}, it is approximately (to first first order) equivalent" \
  " to a disk of radius {R*ANGRAD} at some sufficiently large" \
  " distance {R} from the scene --- large enough for the" \
  " resulting illumination to be considered uniform all over the scene." \
  "\n" \
  "  As special cases, if {ANGRAD == 0}, the lamp is" \
  " a point light source; if {ANGRAD == PI/2}, the illumination is" \
  " is equivalent to that provided by a uniformly bright infinite" \
  " light wall perpedicular to {DIR}; if {ANGRAD == PI}, the" \
  " illumination is that of a spherical lamp surrounding the" \
  " scene, which provides an isotropic (``ambient'') light field. (In" \
  " the latter case, the {DIR} parameter is irrelevant).  For" \
  " {ANGRAD} between {PI/2} and {PI}, the effect is equivalent" \
  " to an uniform ambient light of intensity {PWR} minus the" \
  " illumination provided by a lamp at {-DIR} with angular" \
  " radius {PI/2-ANGRAD}, with smaller intensity.\n" \
  "\n" \
  "  In any case, the the parameter {ANGRAD} is relevant" \
  " only for surfaces that are nearly parallel to" \
  " {DIR}.  Specifically, when the angle {THETA} between" \
  " {DIR} and the surface's tangent plane decreases from {+ANGRAD} to" \
  " {-ANGRAD}, the apparent brightness of the surface falls" \
  " smoothly from 1 to a minimum value.  (More precisely, it" \
  " varies as a quadratic function of {cos(THETA)}, chosen" \
  " so as to interpolate between" \
  " the lit and unlit parts with first-order continuity).\n" \
  "\n" \
  "  Thus, the {ANGRAD} parameter affects only the" \
  " softness of the terminator (the boundary between the lit" \
  " and unlit parts of the object).  A positive {ANGRAD}" \
  " leads to a softer penumbra in the parts of the image" \
  " where the lamp is only partly visible.  The {ANGRAD}" \
  " parameter has no effect on the apparent color of surfaces" \
  " that are fully lighted or completely inside the proper shadow."

vec_typedef(pst_lamp_vec_t,pst_lamp_vec,pst_lamp_t *);
  /* A list of pointers to lamp descriptors. */
  
pst_lamp_t *pst_lamp_new(r3_t *dir, frgb_t *pwr, double crad);
  /* Creates a new lamp with the specified attributes. 
  
    If {dir} not {NULL} and {*dir} is a non-zero vector, the lamp's
    direction is set to that vector,normalizd to unit length; otherwise
    direction is set o {(NAN,NAN,NAN)},
    
    If {pwr} is anot {NULL} and {*pwr} is a color with finite
    components, the lamp's power is set to that color; otherwise it is
    set to {(NAN,NAN,NAN)},
    
    The {crad} parameter is the co-sine of the angular radius {ANGRAD}
    of the lamp, and must be either {NAN} or between {-1} and {+1}. */

/* SHADING */

double pst_lamp_geom_factor(r3_t *nrm, r3_t *dir, double crad);
  /* Computes the geometric factor (relative amount of light falling
    per unit area) for a surface with normal {nrm} and a single
    distant round lamp with angular radius {rad = acos(crad)}
    located in direction {dir}.
    
    The result is a number between 0 and 1. It is 1 iff {nrm == dir},
    minimal if {nrm == -dir}, and transitions gradually between the two in
    a range of directions approximately orthogonal to {dir}.
 
    When a perfectly white and Lambertian surface with normal vector
    {nrm} is illuminated by a monochromatic distant round lamp with
    direction {dir}, intensity vector {pwr}, and co-sine-radius {crad},
    the apparent lightness of that surface in channel {c} is simply
    {pwr[c]*pst_lamp_geom_factor(nrm,dir,crad)}. */

/* COMMAND LINE PARSING */

void pst_lamp_spec_write(FILE *wr, pst_lamp_t *src);
  /* Writes to {wr} the parameters of lamp {src}, in a format
    compatible with {pst_lamp_spec_parse}. */

pst_lamp_t *pst_lamp_spec_parse(argparser_t *pp, uint32_t *N_P);
  /* Parses a lamp specification from the command line, starting 
    with "lamp", "wall", or "ambient".  If there is none, 
    returns {NULL}.
  
    The syntax is described by {pst_lamp_spec_HELP} and
    {pst_lamp_spec_HELP_INFO}.  All parameters of the same lamp
    must appear together in the command line, starting at the 
    current position of the parser {pp}.
    
    The radius parameter, if given, is converted from degrees to the
    co-sine {src.crad}. The vector specified with "direction" is
    normalized to unit length, and the azimuth and elevation specified
    with "angles" are converted to a unit direction vector.
    
    If {N_P} is not {NULL}, allows also "array {N}" instead of "lamp",
    and returns the number {N} in {*N_P}. 
    
    If unspecified, the direction {src.dir} is set to {(0,0,0)} for
    "ambient" and {(NAN,NAN,NAN)} otherwise; the power {src.pwr} is set
    to {(NAN,NAN,NAN)}; and the radius co-sine {src.crad} is set to
    {NAN} for "lamp" and "array", 0 for "wall", and {-1} for "ambient".
    The caller must check for {NAN} values and provide appropriate
    defaults. */
  
/* Help/info for full lamp specification: */

#define pst_lamp_spec_params_list_INFO \
  "\"angles\", \"direction\", \"radius\", \"ambient\", \"wall\", and \"power\""
 
#define pst_lamp_spec_HELP_INFO(more_HELP_INFO,more_direction_INFO,more_radius_INFO,more_power_INFO) \
  "A light source description may be one of the following:\n" \
  "\n" \
  "      lamp {SOURCE_PARAMETERS}\n" \
  "        Specifies a very distant compact source in a specific" \
  " direction with the given parameters (see below).\n" \
  "\n" \
  more_HELP_INFO \
  "      wall {SOURCE_PARAMETERS}\n" \
  "        Specifies a light source that is an infinite luminous plane with the given parameters.  The plane will be perpendicular to the spefified direction.\n" \
  "\n" \
  "      ambient {SOURCE_PARAMETERS}\n" \
  "        Specifies an ambient (isotropic) light source.\n" \
  "\n" \
  "    The {SOURCE_PARAMETERS} are:\n" \
  "\n" \
  "      direction {DX} {DY} DZ}\n" \
  "        Specifies the direction of the lamp, in Cartesian" \
  " coordinates.  Only the direction of the vector is important," \
  " not its length.  This parameter excludes" \
  " the \"angles\" parameter, and is not allowed for \"ambient\" sources." more_direction_INFO "\n" \
  "\n" \
  "      angles {AZIMUTH} {ELEVATION}\n" \
  "        Specifies the spherical coordinates of the lamp" \
  " relative to the scene.  The {AZIMUTH} is measured" \
  " from the {X} axis towards the {Y} axis. The {ELEVATION} is" \
  " measured from the {XY} plane towards the Z axis.  Both" \
  " angles should be given in degrees.  This parameter excludes" \
  " the \"direction\" parameter, and is not allowed for \"ambient\" sources." more_direction_INFO "\n" \
  "\n" \
  "      radius {ANGRAD}\n" \
  "        Specifies the angular radius of the lamp, as seen from the" \
  " scene.  The value should be in degrees; " \
  " it should be non-negative and between 0 and 180.  A \"lamp\" with zero radius" \
  " is a distant point-like light source.   A \"lamp\" with radius 90 covers a whole hemisphere of the sky sphere, and is equivalent to a \"wall\" source.  A \"lamp\" with \"radius\" 180 is equivalent to an \"ambient\" source." more_radius_INFO "  This parameter is not allowed for \"ambient\" and \"wall\" sources.\n" \
  "\n" \
  "      power " frgb_parse_color_HELP "\n" \
  "        Specifies the intensity and color of the light produced" \
  " by the lamp. It consists of " frgb_parse_color_INFO ".  This is the apparent color that results when the light source illuminates a white Lambertian (matte) surface whose normal direction is the light source's direction." more_power_INFO ""

#endif
