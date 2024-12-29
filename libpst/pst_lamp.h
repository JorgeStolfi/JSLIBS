#ifndef pst_lamp_H
#define pst_lamp_H

/* pst_lamp.h -- modeling lamps (light sources). */
/* Last edited on 2024-12-28 21:03:40 by stolfi */

#include <bool.h>
#include <r3.h>
#include <frgb.h>
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

pst_lamp_t *pst_lamp_spec_parse(argparser_t *pp, bool_t next);
  /* Parses a lamp specification from the command line.
  
    The syntax is described by {pst_lamp_spec_HELP} and
    {pst_lamp_spec_HELP_INFO}. All parameters of the same lamp
    must appear together in the command line.
    
    If {next} is true, checks only the next command line argument for
    the \"-lamp\" keyword, else looks for it anywhere in the command
    line. If the keyword \"-lamp\" is not found, returns NULL.
    
    See {argparser.h} for an explanation of the {pp} parameter. */
  
/* Help/info for full lamp specification: */

#define pst_lamp_spec_params_list_INFO \
  "\"angles\", \"direction\", \"radius\", \"ambient\", \"wall\", and \"power\""

#define pst_lamp_spec_HELP \
  "-lamp \\\n" \
  pst_lamp_spec_params_HELP
 
#define pst_lamp_spec_INFO \
  "Begins the description of a lamp (light source).  May be" \
  " followed by one or more lamp" \
  " parameters, " pst_lamp_spec_params_list_INFO "."
 
#define pst_lamp_spec_HELP_INFO \
  "  -lamp {LAMP_PARAMETERS}\n" \
  "    " pst_lamp_spec_INFO "\n" \
  "\n" \
  pst_lamp_spec_params_HELP_INFO

/* Parsing the lamp parameters. */

pst_lamp_t *pst_lamp_spec_params_next_parse(argparser_t *pp);
  /* Like {}, but assumes that the "-lamp" keyword has just been parsed.
    Parses only the lamp parameters, starting at the next command line
    argument. If unspecified, the direction is set to the null vector,
    the power is set to an empty vector, and {crad} is set to {+INF}. */

#define pst_lamp_spec_params_HELP \
  "      [ " pst_lamp_spec_ambient_HELP " \\\n" \
  "      | [ " pst_lamp_spec_angles_HELP " | " pst_lamp_spec_direction_HELP " ] \\\n" \
  "        [ " pst_lamp_spec_radius_HELP " | " pst_lamp_spec_wall_HELP "] \\\n" \
  "      ] \\\n" \
  "      [ " pst_lamp_spec_power_HELP " ]"

#define pst_lamp_spec_params_HELP_INFO \
  pst_lamp_spec_angles_HELP_INFO "  Excludes" \
  " the \"direction\" and \"ambient\" parameter.\n" \
  "\n" \
  pst_lamp_spec_direction_HELP_INFO "  Excludes" \
  " the \"angles\" and \"ambient\" parameters.\n" \
  "\n" \
  pst_lamp_spec_radius_HELP_INFO "  Excludes the" \
  " \"ambient\" and \"wall\" parameters.\n" \
  "\n" \
  pst_lamp_spec_wall_HELP_INFO "  Excludes the" \
  " \"radius\" and \"ambient\" parameters.\n" \
  "\n" \
  pst_lamp_spec_ambient_HELP_INFO "  Excludes the \"angles\"," \
  " \"direction\", \"wall\", and \"radius\" parameters.\n" \
  "\n" \
  pst_lamp_spec_power_HELP_INFO "  If not specified, assumes \"power 0 0 0\"."
  
/* Help/info for individual lamp parameters: */

#define pst_lamp_spec_angles_HELP \
  "angles {AZIMUTH} {ELEVATION}"

#define pst_lamp_spec_angles_INFO \
  "Specifies the spherical coordinates of the lamp" \
  " relative to the scene.  The {AZIMUTH} is measured" \
  " counterclockwise from the X axis. The {ELEVATION} is" \
  " measured from the XY plane towards the Z axis.  Both" \
  " angles should be given in degrees."

#define pst_lamp_spec_angles_HELP_INFO \
  "      " pst_lamp_spec_angles_HELP "\n" \
  "        " pst_lamp_spec_angles_INFO 

  
#define pst_lamp_spec_direction_HELP \
  "direction {DX} {DY} DZ}"

#define pst_lamp_spec_direction_INFO \
  "Specifies the direction of the lamp, in Cartesian" \
  " coordinates.  Only the direction of the vector is important," \
  " not its length."

#define pst_lamp_spec_direction_HELP_INFO \
  "      " pst_lamp_spec_direction_HELP "\n" \
  "        " pst_lamp_spec_direction_INFO 
  

#define pst_lamp_spec_radius_HELP \
  "radius {ANGRAD}"
 
#define pst_lamp_spec_radius_INFO \
  "Specifies the angular radius of the lamp, as seen from the" \
  " scene.  The value should be in radians" \
  " (not degrees!); it should be non-negative and at most" \
  " {PI/2}.  The default is \"radius 0\"," \
  " meaning a point-like light source."

#define pst_lamp_spec_radius_HELP_INFO \
  "      " pst_lamp_spec_radius_HELP "\n" \
  "        " pst_lamp_spec_radius_INFO 


#define pst_lamp_spec_wall_HELP \
  "wall"

#define pst_lamp_spec_wall_INFO \
  "Specifies a light source that covers a whole hemisphere" \
  " of directions with uniform apparent brightness; equivalent" \
  " to \"radius 1.5707963\"."

#define pst_lamp_spec_wall_HELP_INFO \
  "      " pst_lamp_spec_wall_HELP "\n" \
  "        " pst_lamp_spec_wall_INFO
 

#define pst_lamp_spec_ambient_HELP \
  "ambient"

#define pst_lamp_spec_ambient_INFO \
  "Specifies an isotropic `ambient light' field; equivalent" \
  " to \"radius 3.1415926\"."

#define pst_lamp_spec_ambient_HELP_INFO \
  "      " pst_lamp_spec_ambient_HELP "\n" \
  "        " pst_lamp_spec_ambient_INFO


#define pst_lamp_spec_power_HELP \
  "power " frgb_parse_color_HELP
  
#define pst_lamp_spec_power_INFO \
  "Specifies the intensity and color of the light produced" \
  " by the lamp. It consists of " frgb_parse_color_INFO "" 

#define pst_lamp_spec_power_HELP_INFO \
  "      " pst_lamp_spec_power_HELP "\n" \
  "        " pst_lamp_spec_power_INFO

#endif
