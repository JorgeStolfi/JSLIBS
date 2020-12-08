#ifndef pst_light_H
#define pst_light_H

/* pst_light.h -- light field modeling. */
/* Last edited on 2006-05-02 16:55:41 by stolfi */

#include <bool.h>
#include <r3.h>
#include <r3x3.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_lamp.h>
     
typedef struct pst_light_t 
  { pst_lamp_vec_t lmpv; /* List of lamps. */
  } pst_light_t;
  /* A light field. Presently the model is just a list of lamps
    (which may include ambient light fields or light walls). */

#define pst_light_model_INFO \
  "The light field model consists of one or more round" \
  " light sources (/lamps/), described below.  Each" \
  " lamp contributes a specific term to the total light" \
  " field reaching the scene from the outside."

pst_light_t pst_light_new(int NS, int NC);
  /* Creates a light field model with {NS} lamps (at least one)
    and {NC} color channels (possibly zero).  The
    {dir} of each lamp is set to {(0,0,0)} and the {crad} 
    is set to {+1} The {pwr} vectors are allocated
    and initialized with zeros. */

pst_light_t pst_light_from_lamps(pst_lamp_vec_t lmpv);
  /* Creates a light field model with the given list of lamps. */

pst_light_t pst_light_copy(pst_light_t *lht);
  /* Creates a new light field model with the same contents as {lht},
    but with newly allocted storage at all levels. Changing in any
    parameter of {lht} will not affect the copy, and vice-versa. */

vec_typedef(pst_light_vec_t,pst_light_vec,pst_light_t);
  /* A list of light fields. */

void pst_light_regularize_channels(pst_light_t *lht, int NC);
  /* Makes sure that all lamps have {NC} channels in their power
    vectors.  Vectors with zero elements are expanded {NC} elements
    and filled with {1/NS}, where {NS} is the number of 
    lamps.  Vectors with one channel are expanded to {NC} 
    by replicating the first element.  Vectors with {NC}
    elements are left untouched.  A vector with any other
    length causes the procedure to bomb out. */

void pst_light_ensure_one_lamp(pst_light_t *lht, double cmin, pst_lamp_t **src);
  /* Makes sure that the light field {lht} contains at least one
    source with cosine-of-radius {cmin} or greater, that is, with
    angular radius not exceeding {acos(cmin)}. If there is no such
    lamp in {lht}, adds one, with {src.crad == cmin}. In any case,
    returns in {*src} the address of the lamp with smaller radius,
    i.e. greatest {src.crad} --- which will be {cmin} or greater. */

void pst_light_alignment_matrix(r3_t *u, r3_t *v, r3x3_t *M);
  /* If {u} is not NULL and {*u} is not {(0,0,0)},sets {*M} to 
    a {3×3} rotation matrix that maps the direction vector {*u} to {(1,0,0)}.
    In that case, if if {v} is not NULL and {*v} is not {(0,0,0)},
    the matrix will map {*v} to the {X,Y} plane, with positive {Y}.
    If {u} is NULL, or {*u} is the null vector, set {M} to the
    identity matrix. */

/* CREATING ARRAYS OF LAMPS */

  /* The following procedures append sets of lamps to the argument
  vector {lmpv}. They assume {lmpv[0..*NSP-1]} are the existing
  lights. They expand {lmpv} as needed and increment {*NCP} by the
  number of lights added. Note that they may leave {lmpv.nel > *NCP},
  so the client must eventually trim the vector. */

void pst_light_add_uniform_array
  ( pst_lamp_vec_t *lmpv, 
    int *NSP, 
    int NA,
    r3_t *dir0, 
    r3_t *dir1, 
    double_vec_t *pwr,
    double crad
  );
  /* Adds to the lamp vector {lmpv} a set of of {NA} lamps
    uniformly distributed over the sphere of directions. Works only
    for certain values of {NA}: 1, 2, 3, 4, 6, 8, 12, 14, 20, 32, 60.
    
    When {NA=1}, the result is a single ambient lamp. When {NA=2}, the
    result is a pair of opposite wall lamps. When {NA=3}, the result
    is that pair plus an ambient lamp. When {NA} is 4 or more, the
    lamps are arranged at the corners of a regular or semi-regular
    polyhedron, and their radii are all less than {PI/2}.
    
    If {dir0} is not NULL, {*dir0} is not zero, and the first lamp of
    the array has a nonzero direction, the lamp array is rotated so
    that its first lamp is located in the direction {*dir0}. In that
    case, if {dir1} is not NULL, {*dir1} is not zero, and the second
    lamp exists and has nonzero direction, the array is then turned
    around {dir0} so that its second lamp lies on the {*dir0,*dir1}
    plane.
    
    If {pwr} is not NULL, the power of each lamp is set to a fresh 
    copy of {*pwr}; otherwise the powers are all set to empty vectors.
    
    If {crad} is finite, it is used as the {crad} of all lamps in the
    array. Otherwise, the angular radii are set to appropriate values
    depending on {NA}. */

void pst_light_add_single(pst_lamp_vec_t *lmpv, int *NSP, double crad);
  /* Adds a single lamp with  given {crad} and direction (+1,0,0). */
  
void pst_light_add_pair(pst_lamp_vec_t *lmpv, int *NSP, double crad);
  /* Adds a pair of lamps with given {crad} and directions (+1,0,0), (-1,0,0). */
  
void pst_light_add_tetra(pst_lamp_vec_t *lmpv, int *NSP, bool_t dual);
  /* Adds four lamps at the corners of a tetrahedron.  If {dual},
    uses the dual directions. */
  
void pst_light_add_octa(pst_lamp_vec_t *lmpv, int *NSP);
  /* Adds six lamps at the vertices of a octahedron, on the 
    coordinate axes. */
  
void pst_light_add_icosa(pst_lamp_vec_t *lmpv, int *NSP);
  /* Adds 12 lights at the vertices of an icosahedron. 
    The axes bisect pairs of opposite edges. */
  
void pst_light_add_dodeca(pst_lamp_vec_t *lmpv, int *NSP);
  /* Adds 20 lights at the vertices of a dodceahedron,
    dual to the icosahedron above. The axes bisect pairs 
    of opposite edges. */

/* COMMAND LINE PARSING */

pst_light_t pst_light_spec_parse(argparser_t *pp, bool_t next, int *NCP);
  /* Parses the description of a light field, consisting of a list of
    lamp descriptions. If {next} is true, the lamp list must start
    at the next command line argument; else looks for it anywhere in
    the command line. If no lamps are found, returns a model with zero
    lamps.  The syntax is described by {pst_light_spec_HELP}
    and {pst_light_spec_HELP_INFO}. */

#define pst_light_spec_HELP \
  "{ { -lamp | -array {NLAMPS} } \\\n" \
  "      " pst_lamp_spec_params_HELP " } ..."
  
#define pst_light_spec_HELP_INFO \
  pst_light_spec_lamp_HELP_INFO "\n" \
  pst_light_spec_array_HELP_INFO
  

/* Help/info for single lamp specs: */

#define pst_light_spec_lamp_HELP \
  "-lamp {LAMP_PARAMETERS}"

#define pst_light_spec_lamp_INFO \
  "Introduces an additional lamp.  Should be followed" \
  " by the lamp's parameters.  This keyword must" \
  " appear at least once."
  
#define pst_light_spec_lamp_HELP_INFO \
  "  " pst_light_spec_lamp_HELP \
  "    " pst_light_spec_lamp_INFO


/* Help/info for lamp array specs: */

#define pst_light_spec_array_HELP \
  "-array {NLAMPS} {LAMP_PARAMETERS}"
  
#define pst_light_spec_array_INFO \
  "Adds to the lighting model an array of {NLAMPS} lamps,"\
  " with the given parameters"\
  " uniformly distributed over the direction sphere.  Currenly"\
  " implemented only for a few values of {NLAMPS}, including"\
  " 1,2,3,4,6,8,12,14,20,32."

#define pst_light_spec_array_HELP_INFO \
  "  " pst_light_spec_array_HELP \
  "    " pst_light_spec_array_INFO

/* Help/info for lamp parameter specs (single or arrays): */

#define pst_light_spec_lamp_params_HELP_INFO \
  pst_lamp_spec_angles_HELP_INFO "  Each lamp description must have" \
  " either \"angles\" or \"direction\" (not both).  For a lamp" \
  " array, it specifies the direction of the first lamp.\n" \
  "\n" \
  pst_lamp_spec_direction_HELP_INFO "\n" \
  "\n" \
  pst_lamp_spec_radius_HELP_INFO "  The default is \"radius 0\"," \
  " meaning a point-like light source.\n" \
  "\n" \
  pst_lamp_spec_power_HELP_INFO "  The default is" \
  " \"power 1 / {NS}\", where {NS} is the number of lamps" \
  " specified for the scene.  Thus, if this parameter is" \
  " omitted for all lamps, the apparent intensity of a matte white" \
  " surface will be 1.0 or less."

void pst_light_spec_write(FILE *wr, pst_light_t *lht);
  /* Writes to {wr} the parameters of the light field {lht}, in a format
    compatible with {pst_light_spec_parse}. */

#endif
