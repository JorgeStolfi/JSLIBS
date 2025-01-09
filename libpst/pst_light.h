#ifndef pst_light_H
#define pst_light_H

/* pst_light.h -- light field modeling. */
/* Last edited on 2025-01-05 10:17:27 by stolfi */

#include <bool.h>
#include <r3.h>
#include <r3x3.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_lamp.h>
     
typedef struct pst_light_t 
  { uint32_t NL;         /* Number of lamps. */
    pst_lamp_vec_t lmpv; /* List of lamps. */
  } pst_light_t;
  /* A light field. Presently the model is just a list of lamps
    { lmpv.e[0..NL-1]} (which may include ambient light fields 
    or light walls). */

#define pst_light_model_INFO \
  "The light field model consists of one or more round" \
  " light sources (/lamps/), described below.  Each" \
  " lamp contributes a specific term to the total light" \
  " field reaching the scene from the outside."

pst_light_t *pst_light_new(uint32_t NL);
  /* Creates a light field model with {NL} lamps (at least one).  The
    {dir} of each lamp is set to {(NAN,NAN,NAN)},
    the power {pwr} is set to {(0,0,0)}, and the {crad} 
    is set to {+1}. */

pst_light_t *pst_light_from_lamps(uint32_t NL, pst_lamp_vec_t lmpv);
  /* Creates a light field model with lamps
    {lmpv.e[0..NL-1]}. */

pst_light_t *pst_light_copy(pst_light_t *lht);
  /* Creates a new light field model with the same contents as {lht},
    but with newly allocted storage at all levels. Changing in any
    parameter of {lht} will not affect the copy, and vice-versa. */

void pst_light_ensure_one_lamp(pst_light_t *lht, double cmin, pst_lamp_t **src);
  /* Makes sure that the light field {lht} contains at least one
    source with cosine-of-radius {cmin} or greater, that is, with
    angular radius not exceeding {acos(cmin)}. If there is no such
    lamp in {lht}, adds one, with {src.crad == cmin}. In any case,
    returns in {*src} the address of the lamp with smaller radius,
    i.e. greatest {src.crad} --- which will be {cmin} or greater. */

/* CREATING ARRAYS OF LAMPS */

  /* The following procedures append sets of lamps to the argument
  vector {lmpv}. They assume {lmpv[0..*NSP-1]} are the existing
  lights. They expand {lmpv} as needed and increment {*NCP} by the
  number of lights added. Note that they may leave {lmpv.nel > *NCP},
  so the client must eventually trim the vector. */

void pst_light_add_uniform_array
  ( pst_light_t *lht, 
    uint32_t NA,
    r3_t *dir0, 
    r3_t *dir1, 
    frgb_t *pwr,
    double crad
  );
  /* Adds to the lamp vector {lmpv} a set of of {NA} lamps
    uniformly distributed over the sphere of directions. Works only
    for certain values of {NA}: 1, 2, 3, 4, 6, 8, 12, 14, 20, 32, 60.
    
    When {NA=1}, the result is a single ambient lamp. When {NA=2}, the
    result is a pair of opposite wall lamps. When {NA=3}, the result
    is that pair plus an ambient lamp. When {NA} is 4 or more, the
    lamps are arranged at the corners of a regular or semi-regular
    polyhedron, and their agular radii are all less than {PI/2}.
    
    If {dir0} is not NULL, and {*dir0} is a valid direction, and the
    first lamp of the array has a nonzero direction, the lamp array is
    rotated so that its first lamp is located in the direction {*dir0}.
    In that case, if {dir1} is not NULL, and {*dir1} is a valid
    direction, and the second lamp exists and has nonzero direction, the
    array is then turned around {dir0} so that its second lamp lies on
    the {*dir0,*dir1} plane.
    
    If {pwr} is not NULL, the power of each lamp is set to {*pwr};
    otherwise the powers are all set to {(0,0,0)}.
    
    If {crad} is finite, it is used as the {crad} of all lamps in the
    array. Otherwise, the angular radii are set to appropriate values
    depending on {NA}. */
  
void pst_light_add_single(pst_light_t *lht, double crad);
  /* Adds to {lmpv} a single lamp with  given {crad} and direction (+1,0,0). */

void pst_light_add_pair(pst_light_t *lht, double crad);
  /* Adds a pair of lamps with given {crad} and directions (+1,0,0), (-1,0,0). */
  
void pst_light_add_tetra(pst_light_t *lht, bool_t dual);
  /* Adds four lamps at the corners of a tetrahedron.  If {dual},
    uses the dual directions. */
  
void pst_light_add_octa(pst_light_t *lht);
  /* Adds six lamps at the vertices of a octahedron, on the 
    coordinate axes. */
  
void pst_light_add_icosa(pst_light_t *lht);
  /* Adds 12 lights at the vertices of an icosahedron. 
    The axes bisect pairs of opposite edges. */
  
void pst_light_add_dodeca(pst_light_t *lht);
  /* Adds 20 lights at the vertices of a dodceahedron,
    dual to the icosahedron above. The axes bisect pairs 
    of opposite edges. */

/* SHADING */

frgb_t pst_light_shading(pst_light_t *lht, r3_t *nrm);
  /* Computes the shading factor 
    from the light field {lht} for a surface with normal {nrm}.
    That is the amount of light per channel
    falling per unit area of the surface.  If the surface is 
    Lambertian with white intrinsic color (albedo), that will
    be its apparent color, viewed from any direction.
    
    The result is the sum of the intensity of each lamp in {lht}
    times its geometric factor for normal {nrm} (See {pst_lamp_geom_factor}).
    Depending on the lamp intensities, the final value may be greater than 1. */

/* LIGHTING SPEC I/O */

void pst_light_spec_write(FILE *wr, pst_light_t *lht);
  /* Writes to {wr} the parameters of the light field {lht}, in a format
    compatible with {pst_light_spec_parse}. */

/* COMMAND LINE PARSING */

pst_light_t *pst_light_spec_parse(argparser_t *pp, bool_t next);
  /* Parses the description of a light field, consisting of a list of
    lamp descriptions. If {next} is true, the lamp list must start
    at the next command line argument; else looks for it anywhere in
    the command line. If no lamps are found, returns a model with zero
    lamps.  The syntax is described by {pst_light_spec_HELP}
    and {pst_light_spec_HELP_INFO}. */

/* HELP/INFO FOR MULTIPLE LAMPS AND/OR LAMP ARRAYS */

#define pst_light_spec_HELP \
  "{ { lamp | array {NLAMPS} | wall | ambient } {SOURCE_PARAMETERS} }..."
 
#define pst_light_spec_HELP_INFO \
  "The light field is described as a list of one or more light" \
  " sources.  "  pst_lamp_spec_HELP_INFO(pst_light_array_spec_HELP_INFO,pst_light_array_direction_INFO,pst_light_array_radius_INFO,pst_light_array_power_INFO) ""
  
#define pst_light_array_spec_HELP_INFO \
  "      array {NLAMPS}\\\n" \
  "        Specifies an array of {NLAMPS} light sources with the given parameters," \
  " whose directions are"\
  " uniformly distributed over the direction sphere.  Currenly"\
  " implemented only for a few values of {NLAMPS}, including"\
  " 1,2,3,4,6,8,12,14,20,32."
  
#define pst_light_array_direction_INFO \
  "  For an \"array\" source, the direction is that of the first lamp in the array."
  
#define pst_light_array_radius_INFO \
  "  For an \"array\" source, if \"radius\" is not specified, the" \
  " radii of the lamps in the array are adjusted so that each lamp" \
  " extends to about halfway of the distance do the nearest lamp in the array."
  
#define pst_light_array_power_INFO \
  "  Whenever \"power\" is not specified, the default in each color channel will be" \
  " \"power {1/NL}\", where {NL} is the total number of lamps" \
  " specified for the scene.  Thus, if this parameter is" \
  " omitted for all lamps, the apparent intensity of a matte white" \
  " surface will be 1.0 or less.   If \"power {POWER}\" is specified" \
  " for an \"array\" with {NLAMPS} lamps, the power of each of" \
  " those lamps will be {POWER/NLAMPS}."
  
#endif
