#ifndef pst_geom_sphere_H
#define pst_geom_sphere_H

/* pst_geom.h -- geometric tools for images of spheres. */
/* Last edited on 2025-01-21 19:45:23 by stolfi */

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <argparser.h>
#include <ellipse_crs.h>

#include <pst_basic.h>

/* UTILITIES */

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double rad);
  /* Clips the direction vector {udir} to a spherical cap with center
    on the direction vector {sdir} and angular radius {rad. Namely, if
    the angle between {udir} and {sdir} is greater than {rad} radians,
    sets {udir} to the unit vector that lies on the shortest arc
    between the two, at {rad} radians from {sdir}. */

/* GEOMETRY OF SPHERE PROJECTIONS */

  /* A {ellipse_crs_t} describes the perspective projection of a
    spherical object. It is assumed that the sphere radius is small
    compared to the camera-to-object distance, so that the perspective
    projection is well approximated by a parallel oblique projection.
    This assumption should be harmless for a small sphere taken at
    typical camera geometries.
    
    The sphere's projection is then an ellipse, and the center of the
    sphere projects onto the center {ctr} of the ellipse. The length
    of the smallest radius (minor semidiameter) of the ellipse is
    {rad}, while the longest radius (major semidiameter) is parallel
    to {str} and has length {rad+len}, where {len} is the Euclidean
    length of {str}.

    All dimensions and coordinates are in pixels. The origin is
    assumed to be at the lower left corner of the image, with the
    Y axis pointing up. */

r3_t pst_geom_sphere_compute_normal(r2_t *uv);
  /* Computes the outwards normal vector to the sphere's surface at
    the point that projects onto the point {uv}. Assumes a {U,V,W}
    coordinate system where the sphere has radius 1 and center at the
    origin, and {W} points towards the camera. */

void pst_geom_sphere_view_matrices
  ( ellipse_crs_t *geo, /* Geometry of sphere's projection. */
    r3x3_t *xym_to_uvm, 
    r3x3_t *uvw_to_xyz
  );
  /* Computes the projection matrices {xym_to_uvm} and
    {uvw_to_xyz} needed by {pst_normal_map_from_proc}.
    
    The first matrix, {xym_to_uvm}, transforms the homogeneous
    coordinates {[x,y,1]} in the Image Coordinate System (ICS) to
    coordinates {[u,v,1]} on a two-dimensional homogeneous Object
    Coordinate System (OCS) where {[0,0,1]} is the sphere's center,
    {[1,0,1]} is the tip of the minor semidiameter, and {[0,1,1]} is
    the point of the sphere that projects onto the tip of the major
    semidiameter.  These points define a plane that goes through 
    the sphere's center, perpendicular to the line from the center
    to the camera.
    
    The OCS system has a third coordinate {w} such that {u=0,v=0,w=1}
    is the point on the sphere's surface that is closest to the
    camera. It follows that the {w} coordinate of the visible point of
    the sphere with given {u,v} coordinates is {w = sqrt(1-u^1-v^2)}.
    
    The second matrix, {uvw_to_xyz}, is a linear map that converts
    vectors from the OCS coordinates back to ICS coordinates {x,y,z},
    where {z} is displacement perpendicular to the image plane (in
    pixels). */

/* COMMAND LINE ARGUMENT PARSING */
  
ellipse_crs_t *pst_geom_sphere_parse
  ( argparser_t *pp,
    bool_t next,
    r2_t *ctrdef, 
    double *ctrAdj, 
    double *radAdj, 
    double *strAdj,
    double *adjStep
  );
  /* Parses from the command line for tge keyword "-sphere" possibly
    followed the geometric parameters (center, radius and stretch
    vector) of the projection of a sphere, and returns them packaged
    as a {ellipse_crs_t}.
    
    If {next} is TRUE, looks for "-sphere" only at the next command
    line argument, else looks for it anywhere in the command line.
    
    The center defaults to {ctrdef}; if {ctrdef} is NULL, the center
    is mandatory. The radius is mandatory. The stretch vector defaults
    to {(0,0)}.
    
    If {ctrAdj} is not NULL, the center coordinates may be followed by
    "adjust {AMOUNT}". In that case, the procedure sets
    {*ctrAdj=AMOUNT} if the option is present, and {*ctrAdj=0}
    otherwise. The {radAdj} and {strAdj} provide the same alternative
    for the radius and stretch vector, respectively.
    
    If {adjStep} is not NULL, the procedure also allows "step {STEP}"
    after all parameter and adjustment specs. In that case, the
    procedure sets {*adjStep=STEP} if that option is present, and to
    {*adjStep=+INF} otherwise.
    
    The syntax is described by {pst_geom_sphere_XXX_HELP} and
    {pst_geom_sphere_XXX_INFO}, where {XXX} is {center}, {radius}, or
    {stretch}. All the parameters of the same sphere must appear
    together in the command line. See {argparser.h} for an explanation
    of the {pp} parameter.  */

#define pst_geom_sphere_center_HELP \
  "center {CTRX} {CTRY}"
  
#define pst_geom_sphere_center_INFO \
  "Specifies the center of the sphere's projetion, in pixels" \
  " from the bottom left corner."

#define pst_geom_sphere_center_HELP_INFO \
  "      " pst_geom_sphere_center_HELP "\n" \
  "        " pst_geom_sphere_center_INFO 

#define pst_geom_sphere_radius_HELP \
  "radius {RAD}"
  
#define pst_geom_sphere_radius_INFO \
  "Specifies the radius of the sphere's projection on" \
  " the image, in pixels.  If the gauge's projection" \
  " is stretched to an ellipse, {RAD} is the smallest radius," \
  " i.e. the minor semidiameter."

#define pst_geom_sphere_radius_HELP_INFO \
  "      " pst_geom_sphere_radius_HELP "\n" \
  "        " pst_geom_sphere_radius_INFO 
  
#define pst_geom_sphere_stretch_HELP \
  "stretch {STRX} {STRY}"

#define pst_geom_sphere_stretch_INFO \
  "Specifies the direction and and amount of perspective stretching" \
  " of the gauge's elliptical image.  The longest radius (the major" \
  " semidiameter) of the ellipse will parallel to the vector" \
  " {(STRX,STRY)}, and its length will be {RAD} plus the length of that" \
  " vector."

#define pst_geom_sphere_stretch_HELP_INFO \
  "      " pst_geom_sphere_stretch_HELP "\n" \
  "        " pst_geom_sphere_stretch_INFO 
  
#define pst_geom_sphere_center_adjust_HELP \
  "adjust {CTR_ADJUST}"

#define pst_geom_sphere_radius_adjust_HELP \
  "adjust {RAD_ADJUST}"

#define pst_geom_sphere_stretch_adjust_HELP \
  "adjust {STR_ADJUST}"

#define pst_geom_sphere_all_adjust_INFO \
  "These options can be specified after the values of the" \
  " \"center\", \"radius\", and/or \"stretch\" parameters," \
  " respectively.  Each specifies the maximum amount by which" \
  " the preceding parameter may be adjusted by the program.  The" \
  " {CTR_ADJUST} amount applies to each center coordinate" \
  " independently.  The {STR_ADJUST} amount applies only to" \
  " the length of the stretch vector {(STRX,STRY)}; its direction" \
  " is preserved."

#define pst_geom_sphere_all_adjust_HELP_INFO \
  "      " pst_geom_sphere_center_adjust_HELP "\n" \
  "      " pst_geom_sphere_radius_adjust_HELP "\n" \
  "      " pst_geom_sphere_stretch_adjust_HELP "\n" \
  "        " pst_geom_sphere_all_adjust_INFO

#define pst_geom_sphere_adjust_step_HELP \
  "step {ADJ_STEP}"

#define pst_geom_sphere_adjust_step_INFO \
  "This option may be specified once, after all values of" \
  " \"center\", \"radius\", and/or \"stretch\" parameters" \
  " and their \"adjust\" amounts (if any).  It specifies the" \
  " ideal granularity by which all parameters may be adjusted by" \
  " the program.  If no \"adjust\" was specified, this parameter" \
  " is meaningless.  The actual step used for each parameter may" \
  " be somewhat different, depending on the adjustment inerval, the" \
  " adjustment algorithm used, imposed time limits, etc.."

#define pst_geom_sphere_adjust_step_HELP_INFO \
  "      " pst_geom_sphere_adjust_step_HELP "\n" \
  "        " pst_geom_sphere_adjust_step_INFO
 
void pst_geom_sphere_write(FILE *wr, ellipse_crs_t *geo);
  /* Writes to {wr} the geometric parameters of , in a format compatible
    with {pst_geom_sphere_parse}. */


#endif
