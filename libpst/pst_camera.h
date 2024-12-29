#ifndef pst_camera_H
#define pst_camera_H

/* pst_camera.h -- camera parameters and tools therefor. */
/* Last edited on 2024-12-28 00:42:02 by stolfi */

#include <r2.h>
#include <r3.h> 
#include <hr2.h> 
#include <hr3.h> 
#include <argparser.h>

#include <pst_basic.h>

/* !!! Add radial distortion {kappa} to the parameters. !!! */
/* !!! Allow for non-square pixels. !!! */

/* MATHEMATICAL CAMERA
  
  The procedures in this interface deal with an idealized /camera/,
  which is essentially a mapping that projects any point {P} in some
  3D space to a point {P'} on a certain /image plane/ {IPL}.
  
  For every digital photo or video frame, there is a unique camera
  that describes how points in physical space were mapped to 
  points on the image.
  
  PERSPECTIVE PROJECTION
  
  At present, this interface only supports cameras that perform a
  restricted kind of mappint, a /conical/ (or /perspective/)
  projection onto a plane.
  
  We use a /nominal coordinate system/ (NCS) for 3D space whose
  coordinates {(X,Y,Z)} are measured in the same units (pixels, mm,
  whatever), and such that the image plane {IPL} is the plane with
  equation {Z==0}.
  
  The perspective mapping is entirely determined by its /viewpoint/
  {O}, a point of 3-space (possibly infinite) that has positive 
  {z} coordinate. From {O} we define:
  
    the /optical axis/ as the line through {O} 
    parallel to the {Z} axis;
    
    the /optical center/ {Q} as being the point
    where the optical axis meets the image plane; 
    
    the /focal length/ {F = O.Z = O.z/O.m}; 
    
    the (/angular/) /spread/ {G} of the camera, defined
    as the reciprocal {1/F};  and
    
    the /view plane/ {VPL} as the plane with equation {Z==F}, parallel
    to {IPL} and passing through {O}.
    
  The projection of a point {P} of three-space is the point {P'} of
  {IPL} such that {O}, {P} and {P'} are collinear and {O} is not
  between {P} and {P'}. In particular, if {P} happens to lie on {IPL},
  then it coincides with its projection {P'}. This point exists (and is
  unique) if and only if {P} lies strictly below the view plane {VPL};
  that is, if {P.Z<F}. This condition defines the /visible half-space/
  of the camera.
  
  POINTS AT INFINITY AND HOMOGNEOUS COORDINATES
  
  To refer to points at infinity, we use the concept of /oriented
  projective space/ and /signed homogeneous coordinates/ {[m,x,y,z]},
  as defined in [1]. Namely, the {m} coordinate is the /weight/ or
  scale factor (here always non-negative), and the corresponding
  Cartesian coordinates are {X = x/m}, {Y = y/m}, and {Z = z/m}.
  
  The geometric definition of the perspective map still works if {O}
  or {P} (but not both) are points at infinity, provided that {O} lies
  above the {IPL}, {P} lies below the {VPL}, and neither point is in
  the yonder realm.
  
  When {O} is at infinity, the camera actually performs a
  /cylindrical/ projection, with parallel projection lines. The focal
  length {F} is {+INF}, the spread {G} is zero, and the view plane
  {VPL} is is at infinity. In particular, if {O} is the zenith point
  {[0,0,0,1]}, at infinite distance on the {Z} axis, the map is the
  perpendicular projection onto the image plane.
  
  In general, the optical center {Q} has homogeneous coordinates
  {[O.m,O.x,O.y,0]}, and is a finite point if and only if {O} is
  finite. However, for the perpendincular projection (with
  {O=[0,0,0,1]}) the optical axis and optical center are undefined.
    
  PROJECTION FORMULAS
  
  The projected point {P'} can be represented by its two coordinates
  {(P'.X,P'.Y)} only (the /image coordinates/), since its {Z}
  coordinate is always zero.
  
  If {O} and {P} are finite, then {P' == (O.X + m*(P.X - O.X), O.Y +
  m*(P.Y - O.Y)} where {m = F/(F - P.Z)}.
  
  If {P} is the point at infinity with homogeneous coordinates
  {[0,O.x,O.y,O.z]}, then {P'= (P.X - m*D.x, P.Y - m*D.y)} where {m =
  P.Z/D.z}.
  
  CHOOSING THE NOMINAL COORDINATE SYSTEM
  
  In computer graphics is is customary to choose the image coordinate
  system so that the top left corner of the stored digital image has
  image coordinates {(0,0)}, with the {X} axis pointing to the right
  and the {Y} axis pointing down. Some people prefer to put the orgin
  at the lower left corner, with the {Y} axis pointing up; or at the
  center, with the {Y} axis pointing either way. Coordinates are
  usually measured in pixels, so that each pixel is a square of side 1
  whose corners have integer coordinates.
  
  The NCS can be chosen as any of these image systems, complemented
  with a {Z} axis that points away from the image plane, towards the
  viewing side. The only restriction (for now) is that the same
  measurement unit must be used along the three axes, so that the
  Euclidean distance between two points can be computed by the usual
  Pythagorean formula. Thus, if {X} and {Y} are measured in pixels,
  the {Z} axis must be measured in pixels, too.
  
  POSITION OF THE OPTICAL CENTER
  
  The optical axis of a photo or video camera usually goes through the
  center of the (full, uncropped) images it produces. Therefore, in a
  digital image with {NX} by {NY} pixels, the optical center is the
  point {(NX/2,NY/2,F)}.
  
  Note that this is no longer the case when working with an image that
  has been cut out from a larger one.  In that case, the  optical center
  may fall anywhere, even outside the image.
  
  DETERMINING THE FOCAL LENGTH
  
  The focal length {F = O.Z} is the perpendicular distance from the
  viewpoint {O} to the image plane, in the same units used for the {X}
  and {Y} image coordinates. To determine the focal length {F} of a
  camera, take a picture of a ruler of known length {L}, perpendicular
  to the optical axis, placed at a known distance {D} from the
  viewpoint. If the projection of the ruler on the image has length
  {T}, then {F} is {T*(D/L)}. This value of {F} can be used for all
  pictures taken with the same camera, as longa s the zoom setting
  (optical or digital) and image size are not changed.
  
  Varying the lens zoom of a physical camera corresponds to changing
  the focal length {F = O.Z} while preserving the image coordinates
  {O.X,O.Y} of the optical center. Zooming `in' increases {F}, zooming
  `out' decreases it. As the focal length {F} tends to infinity, with
  {O.X,O.Y} fixed, the perspective mapping tends to the cylindrical
  projection, orthogonal to the image plane.
  
  When an image is reduced or enlarged, the Cartesian coordinates of
  the viewpoint (i.e., the focal length and optical center) of the
  associated camera must be multiplied by the magnification factor.
  
  CAMERA SPREAD
  
  The spread {G = 1/F} of a camera is roughly the angle in radians
  subtended by one image pixel near the optical axis. A
  camera with cylindrical projection has spread 0. */

typedef struct pst_camera_t 
  { hr3_point_t O;  /* The view point of the camera. */
  } pst_camera_t;
  /* Representation of a camera. The viewpoint is given as a quadruple
    {O} of signed homogeneous coordinates {[m,x,y,z]}, where {m}
    is the /weight/ (scale factor). The {z} coordinate is always
    positive, and {m} is always non-negative. If {m} is zero, the
    viewpoint is at infinity in the direction of the vector {(x,y,z)}.
    If {m} is positive, the cartesian coordinates {(O.X,O.Y,O.Z)} of
    the viewpoint are {(x/m,y/m,z/m)}. */

double pst_camera_focal_length(hr3_point_t *O);
  /* The focal length of camera {C}. The result is always positive,
    and {+INF} if the camera is at infinity. */
    
double pst_camera_spread(hr3_point_t *O);
  /* The spread of camera {C}.  The result is
    always finite and non-negative, and is 0 
    if the camera is at infinity. */
    
hr2_point_t pst_camera_center(hr3_point_t *O);
  /* The homogeneous coordinates of the optical center of camera {C}.
    Returns {[0,0,0]} for an orthogonal projection, and
    a point at infinity for an oblique cylintrical projection. */
    
r2_t pst_camera_center_cartesian(hr3_point_t *O);
  /* Same as {pst_camera_center}, but returns the Cartesian coordinates
    of {Q}.  The result is meaningless if the viewpoint is infinite. */

hr3_point_t pst_camera_viewpoint_from_center_spread(r2_t *Q, double G);
  /* Computes the viewpoint of a camera given the optical
    center {Q} and the spread {G}. */

r3x3_t pst_camera_normal_correction_matrix(r3_t *dir);
  /* The orthonormal rotation matrix that converts a surface normal vector that is
    relative to view direction {dir} to absoute normal direction. */
 
#define pst_camera_max_angle 80.0
  /* Max angle (in degrees) allowed at the viewpoint from the optical
    axis to any pixel in the image. */

double pst_camera_min_focal_length(r2_t *Q, int32_t NX, int32_t NY);
  /* A reasonable lower limit for the focal length {F} of a camera
    which took an image with {NX} columns and {NY} rows, whose optical
    center has image coordinates {Q}. The result is such that no pixel
    within the image will be more than {pst_camera_max_angle} degrees
    from the optical axis. If {Q} is NULL, assumes that the optical
    center is the center {(NX/2,NY/2)} of the image. */
    
/* CAMERA OPS FOR MULTISCALE */

pst_camera_t pst_camera_shrink(pst_camera_t *C, int32_t dx, int32_t dy, int32_t nw);
  /* Given a camera description appropriate for an image {A},
    returns a camera description for a reduced image {B} 
    obtained by {float_image_mscale_shrink(A, M, dx, dy, nw)}. */

pst_camera_t pst_camera_expand(pst_camera_t *C, int32_t dx, int32_t dy, int32_t nw);
  /* Inverse of {pst_camera_shrink}: given a camera description
    appropriate for the reduced image {B}, returns a camera
    description for the original image {A}. */

/* PARSING CAMERA SPECS FROM THE COMMAND LINE

  Any camera, finite or infinite, can be completely specified on the
  command line by giving the signed homogeneous coordinates its
  viewpoint {O}.
  
  Alternatively, the user may prefer to specify separately the image
  coordinates of the optical center {Q = (O.X,O.Y)} and the spread
  {G = 1/F}. These are converted automatically to the equivalent
  viewpoint coordinates. Specifying {G == 0} with any {Q} gives an
  orthogonal projection. A cylindrical oblique projection cannot be
  specified in this way.
  
  At present, the adjustment amounts can be specified only in 
  the center/spread model. */

void pst_camera_args_parse
  ( argparser_t *pp, /* The command line parsing state. */
    pst_camera_t *C, /* (OUT) The camera parameters. */
    double *OAjd,    /* (OUT) The viewpoint adjustment amount, or NULL. */ 
    double *QAdj,    /* (OUT) The optical center adjustment amount, or NULL. */ 
    double *GAdj     /* (OUT) The spread adjustment amount, or NULL. */ 
  );
  /* Parses from the command line the essential camera parameters 
    and packs them as a {pst_camera_t}.
    
    The syntax is described by {pst_camera_args_XXX_HELP} and
    {pst_camera_args_XXX_INFO}, where {XXX} is {viewpoint}, {center} or
    {spread}. All the parameters of the same camera must appear
    together in the command line. See {argparser.h} for an explanation
    of the {pp} parameter.
    
    The "viewpoint" option is incompatible with the "center" and
    "spread" options. Any cordinates of the viewpoint that are not
    specified with these options is set to {NAN}. This applies also
    when one of these keywords is present but is not followed
    immediately by a numeric value.
    
    If {OAjd} is not NULL, the viewpoint coordinates may be
    followed by the keyword "adjust" and a float value {VPT_ADJUST},
    which is stored in {*OAjd}. If the "adjust" keyword is not
    present, the procedure stores 0 in {*OAjd}. If {OAjd} is
    NULL, the keyword "adjust" is not allowed in the "viewpoint"
    spec.
    
    The {QAdj} and {GAdj} parameters play the same role
    for the "center" and "spread" models. */

/* PRINTOUT */

void pst_camera_print(FILE *wr, pst_camera_t *C, char *fmt);
  /* Prints {C} to {wr}, using format {fmt} for the 
    homogeneous coordinates of the viewpoint. */

void pst_camera_args_print(FILE *wr, pst_camera_t *C, char *fmtp, char *fmtr);
  /* Prints the parameters of {C}, in the center/spread model if
    possible, in a format compatible with {pst_camera_args_parse}.
    Uses {fmtp} for coordinates and {fmtr} for the spread. */

void pst_camera_args_adjust_print_viewpoint
  ( FILE *wr, 
    pst_camera_t *C, 
    double OAjd,
    char *fmtp,
    char *fmtd
  );
  /* Prints the parameters of {C} (and their adjustments, if any) in
    the "viewpoint" model, in a format compatible with
    {pst_camera_args_parse}. Uses {fmtp} for Cartesian
    coordinates of finite points, {fmtd} for coordinates
    of unit direction vectors. */

void pst_camera_args_adjust_print_center_spread
  ( FILE *wr, 
    pst_camera_t *C,
    double QAdj,
    double GAdj,
    char *fmtp, 
    char *fmtr
  );
  /* Prints the parameters of {C} (and their adjustments, if any) in
    the "center"/"spread" model, in a format compatible with
    {pst_camera_args_parse}. Uses {fmtp} for coordinates and {fmtr}
    for the spread. */

/* ---------------------------------------------------------------------- */
/* Documentation for "viewpoint": */

#define pst_camera_args_viewpoint_HELP \
  "viewpoint {MVPT} {XVPT} {YVPT} {ZVPT}"

#define pst_camera_args_viewpoint_adjust_HELP \
  "viewpoint [ {MVPT} {XVPT} {YVPT} {ZVPT} ] [ adjust {VPT_ADJUST} ]"

  
#define pst_camera_args_viewpoint_INFO \
  "This option specifies that the camera's viewpoint" \
  " has the given homogeneous coordinates, which may not" \
  " be all zero.  The {ZVPT} coordinate" \
  " must be positive.  The {MVPT} coordinate is" \
  " the weight (scale factor) and must be" \
  " non-negative.  If {MVPT} is positive, the Cartesian" \
  " coordinates of the viewpoint are" \
  " {(XVPT/MVPT, YVPT/MVPT, ZVPT/MVPT)}.  If {MVPT} is" \
  " zero, the viewpoint is at infinity in the" \
  " direction of the vector {(XVPT, YVPT, ZVPT)}."

#define pst_camera_args_viewpoint_adjust_INFO \
  pst_camera_args_viewpoint_INFO \
  "  If \"adjust\" is present, {VPT_ADJUST} specifies" \
  " the maximum amount by which the viewpoint may be adjusted.  If" \
  " the viewpoint is finite, each Cartesian coordinate" \
  " may be adjusted by that amount. If the viewpoint" \
  " is at infinity, its direction may fe adjusted" \
  " by {VPT_ADJUST} radians."

#define pst_camera_args_viewpoint_HELP_INFO \
  "      " pst_camera_args_viewpoint_HELP "\n" \
  "        " pst_camera_args_viewpoint_INFO 

#define pst_camera_args_viewpoint_adjust_HELP_INFO \
  "      " pst_camera_args_viewpoint_adjust_HELP "\n" \
  "        " pst_camera_args_viewpoint_adjust_INFO

/* ---------------------------------------------------------------------- */
/* Documentation for "center": */

#define pst_camera_args_center_HELP \
  "center {XOPC} {YOPC}"

#define pst_camera_args_center_adjust_HELP \
  "center [ {XOPC} {YOPC} ] [ adjust {OPC_ADJUST} ]"

  
#define pst_camera_args_center_INFO \
  "This option specifies that the camera's" \
  " optical axis goes through the point" \
  " with image coordinates {(XOPC,YOPC)}.  The" \
  " optical center may lie outside the" \
  " image's domain; this is often the" \
  " case when the image was" \
  " cropped from a larger photo."

#define pst_camera_args_center_adjust_INFO \
  pst_camera_args_center_INFO \
  "  If \"adjust\" is present, {OPC_ADJUST} specifies" \
  " the maximum amount by which {XOPC} and/or {YOPC}" \
  " may be adjusted."

#define pst_camera_args_center_HELP_INFO \
  "      " pst_camera_args_center_HELP "\n" \
  "        " pst_camera_args_center_INFO 

#define pst_camera_args_center_adjust_HELP_INFO \
  "      " pst_camera_args_center_adjust_HELP "\n" \
  "        " pst_camera_args_center_adjust_INFO

/* ---------------------------------------------------------------------- */
/* Documentation for "spread": */

#define pst_camera_args_spread_HELP \
  "spread {SPR}"

#define pst_camera_args_spread_adjust_HELP \
  "spread [ {SPR} ] [ adjust {SPR_ADJUST} ]"
  
  
#define pst_camera_args_spread_INFO \
  "This option specifies the" \
  " spread (reciprocal of the focal length) of the camera." \
  " This is approximately the actual size (in meters) of an object" \
  " that is 1 meter away from the camera and whose" \
  " image is one pixel wide.  If {SPR} is zero, the camera" \
  " performs an orthogonal projection."

#define pst_camera_args_spread_adjust_INFO \
  pst_camera_args_spread_INFO \
  "  If \"adjust\" is present, {SPR_ADJUST} specifies the maximum amount" \
  " by which the spread {SPR} may be adjusted."

#define pst_camera_args_spread_HELP_INFO \
  "      " pst_camera_args_spread_HELP "\n" \
  "        " pst_camera_args_spread_INFO 

#define pst_camera_args_spread_adjust_HELP_INFO \
  "      " pst_camera_args_spread_adjust_HELP "\n" \
  "        " pst_camera_args_spread_adjust_INFO

/* REFERENCES 

  [1] "Oriented Projective Geomtry: A Framework for Geometric Computations."
  By Jorge Stolfi.  Academic Press (1989). */

#endif
