#ifndef pst_sphere_H
#define pst_sphere_H

/* pst_sphere.h -- geometric tools for images of spheres. */
/* Last edited on 2024-11-04 07:26:13 by stolfi */

#include <ellipse_crs.h>
#include <r2.h>
#include <r3.h>
#include <hr3.h> 
#include <hr3_pmap.h> 
#include <argparser.h>

#include <pst_basic.h>
#include <pst_camera.h>

/* TOOLS FOR IMAGES OF SPHERES
  
  The procedures in this interface deal with spherical objects
  as they appear in digital images. 
  
  In this interface we use the 3D /nominal coordinate system/ (NCS)
  which is tied to the image as described in {pst_camera.h}. All three
  NCS coordinates {(X,Y,Z)} are measured in the same unit. The /image
  plane/ {IPL} of the projection is the plane {Z==0} of the NCS, and
  points of three-space are projected onto it from some /viewpoint/
  {O} with {O.Z > 0}.
  
  ANCHORED SPHERES
  
  One cannot determine the actual radius of a sphere {S} from its
  appearance on a perspective image. A sphere with twice the size of
  {S} that is placed twice as far from the viewpoint {O}, in the same
  direction, will have exactly the same projection on the image plane.
  
  Therefore, we may as well assume that all spheres are /anchored/,
  that is, have their center on the image plane. Said another way, we
  represent the set of all spheres that have a given projection by the
  (unique) anchored sphere in that set.
  
  An anchored sphere {S} is defined by three parameters: its radius
  {S.R} and the X and Y coordinates of its center {S.K}. If a general
  sphere {S} (not necessarily anchored) has radius {S.R} and center
  {S.K} anywhere in 3-space, then the equivalent anchored sphere {S'}
  has center {S'.K = m*(S.K.X, S.K.Y)} and radius {S'.R = m*R}, where
  {m = O.Z/(O.Z - S.K.Z)}. */
  
typedef struct pst_sphere_t
  { r2_t K;     /* Nominal center on the image plane. */
    double R;   /* Nominal radius (in NCS units). */
  } pst_sphere_t;


/* COMPUTING A SPHERE'S PROJECTION FROM A KNOWN VIEWPOINT
  
  At present, the procedures below consider only /proper/ anchored
  spheres, whose projection on the image plane is finite and
  non-empty. An anchored sphere {S} is proper only if {S.R < O.Z}. 
  
  In general, the perspective projection of a proper sphere {S}
  is an ellipse {E}.  We use here the representation of ellipses
  defined in {ellipse_crs_t}, consisting of the center {E.ctr},
  the transverse radius (minor semidiameter) {E.rad}, and 
  a stretch vector {E.str}.
  
  If {O} is {[0,0,0,1]}, the ellipse becomes a circle with center
  {E.ctr == S.K}, radius {E.rad == S.R}, and null stretch vector. If
  the camera is at finite distance ({O.m != 0}), but {S.K} lies on the
  optical axis ({S.K == Q}), the projection is still a circle with
  {E.ctr == S.K}, but with apparent radius {E.rad} somewhat larger
  than {S.R}.
  
  In all other cases the projection is a non-circular ellipse {E} with
  a nonzero stretch vector {E.str}. The major axis of the ellipse
  is parallel to {E.str} and goes through {S.K}, {E.ctr}, and the optical
  center {Q} of the image. If {O} is finite, {Q} is finite, the
  centers {E.ctr} and {S.K} are distinct, and the minor semidiameter
  {E.rad} is somewhat larger than {S.R}. If {O} is infinite, {E.ctr}
  coincides with {S.K}, and {E.rad} is equal to {S.R}. */

ellipse_crs_t pst_sphere_to_ellipse(pst_sphere_t *S, hr3_point_t *O);
  /* Returns the ellipse that is the projection of a sphere
    that has nominal center {K} and nominal radius {R}.
    The viewpoint {*O} must be completely determined 
    (no {NAN} components). */

/* DETRMINING A SPHERE FROM ITS PROJECTION */

pst_sphere_t pst_sphere_from_ellipse(hr3_point_t *O, ellipse_crs_t *E);
  /* Returns the sphere {S} that has projection {E} from 
    viewpoint {*O}.
    
    The problem is overdetermined in general, and there may be no
    sphere whose projection from {*O} is {*E}. In those cases, the
    procedure finds a sphere {S} whose projection is close to {*E} in
    some sense.  The client may want to compute an ellipse {E'}
    with{} and compare it with {*E}.
    
    The client may specify that some coordinates of the viewpoint {*O}
    are unknown, by setting them to {NAN}. The procedure will then replace
    any such {NAN} in {*O} by an apropriate value, deduced from {E}
    and the other coordinates of {O}. The procedure does not 
    directly change coordiantes of {Q} that are not {NAN};
    but it may rescale {*O} by an arbitrary positive factor.
    Since the factor may be arbitrarily large, some non-zero elements 
    of {Q} may become zero.
    
    Let {[O.m,O.x,O.y,O.z]} be the homogeneous coordinates of {O}. The
    weight {O.m} must be defined and non-negative, and {O.z} must be
    either {NAN} or positive. Let 'N' denote {NAN}, '*' denote a nonzero,
    non-{NAN] value, and '@' denote zero or '*'. Then the following
    are the acceptable cases for {*O}, with the relevant constraints on {E}:
    
      A. Cylindrical projection ({O.m == 0}):
    
      A2. {[0,*,*,*]}: means oblique cylindrical projection from a known
          direction. The ellipse {*E} should have the proper aspect 
          and orientation.

          {[0,0,*,*]}: same as above.
          {[0,*,0,*]}: same as above.
          {[0,0,0,*]}: same as above (orthogonal).
          {[0,0,0,N]}: same as above (orthogonal).

      A3. {[0,*,*,N]}: cylindrical projection from a specific azimuth
          (or its opposite) and unknown elevation, orthogonal in the
          limit. The procedure finds the elevation from {*E} (possibly
          90 degrees). The stretch vector {*E} should be null or aligned
          with that azimuth.

          {[0,*,0,N]}: same as above.  
          {[0,0,*,N]}: same as above.  
          {[0,N,0,*]}: same as above. 
          {[0,0,N,*]}: same as above. 
          {[0,0,N,N]}: same as above. 
          {[0,N,0,N]}: same as above. 

      A4. {[0,N,*,*]}: means cylindrical projection from a direction
          orthogonal to the vector {(0,-O.z,+O.y)}. The procedure
          finds that direction from {*E}. The ellipse {E} should have 
          proper aspect and orientation.

          {[0,*,N,*]}: analogous to the above for the vector
          {(-O.z,0,+O.x)}. 
     
      A5. {[0,N,N,*]}: means cylindrical projection from an unknown
          direction.  The procedure finds the direction from {*E};
          inn general there are two solutions; one is chosen arbitrarily.
          There are no contraints on {*E}.

          {[0,N,*,N]}: same as above.
          {[0,*,N,N]}: same as above.
          {[0,N,N,N]}: same as above.

      B. Conical projection ({O.m > 0}): 

      B1. {[1,@,@,*]}: conical projection from a finite known point.
          The ellipse {*E} should have the correct
          aspect and orientation.

      B2. {[1,N,@,*]}: the focal length {F} is finite and fixed, and the
          optical center {Q} lies on a specific line parallel to the {X}
          axis. The procedure finds {Q.X} from {E}. The ellipse {E} 
          should have the right aspect and orientation.

          {[1,@,N,*]}: similar to the above, for the {Y} axis. 

      B3. {[1,N,N,*]}: the focal length {F} is finite and fixed, but
          the optical center {Q} is unknown. The procedure finds {Q}
          from {E}; there are usually two solutions, one is chosen
          arbitrarily. There are no constraints on {*E}.

      B4. {[1,@,@,N]}: the optical center {Q} is finite and fixed, but
          the focal length {F} is unknown and may be infinite in the
          limit. The procedure finds {F} (possibly {+INF}) from {E}.
          The ellipse {*E} should have the proper  orientation.
          
      B5. {[1,N,N,N]}: the viewpoint is completely unknown and may be
          infinite in the limit. The viewing dstence cannot be
          determined, so the procedure assumes {O} is infinite: it
          sets {O} to {[0,N,N,N]} and proceeds as in case A5 above.

      B6. {[1,N,@,N]}: The viewpoint lies anywhere on a plane 
          parallel to the {X,Z} plane, and may be infinite in the
          limit. If the major axis of {E} meets the plane
          at a point {(A.X,A.Y)}, the procedure sets 
          {O.m=1}, {O.x=A.X}, {O.y=A.Y}, and falls into case B4.
          Otherwise the procedure sets {O.m=0}, {O.y=0}, and falls 
          into case A3 above.
          
          {[1,@,N,N]}: similar to the above, swapping {X} and {Y}.
          
     */

/* MAPPING BETWEEN 2D AND 3D COORDINATES */

hr3_pmap_t pst_sphere_perspective_view_map(pst_camera_t *C, ellipse_crs_t *E);
  /* Computes a homogeneous projective map {M} needed by
    {pst_normal_map_from_proc}. See {pst_sphere_basic.h} for
    details. */
    
    
/*
    
    The coordinate transformation map {M} computed by this procedure
    can be used to compute the ICS coordinates {(X,Y,Z)} of a point
    {P} on the visible part of the sphere, given only the image
    coordinates {(X',Y')} of its projection {P'}. From that one can
    easily compute the ICS coordinates of the unit vector normal to
    the sphere at that point. See
    {pst_sphere_perspective_compute_point_and_normal}. */
    

#endif
