#ifndef pst_sphere_basic_H
#define pst_sphere_basic_H

/* pst_sphere_basic.h -- basic geometric tools for images of spheres. */
/* Last edited on 2024-11-04 07:26:30 by stolfi */

#include <ellipse_crs.h>
#include <r2.h>
#include <r3.h>
#include <hr3.h> 
#include <hr3_pmap.h> 
#include <argparser.h>

#include <pst_basic.h>
#include <pst_camera.h>

/* See {pst_sphere.h} for definition of sphere parameters. */
/* See {ellipse_crs.h} for definition of ellipse parameters. */
/* Canonical parameter order: {G,D,R,dst,rad,maj}. */

void pst_sphere_G_D_R_to_dst_rad_maj
  ( double G,      /* Angular spread {1/F} of camera (possibly zero). */
    double D,      /* Dist from optical center to sphere center (finite). */
    double R,      /* The radius of the sphere. */
    double *dst,   /* (OUT) Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
    double *rad,   /* (OUT) Length of minor semidiameter of ellipse. */
    double *maj    /* (OUT) Length of major semidiameter of ellipse. */
  );
  /* Computes the elements of the image {E} of an anchored sphere
    {S} under a projection with finite optical center {Q}.
    
    The procedure is given the angular spread {G} of the camera, the
    sphere radius {R}, and the distance {D} from the optical center of
    the image to the sphere's center {K}. The distance {D} should be
    finite, meaning that the projection is either conical, from a
    finite viewpoint, or orthogonal cylindrical.
    
    The procedure computes the minor and major semidiameters
    {*rad,*maj} of the projection of a sphere, and the distance {dst}
    from {Q} to the center of the ellipse {E} (which is usually
    different from {K}). 
    
    If {G} is zero, the parameter {dst} is irrelevant, and {*D} will
    be set to {dst}. The parameters {rad,maj,dst} may be NULL. */

void pst_sphere_G_dst_rad_to_D_R_maj
  ( double G,      /* Angular spread {1/F} of camera. */
    double dst,    /* Dist from {Q} to center {ctr} of ellipse (finite). */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double *D,     /* (OUT) Dist from optical center to sphere center (finite). */
    double *R,     /* (OUT) The radius of the sphere. */
    double *maj    /* (OUT) Length of major semidiameter. */
  );
  /* Computes the elements of a sphere {S} given those of its
    elliptical projection {E}, when the viewpoint is known.
  
    The procedure is given the camera's angular spread {G}, the minor
    semidiameter {rad} of {E}, and the distance {dst} from the optical
    center {Q} of the image to the center {ctr} of {E}.
    
    The procedure computes the radius {*R} of the sphere, and the
    distance {*D} from {Q} to the sphere's center. It also computes
    the major radius {*maj} of the ellipse, as deduced from
    {G,dst,rad}.
    
    The output {*D} may be infinite if {dst} is infinite or {G} is zero.
    The parameters {R,D,maj} may be NULL. */

void pst_sphere_dst_rad_maj_to_G_D_R
  ( double dst,    /* Dist from {Q} to center {ctr} of ellipse (finite). */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double maj,    /* Length of major semidiameter of ellipse. */
    double *G,     /* (OUT) Angular spread {1/F} of camera. */
    double *D,     /* (OUT) dist from {Q} to center {K} of sphere. */
    double *R      /* (OUT) The radius of the sphere. */
  ); 
  /* Computes the elements of a sphere {S} from those of its
    elliptical projection {E}, when the optical center {Q} of the
    image is known but the focal length {F} is not.
    
    The procedure is given the minor and major semidiameters {rad,maj}
    of its elliptical projection {E}, and the distance {dst} from the
    optical center {Q} of the image to the center {ctr} of the
    projection.
    
    The procedure computes the angular spread {*G = 1/F} of the
    camera, the radius {*R} of the sphere, and the distance {*D} from
    {Q} to the sphere's center {K}.
    
    If {rad == maj}, then {*G} cannot be determined from the
    data: either {dst} is zero and {G} can be anything, or {G} is zero
    and {dst} is meaningless. In that case the procedure will
    arbitrarily set {*G} to 0 (meaning an orthogonal cylindrical
    projection) and {*D} to {dst}.
    
    If {rad != maj} and {dst} is infinite, then {*D} will be infinite, {*G}
    will be 0, and {*R} will be {rad} will be computed from the ratio
    {maj/rad}. The parameters {G,R,D} may be NULL. */

void pst_sphere_G_rad_maj_to_R_D_dst
  ( double G,      /* Angular spread {1/F} of camera. */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double maj,    /* Length of major semidiameter of ellipse. */
    double *R,     /* (OUT) The radius of the sphere. */
    double *D,     /* (OUT) dist from {Q} to center {K} of sphere, or {+INF}. */
    double *dst    /* (OUT) Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
  ); 
  /* Computes the elements of a sphere {S} from those of its
    elliptical projection {E}, when the the focal length {F} is
    known but optical center {Q} of the image is not.
    
    The procedure is given the minor and major semidiameters {rad,maj}
    of its elliptical projection {E}, and the angular spread {*G =
    1/F} of the camera.
    
    The procedure computes the radius {*R} of the sphere, and the
    distances {*D,*dst} from the optical center {Q} to the centers
    {K,ctr}, of {S} and {E}, respectively.
    
    If {rad == maj}, then {*D} and {*dst} are set to 0, and {*R} is
    computed from {rad} and {G}.
    
    If {rad != maj} and {G} is zero, then {*D} and {*dst}
    are set to {+INF} and {*R} is set to {rad}.
    The parameters {R,G,dst} may be NULL. */

/* PERSPECTIVE VIEW MATRICES */    
    
hr3_pmap_t pst_sphere_perspective_view_map(hr3_point_t *O, r2_t *K, double R);
  /* Computes a homogeneous projective map {M} needed by
    {pst_normal_map_from_proc}.
    
    The map {M} takes the homogeneous coordinates of a point in the 
    image coordinate system (ICS) and yields its homogeneous coordinates
    in another system, called the sphere's /view coordinate system/ (VCS).

    The VCS coordinates are {[t,u,v,w]} where {t} is the weight; so
    that the corresponding Cartesian coordinates are {(u/t,v/t,w/t)}.
    The VCS is defined in terms of the /horizon/, which is the circle
    on the sphere that corresponds to the silhouette of the
    projection. Note that the center {K'} of the horizon is located
    inside the sphere, somewhere between the sphere's center {K} and
    the point {K''} of the sphere's surface that is closest to the
    camera. Note also that {K'} and {K''} project onto {K}, which is
    not always the center {E.ctr} of the sphere's projection.
    
    Let {B'} be any point on the horizon such that {B'.Z} is maximum.
    Note that {B'} projects to a point on the silhouette which lies at
    maximum distance from the optical axis. Let {A'} be another point
    of the horizon, located 90 degrees clockwise from {B'} as seen
    from the camera. Note that the segment {A--K} is parallel to the
    image plane ({A'.Z == K'.Z}), and that the segments {K'--A'} and
    {K'--B'} are perpendicular in three-space and on the image. Note
    however that the projections of these segments are not the
    semidiameters of {E}.
    
    The homogeneous VCS coordinate system is such that:
    
      (1) the sphere in question is the unit-radius origin-centered
      sphere of the VCS, with implicit equation {- t^2 + u^2 + v^2 +
      w^2 = 0};
      
      (2) the center {K'} of the horizon is the origin {[1,0,0,0]} of
      the VCS;
      
      (3) the closest point {K''} has VCS coordinates {[1,0,0,1]};
      
      (4) the points {A'} and {B'} have VCS coordinates {[1,1,0,0]}
      and {[1,0,1,0]}, respectively.
      
    It follows from (1) and (2) that the optical axis in the ICS  
    corresponds to the {W} axis of the VCS; and that the viewpoint {O}
    has VCS coordinates {[0,0,0,1]}, i.e. the point at {W == +INF}
    on the {W} axis.  Moreover, the plane of the horizon 
    corresponds to the plane {W == 0} of the VCS.  
    
    Thus, the perspective projection performed by the camera,
    from the visible half-space onto the image plane {Z == 0}, 
    corresponds in VCS coordinates to a vertical projection of the 
    whole space onto some plane {IPL*} in VCS space.  The plane {IPL*}
    usually cuts obliquely through the unit sphere, and does not
    pass through the origin. */

#endif
