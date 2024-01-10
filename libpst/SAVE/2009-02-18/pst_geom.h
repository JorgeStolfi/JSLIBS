#ifndef pst_geom_H
#define pst_geom_H

/* pst_geom.h -- geometric tools for images of spheres. */
/* Last edited on 2009-02-28 05:57:32 by stolfi */

#include <ellipse_crs.h>
#include <r2.h>
#include <r3.h>
#include <hr3.h> 
#include <argparser.h>

#include <pst_basic.h>
#include <pst_camera.h>

/* TOOLS FOR IMAGES OF SPHERES
  
  The procedures in this interface deal with a perspective 
  (conical) /projection/ of spherical objects in 3D space 
  onto an /image plane/, specified by a camera description {C}. 
  
  IMAGE COORDINATE SYSTEM
  
  These procedures all use an orthogonal /image coordinate system/
  (NCS) associated to that plane. The {X} and {Y} axes of the NCS lie
  on the image plane, and the {Z} axis points out of the plane, into
  the half-space that contains the /viewpoint/, the center of the
  perspective projection. The image coordinates {(X,Y,Z)} are measured
  in the same unit.
  
  (In computer graphics one frequently uses an NCS whose origin is at
  the top left corner of the image, with the {X} axis pointing left
  and the {Y} axis pointing down, where each pixel is a unit square.
  Sometimes the origin is at the bottom left corner, or at the center
  of the image, and the Y axis may point up. These details are not
  relevant for the procedures in this interface. However, it is
  important that the three image coordinates have the same length
  unit, so that Euclidean distances can be computed by the Pythagorean
  formula.)
  
  See {pst_camera.h} for further details on the NCS and related
  concepts.
  
  PROPER AND ANCHORED SPHERES
  
  In this interface we consider only /proper spheres/, whose
  projection on the image plane is finite and non-empty. Futhermore
  we assume that the spheres are /anchored/, meaning that their
  center lies on the image plane. Note that for any proper sphere {S}
  there is an anchored proper sphere {S0} which has exactly the same
  projection as {S}.
  
  NOMINAL SPHERE PARAMETERS
  
  For the purposes of this interface, we define the /nominal center/
  {K} and /nominal radius/ {R} of a proper sphere as being the center
  and radius of the equivalent anchored sphere. The nominal center can
  be represented by two NCS coordinates {(K.X,K.Y)} since the {Z}
  coordinate is always zero.
  
  If a non-anchored
  sphere has center {K=(K.X,K.Y,K.Z)} and radius {R}, then its nominal
  center and radius are {m*(K.X,K.Y)} and {m*R}, where {m =
  F/(F-K.Z)}.  Note that the sphere is proper only if {K.Z+R < F}.
  
  SPHERE PROJECTION
  
  The perspective projection of a proper sphere is an ellipse.
  For the camera at infinity, the projection is a circle
  whose center and radius coincide with the nominal center
  and nominal radius of the sphere.  
  
  If the camera is not at infinity, but the nominal center is on the
  optical axis, the projection is a circle centered on the axis, with
  an apparent radius {rad} somewhat larger than the true radius.
  
  In the remaining cases the projection is an ellipse whose major axis
  goes trough the optical axis. The center {ctr} of that projection
  is somewhar displaced from the nominal center, along the major axis.
  The minor semidiameter {rad} is somewhat larger than the nominal radius.
  The ellipse has a nonzero stretch vector {str}, also pointing
  away from the optical axis. */

void pst_geom_sphere_compute_stretch_from_apparent_params
  ( pst_camera_t *C,  /* Camera projection. */
    r2_t *ctr,        /* The center of the sphere's projection on the image. */
    double rad,       /* The minor semidiameter of the projection. */
    r2_t *strP        /* (OUT) The major semidiameter of the projection. */
  );
  /* Computes the expected stretch vector {*strP} of a sphere's
    projection, given the center {ctr} of the projection and the
    transverse radius {rad} (minor semidiameter) of the projection.
    Note that {ctr} and {rad} are apparent parameters of the
    projection, not the nominal parameters of the sphere. */

/* GEOMETRY OF SPHERE PROJECTION FROM NOMINAL PARAMETERS  */

void pst_geom_sphere_compute_projection_from_nominal_params
  ( pst_camera_t *C,   /* Camera projection. */
    r2_t *tct,         /* The nominal center of the sphere. */
    double trd,        /* The nominal radius of the sphere. */
    ellipse_crs_t *EP  /* (OUT) The sphere's projection on the image. */
  );
  /* Returns the ellipse that is the projection of a sphere
    that has nominal center {tct} and nominal radius {trd}. */

void pst_geom_sphere_compute_apparent_params_from_nominal_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *tct,       /* The nominal center of the sphere. */
    double trd,      /* The nominal radius of the sphere. */
    r2_t *ctrP,      /* (OUT) The center of the sphere's projection on the image. */
    double *radP,    /* (OUT) The minor semidiameter of that projection. */
    r2_t *strP       /* (OUT) The stretch vector of the projection. */
  );
  /* Returns the center {*ctrP}, minor semidiameter {*radP}, and
    stretch vector {*strP} of the projection of a sphere that has nominal
    center {tct} and nominal radius {trd}.  If any of {ctrP,radP,strP}
    is {NULL}, the corresponding parameter is not computed. */

void pst_geom_sphere_compute_stretch_from_nominal_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *tct,       /* The nominal center of the sphere on the image. */
    double trd,      /* The nominal radius of the sphere on the image. */
    r2_t *strP       /* (OUT) the stretch vector of the sphere's projection. */
  );
  /* Computes the stretch vector {*strP} of the projection
    of the sphere with nominal center is {tct} and nominal radius {trd}. */

void pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *ctr,       /* The center of the sphere's projection on the image. */
    double rad,      /* The minor semidiameter of that projection. */
    r2_t *tctP,      /* (OUT) The nominal center of the sphere. */
    double *trdP,    /* (OUT) The nominal radius of the sphere. */
    r2_t *strP       /* (OUT) The correct stetch vector of the projection. */
  );
  /* Stores into {*tctP} and {*trdP} the nominal center and nominal radius
    of a sphere whose projected outline has center {ctr} and minor
    semidiameter {rad}. */

/* MAPPING BETWEEN 2D AND 3D COORDINATES */

hr3_pmap_t pst_geom_sphere_perspective_view_map(pst_camera_t *C, ellipse_crs_t *E);
  /* Computes a homogeneous projective map {M} needed by
    {pst_normal_map_from_proc}.
    
    The map {M} takes the homogeneous coordinates of a point in the 
    coordinate system (NCS) and yields its homogeneous coordinates
    in another system, called the sphere's /view coordinate system/ (VCS).

    The VCS coordinates are {[t,u,v,w]} where {t} is the weight; so
    that the corresponding Cartesian coordinates are {(u/t,v/t,w/t)}.
    The VCS is defined in terms of the /horizon/, which is the circle
    on the sphere that corresponds to the silhouette of the
    projection. Note that the center {K'} of the horizon is located
    inside the sphere, somewhere between the sphere's center {K} and
    the point {K''} of the sphere's surface that is closer to the
    camera.  Note also that {K}, {K'} and {K''} project onto the
    same point, which is not always the center {E.ctr} of the 
    sphere's projection.
    
    Let {B'} be the on the horizon such that {B'.Z} is maximum; or, if
    the whole horizon is parallel to the image, any point of the
    horizon). Note that {B'} projects to the point of the silhouette
    which is fathest from the optical axis. Let {A'} be another point
    of the horizon, located 90 degrees clockwise from {B'} as seen
    from the camera. Note that {A'.Z == K'.Z}, and that the segments
    {K'--A'} and {K'--B'} are perpendicular in three-space and on the
    image. Note however that the projections of these segments are not
    the semidiameters of {E}.
    
    The homogeneous VCS coordinate system is such that:
    
      (1) the sphere in question is the unit-radius origin-centered
      sphere of the VCS, with implicit equation {- t^2 + u^2 + v^2 +
      w^2 = 0};
      
      (2) the center {K'} of the horizon is the origin {[1,0,0,0]} of
      the VCS;
      
      (3) the closest point {K''} has VCS coordinates {[1,0,0,1]};
      
      (4) the points {A'} and {B'} have VCS coordinates {[1,1,0,0]}
      and {[1,0,1,0]}, respectively.
      
    It follows from (1) and (2) that the optical axis in the NCS  
    corresponds to the {W} axis of the VCS; and that the viewpoint {vpt}
    has VCS coordinates {[0,0,0,1]}, i.e. the point at {W == +INF}
    on the {W} axis.  Moreover, the plane of the horizon 
    corresponds to the plane {W == 0} of the VCS.  
    
    Thus, the perspective projection performed by the camera,
    from the visible half-space onto the image plane {Z == 0}, 
    corresponds in UCS coordinates to a vertical projection of the 
    whole space onto some plane {IPL*} in VCS space.  The plane {IPL*}
    usually cuts obliquely through the unit sphere, and does not
    pass through the origin.
    
    The coordinate transformation map {M} computed by this procedure
    can be used to compute the NCS coordinates {(X,Y,Z)} of a point
    {P} on the visible part of the sphere, given only the image
    coordinates {(X',Y')} of its projection {P'}. From that one can
    easily compute the NCS coordinates of the unit vector normal to
    the sphere at that point. See
    {pst_geom_sphere_perspective_compute_point_and_normal}. */

void pst_geom_sphere_perspective_compute_point_and_normal
  ( hr3_pmap_t *M, /* NCS-to-VCS coordinate mapping. */
    r2_t *pim,     /* Image coords of the projection {P'} of a sphere point {P}. */
    r2_t *tct,     /* Image coords of nominal center {K} of the sphere. */
    r3_t *psp,     /* (OUT) NCS coordinates of the point {P} on the sphere. */
    r3_t *nrm      /* (OUT) NCS coordinates of the sphere's normal at {P}. */
  );
  /* Computes the (3D) Cartesian NCS coordinates {psp} of a point {P}
    on the visible part of a sphere, given the (2D) image coordinates
    {pim} of its projection {P'}.  Requires the NCS-to-VCS homogeneous
    coordinate transformation matrix computed by
    {pst_geom_sphere_perspective_view_map}.
    
    This computation succeeds only if the given point {P' == pim} lies
    inside the sphere's  projection.  If {P'} falls outside the sphere's
    projection, the computed point {P == *psp} and the normal {N == *nrm}
    are set to NANs. */
    
/* UTILITIES */

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double ard);
  /* Clips the direction vector {udir} to a spherical cap with center
    on the direction vector {sdir} and angular radius {ard. Namely, if
    the angle between {udir} and {sdir} is greater than {ard} radians,
    sets {udir} to the unit vector that lies on the shortest arc
    between the two, at {ard} radians from {sdir}. */

ellipse_crs_t *pst_geom_sphere_new(void);
  /* Alocates a new {elipse_t} record, not initialized. */

#endif
