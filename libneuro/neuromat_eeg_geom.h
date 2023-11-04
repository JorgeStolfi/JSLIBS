#ifndef neuromat_geom_H
#define neuromat_geom_H

/* NeuroMat geometry tools. */
/* Last edited on 2023-10-21 21:47:19 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <r3.h>

#include <neuromat_eeg.h>

/*  SCHEMATIC SCALP COORDINATES

    In the absence of data on actual shape of the subject's scalp and
    on the actual electrode posiitons, we generally work with
    a /schematic scalp/ and /schematic positions/.
    
    The /schematic 3D scalp/ is the unit sphere of {\RR^3}, with the
    zenith point of the scalp represented by the north pole {(0,0,1)}
    of the sphere, the and the nose towards the {(0,1,0)} vector. The
    /3D schematic position/ of a point on the scalp is therefore a
    unit vector of {\RR^3}.
    
    The /2D schematic scalp/ is the plane {\RR^2}. The /2D schematic
    position/ of a point (on the scalp or off it) is its 3D schematic position,
    mapped to the plane by stereographic projection
    from the south pole {(0,0,-1)} onto the equatorial plane {Z=0}. Thus the upper
    hemisphere of the schematic 3D scalp maps to the unit disk of the
    schematic 2D scalp, while the lower hemisphere maps to the
    complement of that disk.  
    
    IDEALIZED SCALP COORDINATES
    
    The /idealized 3D scalp/ is an ellipsoid in {\RR^3} that results
    from stretching the schematic 3D scalp (unit sphere) by factors
    {(rX,rY,rZ)} along the three coordinate axes. These three numbers
    are the /idealized dimensions/ of the scalp, namely the half-width
    from ear to ear, the half-length from nose to occiput, and the
    half-height from ears to top. The /idealized 3D positions/ of a
    point on the scalp are its schematic 3D coordinates stretched by
    those factors.
    
    The /idealized 2D scalp/ is again the plane {\RR^2}. The
    /idealized 2D position/ of apoint on the scalp is the
    stereographic projection of its idealized 3D position from the
    point {(0,0,-rZ)} onto the plane {Z=0}. Thus the upper half of the
    scalp (with {Z>0}) projects to the interior of the ellipse on the
    idealized plane with center {(0,0)} and radii {(rX,rY)}. */

void neuromat_eeg_geom_get_schematic_2D_points(char *capType, int32_t *ne_P, char ***chname_P, r2_t **pos2D_P);
  /*  Returns in {*ne_P} the number {ne} of electrodes in the cap of type {capType}.
    Reaturns in {*chname_P} a vector {chname[0..ne-1]} with the names
    of those electrodes, and in {*pos2D_P} a vector {pos2D[0..ne-1]} with the
    two-dimensional schematic positions of the electrodes.
    
    Currently supports only {capType} equal to "R20", "R128", "R129", and "FN3". 
    The "R129" cap is the same as "R128" plus the voltage reference ("CZ") electrode as {pos2D[128]}. */
    
r2_t *neuromat_eeg_geom_get_schematic_2D_points_by_name(char *capType, int32_t ne, char *chname[]);
  /* Returns a vector {pos2D[0..ne-1]} with the schematic 2D coordinates
    of the electrodes whose names are {chname[0..ne-1]}.
    
    Assumes that those electrodes are a subset of the full set of
    electrodes availabe on the cap {capType}, whose names are
    provided by {neuromat_eeg_get_channel_names(capType,0,NULL)}. Thus
    {capType} must be one of the cap types recognized by that
    function. */

/* STEREOGRAPHIC PROJECTION 2D TO/FROM 3D  */
    
r3_t neuromat_eeg_geom_3D_from_2D(r2_t *p, r2_t *rad2, r3_t *rad3);
r2_t neuromat_eeg_geom_2D_from_3D(r3_t *p, r3_t *rad3, r2_t *rad2);
  /* If {rad2} and {rad3} are {NULL}, these procedures map between schematic 2D and 3D 
    scalp coordinates.
    
    Note that the mapping from {3D} to {2D} is undefined for point
    with {Z = -1}. For points that are not on the scalp, the 2D
    position is NOT the same as projecting radially onto the unit
    sphere and then applying the stereographic projection.
    
    When going from {2D} to {3D} the result is always a unit vector of
    {R^3} distinct from the south pole {(0,0,-1)} (apart from rounding
    errors).
    
    If {rad2} is not {NULL}, it is assumed to be a vector with the
    half-radii {(rX,rY)} of the idealized 2D scalp, and the 2D
    coordinates (input or output) are interpreted as coordinates on
    that idealized 2D scalp. Similarly if {rad3} is not {NULL}, it is
    assumed to be a vector with the half-radii {(rX,rY,rZ)} of the
    idealized 3D skull, and the 3D coordinates (input or output) are
    interpreted as coordinates on that idealized scalp.
    */
   
/* MAPPING UNIT DISK TO/FROM ELLIPSE */

r2_t neuromat_eeg_geom_disk_from_ellipse(r2_t *p, r2_t *ctr, r2_t *rad);
r2_t neuromat_eeg_geom_ellipse_from_disk(r2_t *p, r2_t *ctr, r2_t *rad);
  /* These procedures perform translations and non-uniform scalings
    of schematic {2D} coordinates, so that the unit disk is 
    mapped to the ellipse with center {ctr} and major radii {rad.c[0],
    rad.c[1]}.  */

void neuromat_eeg_geom_map_many_disk_to_ellipse(int32_t np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[]); 
void neuromat_eeg_geom_map_many_ellipse_to_disk(int32_t np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[]); 
  /* These proceedures apply {neuromat_eeg_geom_ellipse_from_disk}
    and {neuromat_eeg_geom_disk_from_ellipse}, respectively, to the points 
    {p[0..np-1]} yielding the points {q[0..np-1]}. */

#endif
