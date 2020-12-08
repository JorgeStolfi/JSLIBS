#ifndef neuromat_geom_H
#define neuromat_geom_H

/* NeuroMat geometry tools. */
/* Last edited on 2013-11-29 01:11:07 by stolfilocal */

#define _GNU_SOURCE

#include <r2.h>
#include <r3.h>

#include <neuromat_eeg.h>

r2_t *neuromat_eeg_geom_get_schematic_2D_points(int ne);
  /* Returns a vector of {ne} points which are the two-dimensional
    schematic positions of the electrodes in the {ne}-electrode
    experiments. Assumes that the head is represented in the schematic
    diagram by the unit disk, as seen from above, with the nose at the top,. */

r2_t *neuromat_eeg_geom_get_20_schematic_2D_points(void);
  /* Same as {neuromat_eeg_geom_get_schematic_2D_points}
    specific to 20-electrode experiments. */

r2_t *neuromat_eeg_geom_get_128_schematic_2D_points(void);
  /* Same as {neuromat_eeg_geom_get_schematic_2D_points}
    specific to 128-electrode experiments. */
    
/* STEREOGRAPHIC PROJECTION 2D TO/FROM 3D  */
    
r3_t neuromat_eeg_geom_3D_from_2D(r2_t *p);
r2_t neuromat_eeg_geom_2D_from_3D(r3_t *p);
  /* These procedures map between schematic 2D and 3D scalp coordinates.
    The schematic 3D position is a point on the unit sphere of {R^3}, the 
    schematic 2D position is a point of the {XY} plane. 
    
    Assumes stereographic projection from the south pole {(0,0,-1)}. 
    Thus the origin is mapped to/from the north pole {(0,0,+1)}, the 
    {+Z} hemisphere is mapped to/from the unit disk of the {XY} plane, 
    and the edge of that disk (the equator of the unit sphere) is mapped to itself,
    i.e. merely gets an extra coordinate {Z=0}. 
    
    When going from {2D} to {3D} the result is always a  unit vector of {R^3}
    (apart from rounding errors). The mapping from {3D} to {2D} is undefined
    at the south pole. */
   
/* MAPPING UNIT DISK TO/FROM ELLIPSE */

r2_t neuromat_eeg_geom_disk_from_ellipse(r2_t *p, r2_t *ctr, r2_t *rad);
r2_t neuromat_eeg_geom_ellipse_from_disk(r2_t *p, r2_t *ctr, r2_t *rad);
  /* These procedures perform translations and non-uniform scalings
    of schematic {2D} coordinates, so that the unit disk is 
    mapped to the ellipse with center {ctr} and major radii {rad.c[0],
    rad.c[1]}.  */

void neuromat_eeg_geom_map_many_disk_to_ellipse(int np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[]); 
void neuromat_eeg_geom_map_many_ellipse_to_disk(int np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[]); 
  /* These proceedures apply {neuromat_eeg_geom_ellipse_from_disk}
    and {neuromat_eeg_geom_disk_from_ellipse}, respectively, to the points 
    {p[0..np-1]} yielding the points {q[0..np-1]}. */

#endif
