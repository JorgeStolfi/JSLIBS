/* r2_extra.h --- additional operations on points and vectors of {R^2} */
/* Last edited on 2024-12-05 10:27:49 by stolfi */

#ifndef r2_extra_H
#define r2_extra_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <interval.h>
#include <r2.h>
#include <r3x3.h>
#include <r2x2.h>

typedef void r2_map_jacobian_t (r2_t *p, r2x2_t *J);
  /* Type of a procedure that defines a geometric transformation
    from {R^2} to {R^2}.
    
    The procedure should apply the geometric map to point {*p} and
    store the result back into {*p}. If {J} is not NULL, the procedure must
    also multiply {*J} by the Jacobian of the geometric map, evaluated at {*p}.
    In the Jacobian, element {J[i][j]} is the derivative of output {p.c[j]}
    with respect to input {p.c[i]}. 
    
    If the map is not defined at the given point {*p}, the procedure
    should set {*p} to {(NAN,NAN)}; in that case, {*J} need not be
    modified. The procedure should do this, in particular, if either
    coordinate of input {*p} is {NAN}. */

void r2_map_projective(r2_t *p, r3x3_t *M, r2_t *q, r2x2_t *J);
  /* Applies to {p} a projective transformation defined by the matrix {M},
    writing the result to {*q}.  Also post-multiplies the matrix {J}
    by the Jacobian of the map at {p}. */

void r2_map_radial(r2_t *p, r2_t *h, double kappa, r2x2_t *J);
  /* Applies to point {p} a radial deformation map with parameter
    {kappa}, centered at the origin. The inverse map can be obtained
    by negating {kappa}. Also post-multiplies the matrix {J} by the
    Jacobian of the map.
    
    The map is meant to be equivalent to that of R. Tsai's 1985 paper
    for small {kappa}, but it is not quite the same.
    
    The unit of the parameter {kappa} is {1/mm^2}. The procedure
    assumes that {h.c[0]} and {h.c[1]} are the width and height of
    each pixel, respectively, in mm. The pincushion distortion is
    undefined (returns {(NAN,NAN)}) if {p} lies on or outside the
    radius {R = sqrt(1/(2*kappa))} also in millimeters. So {kappa}
    should not exceed {2/d^2} where {d} is the diagonal of the image
    in millimeters. */

void r2_map_twirl(r2_t *p,  r2_t *ctr, double rad, double ang, r2x2_t *J);
  /* Applies to {p} a `twirl' or `whirlpool' map centered at
    {ctr}. Also post-multiplies the matrix {J} by the
    Jacobian of the map.
    
    Every point is rotated around {ctr} by some angle. Points near {ctr}
    are rotated by {ang} radians. Points at distance {rad} are rotated
    by half that amount. As the distance {r} from {ctr} increases, the
    rotation angle decreases like {1/r^2}. */

void r2_map_expand(r2_t *p, double xlo, double xhi, double ylo, double yhi, r2x2_t *J);
  /* Maps the rectangle {[xlo_xhi]×[ylo_yhi]} to the plane {R^2} by
    expanding each coordinate separately with {expand_range}. Also
    post-multiplies the matrix {J} by the Jacobian of the map. */

void r2_map_contract(r2_t *p, double xlo, double xhi, double ylo, double yhi, r2x2_t *J);
  /* Maps the plane {R^2} to the rectangle {[xlo_xhi]×[ylo_yhi]} by
    contracting each coordinate separately with {contract_range}. Also
    post-multiplies the matrix {J} by the Jacobian of the map. */

void r2_map_compute_numeric_jacobian(r2_t *p, r2_map_jacobian_t *map, double step, r2x2_t *K, bool_t debug);
  /* Computes the jacobian {8K} of {map} by numeric differentiation with 
    the given{step}. */

void r2_map_check_jacobian(r2_t *p, r2_map_jacobian_t *map, char *mapname, double eps, bool_t debug);
  /* Checks the Jacobian {J} returned by {map} at {p} against the numerically estimated
    derivatives of {q} relative to {p}, where {q} is computed from {p}
    by the {map} procedure. The estimates are obtained by central
    divided difference with step not exceeding {eps} in the output or
    input domain. The {mapname} */

void r2_get_persp_rectangle_bbox
  ( interval_t tbox[],   /* Rectangle in true coordinates. */
    r3x3_t *T2I,         /* True-to-image projective map matrix. */
    interval_t ibox[]    /* (OUT) bonding box in image coordinates. */
  );
  /* Computes an enclosing box for the image of a rectangle. 
    Requires the rectangle {tbox[0] × tbox[1]}, in
    true coordinates, and the true-to-image matrix {T2I}.
    The box is returned as two intervals {ibox[0..1]},
    one for ech axis. */
  
void r2_get_persp_disk_bbox
  ( r2_t *ctr,         /* Disk center in true coordinates. */
    double rad,        /* Disk radius in true coordinates. */
    r3x3_t *T2I,       /* True-to-image projective map matrix. */
    interval_t ibox[]  /* (OUT) bonding box. */
  );
  /* Computes an enclosing box for the image of a disk. 
    Requires the disk's center {ctr} and radius {rad}, in
    true coordinates, and the true-to-image matrix {T2I}.
    The box is returned as two intervals {ibox[0..1]},
    one for ech axis. */

bool_t r2_pixel_is_inside_persp_rectangle
  ( int32_t x,              /* Pixel column in image. */
    int32_t y,              /* Pixel row in image. */
    double mrg,         /* Safety margin (pixels). */
    r3x3_t *I2T,        /* Image-to-true projective matrix. */
    interval_t tbox[]   /* Rectangle in true coordinates. */
  );
  /* Returns TRUE iff the image pixel in column {x} and row {y},
    lies well inside the image of a rectangle.  The pixel
    is modeled as a unit square with center at {(x+0.5,y+0.5)}
    fattened by {mrg} all around.

    Requires the rectangle {tbox[0] × tbox[1]} in true
    coordinates, and the {3×3} image-to-true homogeneous projective
    matrix {I2T}. */

bool_t r2_pixel_is_inside_persp_disk
  ( int32_t x,        /* Pixel column in image. */
    int32_t y,        /* Pixel row in image. */
    double mrg,   /* Safety margin (pixels). */
    r3x3_t *I2T,  /* Image-to-true projective matrix. */
    r2_t *ctr,    /* Center of disk in true coordinates. */
    double rad    /* Radius of disk in true coordinates. */
  );
  /* Returns TRUE iff the image pixel in column {x} and row {y},
    lies well inside the image of a disk.  The pixel
    is modeled as a unit square with center at {(x+0.5,y+0.5)}
    fattened by {mrg} all around.

    Requires the center {ctr} and radius {rad} of the disk in true
    coordinates, and the {3×3} image-to-true homogeneous projective
    matrix {I2T}. */

void r2_clip_seg_to_unit_disk(r2_t *a, r2_t *b, double *ta, double *tb);
  /* Finds the part {a1--b1} of segment {a--b} that lies inside the
    disk of unit radius and center at the origin.  Returns the 
    fractional positions {*ta,*tb} (in the range {[0 _ 1]}) of {a1} and {b1} 
    in the segment {a--b}.  
    
    In particular, {*ta} will be 0.0 iff {a} is inside the disk; {*tb}
    will be 1.0 ifb {b} is inside the disk; and {0 < *tb < *ta < 1} iff
    the entire segment is outside. */

/* DEBUGGING AIDS */

void r2_debug_point_jacobian(char *label, r2_t *p, r2x2_t *J, char *tail);
  /* Prints the point {(x/w,y/w)} to {stderr}, tagged with the
    string {label} and terminated by the string {tail}. If the
    Jacobian {J} is not NULL, prints it too. */

#endif
