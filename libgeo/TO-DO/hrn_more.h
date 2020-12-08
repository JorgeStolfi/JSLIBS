/* Last edited on 2007-10-14 14:52:54 by stolfi */
/* REGULAR SIMPLEX */

double hr2_pt_pt_diff(hr2_point_t *p, hr2_point_t *q);
  /* Distance between {p} and {q} in the spherical model; that is,
     angle between the vectors {p.c} and {q.c} in {R^3}, in radians. */

sign_t hr2_orient(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r);
  /* Returns the orientation (turning sense) of the triangle {p q r}.  
    Specifically, the triangle {[1 0 0] [0 1 0] [0 0 1]} has positive
    orientation. Returns 0 iff the points are collinear.
    May give inconsistent results for points very close to collinear. */

hr2_line_t hr2_join(hr2_point_t *p, hr2_point_t *q);
  /* Return the line through {p} and {q}. */

hr2_point_t hr2_meet(hr2_line_t *K, hr2_line_t *L);
  /* Return the point common to {K} and {L}. */

r2_t hr2_point_point_dir(hr2_point_t *frm, hr2_point_t *tto);
  /* Direction (a unit-length vector) of point {tto} seen from point {frm}.
    Works even if one of them is at infinity.  Does not work if both
    are at infinity, or coincident, or antipodal. */
    
r2_t hr2_line_dir(hr2_line_t *L);
  /* The direction of line {L}: a unit vector which, on the hither side
    of the plane, is parallel to {L} and travels around its positive
    side in the counterclockwise sense. Assumes {L} is not at
    infinity. */
    
r2_t hr2_line_normal(hr2_line_t *L);
  /* The normal direction of line {L}: a unit vector which,
    on the hither side of the plane, points from {L} into 
    {L}'s positive halfplane.  Assumes {L} is not at infinity. */
    
bool_t hrn_pmap_is_identity(hrn_pmap_t *M);
  /* TRUE iff the domain and counter-domain of map {M} are
    the same space, and {M} is the identity map of that space
    (apart from homogeneous scaling). */


hrn_pmap_t hrn_pmap_from_points(int m, int n, double p[], double u[]);
  /* Returns a projective function that takes the {m+1} cardinal
    points of {T^m} (the unit vectors of {R^{m+1}}) to the {m+1}
    points of {T^n} contained in {p}; and also some point {s} of {T^m},
    of the form {[±1,±1,...,±1]}, to point {u} of {T^n}.
    
    The array {p} is interpreted as a matrix with {m+1} rows
    and {n+1} columns
    
    The point {s} is unique, and is the signature of {u} relative to
    the ordered point tuple {p}. 
    
    The procedure fails if ???. */
