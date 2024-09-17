/* Last edited on 2024-09-04 20:40:57 by stolfi */

/* RATIONAL PROJECTIVE MAPS */

typedef struct hi2_pmap_t { i3x3_t dir; i3x3_t inv; } hi2_pmap_t;
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */

bool_t hi2_pmap_is_identity(hi2_pmap_t *m);
  /* TRUE iff {m} is the identity map (apart from homogeneous scaling). */

hi2_point_t hi2_map_point(hi2_point_t *p, hi2_pmap_t *m);
  /* Applies projective map {m} to point {p}. */

hi2_line_t hi2_map_line(hi2_line_t *L, hi2_pmap_t *m);
  /* Applies projective map {m} to line {L}. */

hi2_pmap_t hi2_comp_map(hi2_pmap_t *m, hi2_pmap_t *n);
  /* Returns the composition of {m} and {n}, applied in that order. */

hi2_pmap_t hi2_inv_map(hi2_pmap_t *m);
  /* Returns the inverse of map {m}. */

hi2_pmap_t hi2_pmap_from_points(hi2_point_t *p, hi2_point_t *q, hi2_point_t *r, hi2_point_t *u);
  /* Returns a projective map that takes the cardinal points {[1,0,0]},
    {[0,1,0]}, and {[0,0,1]} to {p}, {q}, and {r}, respectively; and
    also some point {s} of the form {[±1,±1,±1]} to {u}. 
    
    The point {s} is unique, and is the signature of {u} relative to
    the ordered triple {p,q,r}.
    
    The procedure fails if the set {p,q,r,u} contains three collinear
    points.  */
