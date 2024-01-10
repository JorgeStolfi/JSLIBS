/* General tensor-style Bézier patches. */
/* Last edited on 2009-08-23 13:24:47 by stolfi */

/* ??? TO DO: use <indexing.h>! */

#ifndef bz_patch_H
#define bz_patch_H

#include <psp_basic.h>
#include <bz_basic.h>

#include <interval.h>
#include <box.h>
#include <stdint.h>

/* 
  BÉZIER REPRESENTATION FOR MULTIVARIATE POLYNOMIALS
  
  A polynomial which has maximum degree {g[i]} on coordinate {i}, for
  each {i} in {0..d-1}, can be specified by its /Bézier coefficients/
  relative to some box {B} of {R^d}. These are an array of real
  values, conceptually located at a regular lattice of points spanning
  {B}, with {g[i]+1} points along each spanning axis {i}, and one
  point along each stabbing axis.
  
  BERNSTEIN-BÉZIER TENSOR ELEMENTS
  
  Each Bézier coefficient is specified by a tuple {e[0..d-1]} of {d}
  /Bézier indices/, where each index {e[i]} ranges over {0,.. g[i]}. The
  coefficient with index tuple {e[0..d-1]} multiplies the
  /Bernstein-Bézier tensor element/
  
    { PROD { BB^g_{e[i]}(z[i]) : i in {0,..d-1} } }
    
  where {z[0..d-1]} are the coordinates of the point relative
  to the tile {B}; and {BB^g_e(z)} is the /Bernstein-Bézier
  polynomial/ of degree {g} and index {e}. */

typedef psp_ddim_t bz_patch_ddim_t;
  /* Dimension of the domain of a Bézier patch. */

typedef psp_rdim_t bz_patch_rdim_t;
  /* Dimension of the range of a Bézier patch. */

#define bz_patch_MAX_DDIM (4)
  /* Maximum dimension for the domain of a Bézier patch. Note that
    some vectors are preallocated with this size, and that the number
    of coeffs is exponential in the domain dimension. */

#define bz_patch_MAX_RDIM (10)
  /* Maximum range dimension (number of coordinates in each coeff) for
    a Bézier patch, just to be safe. Could be increased with little
    harm. */

#define bz_patch_MAX_DEGREE (6)
  /* Maximum degre of a Bézier patch along each coordinate,
    just to be safe. Could be increased with modest harm. */

typedef bz_degree_t bz_patch_degrees_t[bz_patch_MAX_DDIM];
  /* The degree of a Bézier patch, along each domain axis. */

typedef bz_degree_t bz_patch_cpindex_t;
  /* Each coeff has a {bz_patch_cpindex_t} along each domain axis,
    which ranges from 0 to the corresponding degree, inclusive. */

typedef struct bz_patch_t /* A tensor-style Bezier patch. */
  { bz_patch_ddim_t m;      /* Dimension of domain (argument space). */
    bz_patch_rdim_t n;      /* Dimension of range (image space). */
    bz_patch_degrees_t g;   /* The degree along axis {i} is {g[i]}. */
    double *c;              /* Bezier coeffs for patch. */
  } bz_patch_t;
  /* A {bz_patch_t} record {b} describes a polynomial map from {R^m} to
    {R^n}.
    
    The map {b} is defined as as a linear combination of Bézier
    elements, which are products of {m} univariate Bernstein-Bézier
    polynomials. Factor {i} of this product is a polynomial of degree
    {g[i]} on coordinate {i} of the argument. The coefficients of this
    linear combination are the /Bézier coefficients/ of the patch,
    /coeffs/ for short. Each coefficient is a vector of {R^n}.
    
    The coeffs are logically organized as an array with {m} subscripts
    {e[0],..e[m-1]}. Each subscript {e[i]} ranges over {0..g[i]},
    therefore the total number {nCoefs} of coeffs is {PROD{(g[i]+1) :
    i = 0..m-1}.
    
    The coeffs are packed in a one-dimensional array {b.c} with {nCoefs} elements,
    each of them consisting of {n} consecutive {double}s. 
    The coeff with indices {e[0],..e[m-1]} is stored
    in {b.c[bix(e[0],..e[m-1])]}, where 
    
      { bix(e[0],..e[m-1]) = SUM{ bstep[i]*e[i] : i = 0..m-1 } }
      
    and
    
      { bstep[i] = PROD{g[j]+1 : j = i+1..m-1} }
    
    Therefore, the value of {b(x)} on a point {x} of {R^m} is the sum
    of
    
      { b.c[bix(e[0],..e[m-1])] * PROD{ B_{e[i]}^{g[i]}(x[i]) : i = 0..m-1 } }
    
    for all coeffs {b.c[e[0],..e[m-1]]}, where {B_k^r(z) ==
    z^k*(1-z)^{r-k}*\choose(r,k)/2^r} is the Bernstein-Bezier
    polynomial of index {k} and degree {r}. Of course, the value of
    {b(x)} can be computed more efectively by successive
    interpolations along each coordinate (the DeCasteljau algorithm).
    
    Although the map is defined over all {R^m}, for most applications
    it is convenient to define the patch proper as the restriction of
    {b} to the unit {m}-cube {U = [0 _ 1]^m} of {R^m}. In particular,
    a Bezier coefficient {b.c[e[0],..e[m-1]]} where every subscript
    {e[i]} is ether {0} or {g[i]} is a `corner' of the patch (the
    image of a corner of {U}). The remaining coeffs (which
    exist only if some {g[i]} is greater than 1) affect the shape of
    the patch between those corners.
    
    The array of coeffs is stored in the one-dimensional
    
    */

bz_patch_t bz_patch_new(bz_patch_ddim_t m, bz_patch_rdim_t n, const bz_degree_t g[]);
  /* Returns a Bézier patch {b} of the specified dimensions and 
    degrees {g[0..m-1]}. The coeeficients {b.c} are newly allocated. */

bz_patch_t bz_patch_uniform_new(bz_patch_ddim_t m, bz_patch_rdim_t n, bz_degree_t g);
  /* Returns a Bézier patch {b} of the specified dimensions and 
    uniform degree vector {(g,..g)}. The coeeficients {b.c} 
    are newly allocated. */

bz_patch_t bz_patch_facet_new
  ( bz_patch_ddim_t m, 
    bz_patch_rdim_t n, 
    const bz_degree_t g[], 
    box_axis_index_t j
  );
  /* Returns a Bézier patch {f} suitable for storing into it a facet
    of another Bézier patch {b}, which has the specified dimensions
    and degrees {g[0..m-1]}. The facet will be perpendicular to axis
    {j} of {b}'s domain. Thus, the new patch {f} will have domain
    dimension {m-1}, range dimension {n}, and its degree vector will
    be {g[0..m-1]} with element {g[j]} deleted. The coefficients
    {f.c} are newly allocated. */

void bz_patch_free(bz_patch_t b);
  /* De-allocates the coefficient vector and other internal storage
    areas of {b} (but not {b} itself). */

double *bz_patch_control_point(bz_patch_t *b, bz_patch_cpindex_t e[]);
  /* Given a Bézier patch {b} of degree {g = b->g} from {R^m} to {R},
    and a vector {e[0],..e[m-1]} of indices, where each {e[i]} lies in
    {0..b->g[i]}, returns the address of the coeff
    {b.c[e[0],..e[m-1]]}. Note that the coeff is neither
    allocated nor copied --- the result is just a pointer into the
    control array {*(b->c)}.. */

void bz_patch_eval(bz_patch_t *b, double x[], double bx[]);
  /* Stores in {bx[0 .. b->n-1]} the image {b(x)} of the  
    point {x[0..m-1]} of {R^m} by the Bézier patch {b}. 
    The vector {bx} must be allocated by the caller, with 
    {b->n} elements. */

void bz_patch_eval_box(bz_patch_t *b, interval_t box[], bz_patch_t *f);
  /* Stores in {*f} a Bézier patch that describes the restriction
    of patch {b} to of the given axis-aligned {box[0..m-1]},
    reparametrized to the unit cube {U}. The patch
    {*f} must be allocated by the caller, with the same dimensions
    and degree as {b}. */

void bz_patch_compute_bbox(bz_patch_t *b, interval_t box[]);
  /* Stores in {box[0..d-1]} an axis-aligned bounding box for the 
    patch {b} --- that is, a box that is guaranteed to contain
    the image point {b(x)}, when the argument point {x} ranges over
    the unit cube {U} of {R^d}. */

void bz_patch_get_facet(bz_patch_t *b, box_axis_t ax, interval_side_t dir, bz_patch_t *t);
  /* Given a Bézier patch {b} with domain {D} of dimension {m}, stores
    in {t} the restriction of {b} to the {m-1}-dimensional face
    (facet) of {D} that is perpendicular to axis {ax}, and is 
    located in the {dir} direction of that axis with respect to {D}. 
    
    The patch {t} must be allocated by the caller, with domain
    dimension {b->m-1}, range dimension {b->n} and degree vector equal
    to {b->g} with element {ax} deleted. */

void bz_patch_get_face( 
    bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  );
  /* Stores in {t} a Bezier patch that describes the shape of an 
    {s}-dimensional face {F} of patch {b}.  The face is selected
    by its signature vector {dir[i]}, as explained above.
    
    The patch {t} must be allocated by the caller, with range dimension
    {t->n = b->n}. Its domain dimension {t->m} must be the 
    dimension of the face, i.e. the number of elements {dir[i]}
    which are equal to {SMD}; and its degree vector {t->g}
    must contain the elements of {b->g} that correspond to 
    those axes. */

void bz_patch_set_face(
    bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  );
  /* Assumes that {*t} is a Bezier patch that describes the shape of an 
    {s}-dimensional face {F} of patch {b}, specified by its signature {dir}. 
    Modifies {b} so that the face {F} has that shape.  The remaining control
    points of {b} are modified too, proportionally to their distance from 
    face {F}. */
    
double bz_patch_try_flatten_face(
    bz_patch_t *b, 
    box_signed_dir_t dir[], 
    double tol
  );
  /* Checks a face {f} of patch {b} described by {dir[]}; if it can be
    approximated by a single multiaffine (degree 1) patch, with error
    at most {tol} in every coordinate, makes it flat. The face is
    selected by its signature vector {dir[i]}, as explained above. In
    any case, returns the final maximum deviation of face {f} from
    flatness.

    The decision to flatten a face depends only on the Bézier coefficients
    of that face.  Thus, if two patches {b1,b2} share a face, they
    will continue to do so after being flattened by this procedure. */ 

void bz_patch_raise_degree(bz_patch_t *b, bz_patch_t *t);
  /* Stores in {t} a Bézier patch of degree vector {t->g} which coincides
    with the patch {b}.  The patch {t} must be allocated by the caller,
    with the same dimensions as {b}, and must satisfy {t->g[i] \geq b->g[i]}
    for all {i}. */

bz_degree_t bz_patch_max_degree(bz_patch_t *b);
  /* The maximum degree of {b} along any domain axis. */

void bz_patch_multiaffine_approx(bz_patch_t *b, bz_patch_t *t);
  /* Stores in {*t} a multiaffine approximation (i.e., a Bézier patch
    of degree {1}) for a given Bézier patch {*b} of arbitrary degree.
    The patch {t} must be allocated by the caller, with the same
    dimensions as {b} and degree 1 in each domain axis.
    
    The approximation has the same corners as the original; thus, if
    two cells share an entire face, their approximations will also
    share that face. (However, that is not necessarily true of cells
    that share only part of a face.) */

double bz_patch_multiaffine_error(bz_patch_t *b, bz_patch_t *t);
  /* Compares the general Bézier patch {b} with the multiaffine
    (degree 1) patch {t}, and returns an upper bound to the difference
    {|b(x) - t(x)|}, for any {x} in the unit cube {U}. Requires 
    {t->m == b->m} and {t->n == b->n}. */

void bz_patch_affine_approx(bz_patch_t *b, double c[], double M[], double *err);
  /* Stores in {c[0..d-1]} and {M[0..d^2-1]} an affine approximation
    {h(p)} for the Bézier patch {*b}. Also stores in {err} an upper
    bound for the approximation error {|b(x) - h(x)|}, for any {x} in
    {U}.
    
    In this approximation, an argument point {x} is mapped to point
    {h(x) = c + 2*(x-m)*M}, where {m = (1/2,.. 1/2)} is the center of
    the unit cube {U}. Note that the approximation may not preserve
    the corners of {b}.
    
    The matrix {M} is stored linearized by rows. Both {c} and {M} must
    be allocated by the caller. */

void bz_patch_grad(bz_patch_t *b, box_axis_t ax, bz_patch_t *t);
  /* Stores in {*t} the derivative of the Bézier patch {*b} along
    coordinate axis {ax}; namely, a Bézier patch from {R^m} to
    {R^{n}}, such that {t(x)[i]} is the derivative of {b(x)[i]} with
    respect to {x[ax]}.
    
    The patch {*t} must be allocated by the caller, with the same
    dimensions as {b}. Its degree {t->g} must be equal to {b->g},
    except that {t->g[ax] = max(0, b->g[ax]-1)}. */

void bz_patch_invert( 
    double *p,
    bz_patch_t *b,
    double tol,
    double x[]
  );
  /* Given a Bézier patch {b} and a point {p[0..d-1]} of {R^d},
    the procedure returns the coordinates {x[0..d-1]} 
    of the point {b^{-1}(p)}.  
    
    More precisely, the procedure finds a point {x} in {R^d} such that
    {|b(x) - p| \leq tol}. Note that {x} may lie outside the reference
    cube {U}, or even be undefined, if {p} lies outside the actual
    cell {b(U)}. */

void bz_patch_split(
    bz_patch_t *b, 
    box_axis_t a,
    double ratio,
    bz_patch_t *bLO,
    bz_patch_t *bHI
  );
  /* Returns descriptions of the two halves of a Bézier patch {b},
    that result from bissecting the unit {m}-cube {U} by a plane
    perpendicular to axis {a}, and reparametrizing each piece over the
    whole cube.
    
    The parameter {ratio} defines the position on the splitting plane:
    {0} means the face of the cube facing towards {-oo}, {1} means the
    face facing towards {+oo}, {0.5} means exact bisection. */ 

void bz_patch_print(FILE *wr, bz_patch_t *b, char *fmt);
  /* Prints the coeff matrix of {b} to file {wr},
    using the {fmt} string for each coordinate. */

/* AUXILIARY PROCEDURES */

void bz_patch_eval_m
  ( bz_patch_ddim_t m,           /* Dimension of parameter space. */
    bz_patch_rdim_t n,           /* Dimension of image space. */
    const bz_degree_t g[], /* Degree of each coordinarte. */
    double *c,           /* Bezier coeffs. */
    double x[],          /* Argument vector (size {m}). */
    double fx[]          /* OUT: Image vector (size {n}). */
  );
  /* Evaluates a Bézier patch of degree {g} from {R^m} to {R^n},
    given by the coeff array {c}, at an argument vector {x[0..m-1]}.
    The result is returned in {fx[0..n-1]}. */
  
void bz_patch_eval_1
  ( bz_patch_rdim_t n,     /* Dimension of image space. */
    bz_degree_t g,   /* Degree of curve. */
    double *c,       /* Bezier coeffs. */
    double x,        /* Argument value. */
    double fx[]      /* OUT: Image vector (size {n}). */
  );
  /* Evaluates a Bézier curve (a Bézier patch with unidimensional
    domain) of degree {g}, given by the coeff array {c}, at an
    argument {x}. The result is returned in {fx[0..n-1]}. */

void bz_patch_compute_steps(
    bz_patch_t *b, /* A Bézier patch. */
    int step[],    /* (OUT) coeff index increment for each axis. */
    int *sz       /* (OUT) total number of coeffs. */
  );
  /* Stores in each {step[i]} the increment in the total coeff index 
    {ix(e[0],..e[m-1])} that corresponds to a unit increment in the 
    index {e[i]}.  Also returns in {*sz} the total number of coeffs. */

bz_patch_t bz_patch_from_box(bz_patch_ddim_t d, interval_t box[]);
  /* Creates a Bézier patch of uniform degree 1 (multi-linear)
    with domain dimension {d} and range dimension {d},
    whose range is the box {box[0..d-1]}. */

bz_patch_t bz_patch_make_test(bz_patch_ddim_t d, bz_patch_rdim_t n, bz_degree_t g);
  /* Creates a Bézier patch of domain dimension {d}, range dimension
    {n}, and uniform degree {g}, sutable for testing. */
    
#endif
