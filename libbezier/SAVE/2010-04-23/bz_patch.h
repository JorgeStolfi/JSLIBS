/* General tensor-style Bézier patches. */
/* Last edited on 2010-04-23 11:06:46 by stolfi */

#ifndef bz_patch_H
#define bz_patch_H

#include <psp_basic.h>
#include <bz_basic.h>
#include <double_array.h>

#include <interval.h>
#include <interval_array.h>
#include <box.h>
#include <stdint.h>

/* 
  BÉZIER REPRESENTATION FOR MULTIVALUED MULTIVARIATE POLYNOMIALS
  
  A real-valued polynomial on {d} variables, which has maximum degree
  {g[i]} on coordinate {i}, for each {i} in {0..d-1}, can be specified
  by its /Bézier coefficients/ relative to box {B} of {R^d}. These are
  an array of real values with {g[i]+1} points along each axis {i}.
  
  The Bézier coefficients are conceptually associated with a regular
  lattice of points spanning the box {B}, Each Bézier coefficient is
  specified by a tuple {e[0..d-1]} of {d} /Bézier indices/, where each
  index {e[i]} ranges over {0,.. g[i]}. The coefficient with index
  tuple {e[0..d-1]} multiplies the /Bernstein-Bézier tensor element/
  
    { PROD { BB^g_{e[i]}(z[i]) : i in {0,..d-1} } }
    
  where {z[0..d-1]} are the coordinates of the point relative
  to the tile {B}; and {BB^g_e(z)} is the /Bernstein-Bézier
  polynomial/ of degree {g} and index {e}. 
  
  The same representation holds for multi-valued polynomials, 
  except that each coefficient is not a single number but an array
  of numbers, with an arbitrary number of indices {r} and an arbitrary
  number of elements along each index.
  
  For example, a Bézier surface patch in 3D (as widely used in
  computer graphics) can be represented as a {bz_patch_t} {b} with
  domain dimension {d=2}. An argument {x} for {b} is a vector with two
  real numbers, the surface parameters {u} and {v}. Each value {b(x)}
  of the patch, and each Bézier coefficient, is a point of 3D space;
  that is, a vector (single-index array) with 3 elements, the spatial
  coordinates {X,Y,Z}.
  
  The gradient of such a surface patch can be represented by a {bz_patch_t}
  {db}, also with domain dimension {d=2}.  Each value {db(x)} of {db}
  is a {2 × 3} matrix containing the partial derivatives of the components
  of {b} at {x}; namely, {((dX/du, dY/du, dZ/du), (dX/dv, dY/dv, dZ/dv))}.
  */

typedef psp_dim_t bz_patch_ddim_t;
  /* Dimension of the domain of a Bézier patch. */

typedef psp_dim_t bz_patch_rinds_t;
  /* Dimension of the range of a Bézier patch. */

#define bz_patch_MAX_DDIM (double_array_MAX_AXES)
  /* Maximum dimension for the domain of a Bézier patch. Note that
    some vectors are preallocated with this size, and that the number
    of coeffs is exponential in the domain dimension. */

#define bz_patch_MAX_RINDS (double_array_MAX_AXES)
  /* Maximum number of axes (indices) in the values of a Bézier patch. */

#define bz_patch_MAX_DEGREE (6)
  /* Maximum degre of a Bézier patch along each coordinate,
    just to be safe. Could be increased with modest harm. */

typedef bz_degree_t bz_patch_degrees_t[bz_patch_MAX_DDIM];
  /* The degree of a Bézier patch, along each domain axis. */

typedef bz_degree_t bz_patch_cpindex_t;
  /* Each coeff has a {bz_patch_cpindex_t} along each domain axis,
    which ranges from 0 to the corresponding degree, inclusive. */

typedef struct bz_patch_t
  { bz_patch_ddim_t d; /* Dimension of domain. */  
    double_array_t C;  /* Bernstein-Bézier coefficients. */
  } bz_patch_t;
  /* A {bz_patch_t} record {b} describes a polynomial map that 
    to each point {x[0..d-1]} of {R^d} associates some array {b(x)}
    of real numbers; where each of these numbers is a polynomial 
    on the coordinates {x}.
    
    The map {b} is defined as as a linear combination of Bézier
    elements, which are products of {d} univariate Bernstein-Bézier
    polynomials. Factor {i} of this product is a polynomial of degree
    {g[i]} on coordinate {i} of the argument. The coefficients of this
    linear combination are the /Bézier coefficients/ of the patch,
    /coeffs/ for short. Each coefficient is a slice
    of the {C} array. The first {d} indices of the {C} array
    select a coefficient, the rest are indices of the coefficient.
    
    The coeffs are logically organized as an array with {d} subscripts
    {e[0..d-1]}. Each subscript {e[i]} ranges over {0..g[i]},
    therefore the total number {nCoefs} of coeffs is {PROD{(g[i]+1) :
    i = 0..d-1}.
    
    The coeffs are stored in an array {b.C} with {d+r} indices, for
    some {r} greater than 0. The Bézier coeff {b.coef[e]} with indices
    {e = e[0..d-1]} is the subarray of {b.C} with {r} indices,
    consisting of all elements whose first {d} indices are
    {e[0..d-1]}. Namely, if {A=b.coef[e[0],..e[d-1]]}, then for any
    index tuple {f[0..r-1]} we have {A[f[0],..f[r-1]] =
    b.C[e[0],..e[d-1],f[0],..f[r-1]}.
    
    Therefore, the value of {b(x)} on a point {x} of {R^d} is the sum
    of
    
      { b.coef[e] * PROD{ B_{e[i]}^{g[i]}(x[i]) : i = 0..d-1 } }
    
    for all Bézier coeffs {b.coef[e]}, where {B_k^r(z) ==
    z^k*(1-z)^{r-k}*\choose(r,k)/2^r} is the Bernstein-Bezier
    polynomial of index {k} and degree {r}. Of course, the value of
    {b(x)} can be computed more efectively by successive
    interpolations along each coordinate (the DeCasteljau algorithm).
    
    Although the map {b()} is defined over all {R^d}, for many
    applications it is convenient to define the patch proper as the
    restriction of {b} to the unit {d}-cube {U = [0 _ 1]^d} of {R^d}.
    In particular, a Bezier coefficient {b.coef[e]} where every
    subscript {e[i]} is ether {0} or {g[i]} is a `corner' of the patch
    (the image of a corner of {U}). The remaining coeffs (which exist
    only if some {g[i]} is greater than 1) affect the shape of the
    patch between those corners. */
/*
  INSPECTION
  
  In all the following, {b.d} denotes the dimension of the domain of
  the patch {b}; that is, the number of indices in its Bézier
  coefficient array. Also {b.r} denotes the number of indices in each
  value {b(x)} of {b}, or in each Bézier coefficient; that is, the
  number of indices of the array {b.C}, minus {b.d}. */

bz_patch_ddim_t bz_patch_get_ddim(bz_patch_t *b);
  /* Returns the dimension of the domain of patch {b},
    namely {b.d}. */

bz_patch_rinds_t bz_patch_get_rinds(bz_patch_t *b);
  /* Returns the number of indices of the values (and coefficients) of
    patch {b}, namely the number of indices of {b.C} minus {b.d}. */

void bz_patch_get_degrees(bz_patch_t *b, bz_degree_t g[]);
  /* For each {i} in {0..d-1}, where {d = bz_patch_get_ddim(b)}, 
    stores into {g[i]} the degree of the patch along
    the domain axis number {i}. */

ix_size_t *bz_patch_get_coeff_sizes(bz_patch_t *b);
  /* Returns the address of an array {dsz[0..d-1]}, where 
    {d = bz_patch_get_ddim(b)}, such that, for each {i} in {0..d-1},
    {dsz[i]} is the number of Bézier coefficients along axis {i} of the
    domain (which is 1 more than the degree of the patch on coordinate
    {i} of the argument). */

ix_size_t *bz_patch_get_range_sizes(bz_patch_t *b);
  /* Returns the address of an array {rsz[0..r-1]}, where 
    {r = bz_patch_get_rdim(b)}, such that, for for each {j} in {0..r-1},
    {rsz[j]} is the number of components along axis {i} of the
    patch's value at a singl epoint of the domain. */
    
/*
  ALLOCATION */

bz_patch_t bz_patch_new(bz_patch_ddim_t d, const bz_degree_t g[], bz_patch_rinds_t r, const ix_size_t rsz[]);
  /* Returns a Bézier patch {b} with {d}-dimensional domain and
    degrees {g[0..d-1]} along the {d} axes of the domain. The value of
    the patch at any point, as well as each Bézier coefficient, will
    be an array with {r} indices and {rsz[j]} elements along each
    index. If {rsz} is NULL, assumes {rsz[j] = 1} for all {j}. The
    coefficient array {b.C} will be newly allocated. */

bz_patch_t bz_patch_uniform_new(bz_patch_ddim_t d, bz_degree_t g, bz_patch_rinds_t r, const ix_size_t rsz[]);
  /* Returns a Bézier patch {b} with {d}-dimensional domain and the
    same degree {g} along each axis of the domain. The parameters {r}
    and {rsz[0..r-1]} have the same meaning as in {bz_patch_new}. The
    coefficient array {b.C} will be newly allocated. */

void bz_patch_free(bz_patch_t *b);
  /* De-allocates the coefficient vector and other internal storage
    areas of {b} (but not {b} itself). */

double *bz_patch_control_point(bz_patch_t *b, bz_patch_cpindex_t e[]);
  /* Given a Bézier patch {b} of degree {g = b->g} from {R^d} to {R},
    and a vector {e[0],..e[d-1]} of indices, where each {e[i]} lies in
    {0..b->g[i]}, returns the address of the coeff
    {b.coef[e]}. Note that the coeff is neither
    allocated nor copied --- the result is just a pointer into the
    control array {*(b->C)}.. */

void bz_patch_eval(bz_patch_t *b, double x[], double_array_t bx);
  /* Stores into the array {bx} the value {b(x)} of the  
    patch {b} at the point {x[0..d-1]} of {R^b.d}. 
    The array {bx} must be allocated by the caller, 
    with . */

void bz_patch_eval_box(bz_patch_t *b, interval_t box[], bz_patch_t *f);
  /* Stores in {*f} a Bézier patch that describes the restriction
    of patch {b} to of the given axis-aligned {box[0..d-1]},
    reparametrized to the unit cube {U} of {R^d}. The patch
    {*f} must be allocated by the caller, with the same dimensions
    and degree as {b}. */

void bz_patch_compute_bbox(bz_patch_t *b, interval_array_t *box);
  /* Stores in each element of the array {box} an interval
    that is guaranteed to contain the corresponding element of
    {b(x)}, when the argument point {x} ranges over
    the unit cube {U} of {R^d}. */

/* 
   FACES AND FACETS */

bz_patch_t bz_patch_facet_new(bz_patch_t *b, box_axis_index_t j);
  /* Returns a Bézier patch {f} suitable for storing into it a facet
    of another Bézier patch {b} that is perpendicular to axis
    {j} of {b}'s domain. Thus, the new patch {f} will have domain
    dimension {d-1}, the same range space as {b}, and its degree vector will
    be {g[0..d-1]} with element {g[j]} deleted. The coefficients
    {f.C} are newly allocated. */

void bz_patch_get_facet(bz_patch_t *b, box_axis_t ax, interval_side_t dir, bz_patch_t *t);
  /* Given a Bézier patch {b} with domain {D} of dimension {d}, stores
    in {t} the restriction of {b} to the {d-1}-dimensional face
    (facet) of {D} that is perpendicular to axis {ax}, and is 
    located in the {dir} direction of that axis with respect to {D}. 
    
    The patch {t} must be allocated by the caller with proper
    parameters, e.g. through {bz_patch_facet_new}. */

void bz_patch_get_face
  ( bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  );
  /* Stores in {t} a Bezier patch that describes the shape of an 
    {s}-dimensional face {F} of patch {b}.  The face is selected
    by its signature vector {dir[i]}, as in {box.h}.
    
    The patch {t} must be allocated by the caller, with range dimension
    {t->r = b->r}. Its domain dimension {t->d} must be the 
    dimension of the face, i.e. the number of elements {dir[i]}
    which are equal to {SMD}; and its degree vector {t->g}
    must contain the elements of {b->g} that correspond to 
    those axes. */

void bz_patch_set_face
  ( bz_patch_t *b, 
    box_signed_dir_t dir[], 
    bz_patch_t *t
  );
  /* Assumes that {*t} is a Bezier patch that describes the shape of an 
    {s}-dimensional face {F} of patch {b}, specified by its signature {dir}. 
    Modifies {b} so that the face {F} has that shape.  The remaining control
    points of {b} are modified too, proportionally to their distance from 
    face {F}. */
    
typedef void bz_patch_visit_t(bz_patch_t *b, interval_t box[]);
  /* Type of a sub-patch visiting function. It receives a patch {b} and
    its corresponding domain box {box[0..b.d-1]}. The {box}
    may be null. */

void bz_patch_enum_faces
  ( bz_patch_t *b, 
    interval_t box[], 
    bz_patch_ddim_t df
  );
  /* Enumerates all faces of {b} with dimension {d}, in order of
    increasing index, and the corresponding faces of the box
    {box[0..df-1]}. Each face will be a patch {bf} with the same
    domain dimension as {b}, but with degree 0 (constant) along each
    axis perpendicular to the face. Each bace of {box} is an array of
    {b.d} intervals {boxf[0..b.d-1]} which are trivial along those
    same axes. Calls {proc(bf, boxt)} for each face {bf} and its
    corrsponding subdomain {boxf}. The {box} may be null, in which
    case each {boxf} will be null. */
 
double bz_patch_try_flatten_face
  ( bz_patch_t *b, 
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
    {t->d == b->d} and {t->r == b->r}. */

void bz_patch_affine_approx(bz_patch_t *b, double c[], double D[], double *err);
  /* Stores in {c[0..d-1]} and {D[0..d^2-1]} an affine approximation
    {h(r)} for the Bézier patch {*b}. Also stores in {err} an upper
    bound for the approximation error {|b(x) - h(x)|}, for any {x} in
    {U}.
    
    In this approximation, an argument point {x} is mapped to point
    {h(x) = c + 2*(x-d)*D}, where {d = (1/2,.. 1/2)} is the center of
    the unit cube {U}. Note that the approximation may not preserve
    the corners of {b}.
    
    The matrix {D} is stored linearized by rows. Both {c} and {D} must
    be allocated by the caller. */

void bz_patch_grad(bz_patch_t *b, box_axis_t ax, bz_patch_t *t);
  /* Stores in {*t} the derivative of the Bézier patch {*b} along
    coordinate axis {ax}; namely, a Bézier patch from {R^d} to
    {R^{r}}, such that {t(x)[i]} is the derivative of {b(x)[i]} with
    respect to {x[ax]}.
    
    The patch {*t} must be allocated by the caller, with the same
    dimensions as {b}. Its degree {t->g} must be equal to {b->g},
    except that {t->g[ax] = max(0, b->g[ax]-1)}. */

void bz_patch_invert
  ( double *p,
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

void bz_patch_split
  ( bz_patch_t *b, 
    box_axis_t a,
    double ratio,
    bz_patch_t *bLO,
    bz_patch_t *bHI
  );
  /* Returns descriptions of the two halves of a Bézier patch {b},
    that result from bissecting the unit {d}-cube {U} by a plane
    perpendicular to axis {a}, and reparametrizing each piece over the
    whole cube.
    
    The parameter {ratio} defines the position on the splitting plane:
    {0} means the face of the cube facing towards {-oo}, {1} means the
    face facing towards {+oo}, {0.5} means exact bisection. */ 

void bz_patch_enum_sub
  ( bz_patch_t *b, 
    interval_t box[],
    bz_patch_visit_t *proc,
    int minRank,
    int maxRank,
    double tol
  );
  /* 
    Enumerates a set of sub-patches of a Bézier patch {bz} with domain
    box {B[0..b.d-1]}. 
    
    The Bézier patch {bz} must be non-null NULL.
    
    The domain rectangle {B[0..b.d-1]} is subdivided recursively in
    half along each axis in turn, cyclically, at least {minRank}
    times. The subdivision continues until {maxRank} splits, or until
    the values of the patch, restricted to the sub-box, deviate less
    than {tol} from the multiaffine patch with the same corners. The
    function {proc} is called for each leaf (unsplit) sub-box. */

void bz_patch_print(FILE *wr, bz_patch_t *b, char *fmt);
  /* Prints the coeff matrix of {b} to file {wr},
    using the {fmt} string for each coordinate. */

/* AUXILIARY PROCEDURES */

void bz_patch_eval_d
  ( bz_patch_ddim_t d,     /* Dimension of parameter space. */
    bz_patch_rinds_t r,    /* Dimension of image space. */
    const bz_degree_t g[], /* Degree of each coordinarte. */
    double *c,             /* Bezier coeffs. */
    double x[],            /* Argument vector (size {d}). */
    double fx[]            /* OUT: Image vector (size {r}). */
  );
  /* Evaluates a Bézier patch of degree {g} from {R^d} to {R^r},
    given by the coeff array {c}, at an argument vector {x[0..d-1]}.
    The result is returned in {fx[0..r-1]}. */
  
void bz_patch_eval_1
  ( bz_patch_rinds_t r,     /* Dimension of image space. */
    bz_degree_t g,   /* Degree of curve. */
    double *c,       /* Bezier coeffs. */
    double x,        /* Argument value. */
    double fx[]      /* OUT: Image vector (size {r}). */
  );
  /* Evaluates a Bézier curve (a Bézier patch with unidimensional
    domain) of degree {g}, given by the coeff array {c}, at an
    argument {x}. The result is returned in {fx[0..r-1]}. */

void bz_patch_compute_steps
  ( bz_patch_t *b, /* A Bézier patch. */
    int step[],    /* (OUT) coeff index increment for each axis. */
    int *sz       /* (OUT) total number of coeffs. */
  );
  /* Stores in each {step[i]} the increment in the total coeff index 
    {e(e[0],..e[d-1])} that corresponds to a unit increment in the 
    index {e[i]}.  Also returns in {*sz} the total number of coeffs. */

bz_patch_t bz_patch_from_box(bz_patch_ddim_t d, interval_t box[]);
  /* Creates a Bézier patch of uniform degree 1 (multi-linear)
    with domain dimension {d} and range dimension {d},
    whose range is the box {box[0..d-1]}. */

bz_patch_t bz_patch_make_test(bz_patch_ddim_t d, bz_patch_rinds_t r, bz_degree_t g);
  /* Creates a Bézier patch of domain dimension {d}, range dimension
    {r}, and uniform degree {g}, sutable for testing. */
    
#endif
