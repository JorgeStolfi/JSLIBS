#ifndef hr2_pmap_encode_H
#define hr2_pmap_encode_H

/* Tools for minimal encoding of projective maps. */
/* Last edited on 2024-11-20 12:00:35 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <hr2_pmap.h>
#include <hr2.h>

/* Every projective map of two-dimensional oriented projective space
  {\RT^2} is uniquely determined by a 3x3 real matrix with nonzero
  determinant, with the proviso that two matrices define the same map if
  and only if thet differ by a positive scale factor. Moreover,
  infinitesimal changes in the matrix that do not change the sign of the
  determinant correspond to infinitesimal changes in the map, in the
  sense that images of finite points suffer only infinitesimal changes.
  
  Thus we can also say that the set {\RM^2} of all projective maps of
  {\RT^2} is homeomorphic to the set {A} of 3x3 real matrices with unit
  norm (whose elements squared add to 1) minus the subset {B} of
  matrices of {A} with zero determinant. The set {A} is homeomorphic to
  the sphere {\RS^8}, and the set {B} is a one-dimensional submanifold
  that separates the projective maps in two sets {\RM^2_{+}} and
  {{\RM^2_{-}}, /psotive/ and /negative/, according to the sign of the
  determinant. Maps of the first set preserve all relative positions (as
  per {hr2_side}) of points and lines, while those of the second set
  reverse them.
  
  Projective maps of any special type (such as translations, isometries,
  etc) are also divided into positive and negative ones, each set being
  a subset of {\RM^2_{+}} or {{\RM^2_{-}},  respectively, with dimension
  smaller than 8.
  
  For all types except {GENERIC}, the procedures in this interface
  define a mapping between the set of positive projective maps
  {\RN(t)_{+}} of a given type {t} and the Cartesian space {\RR^d(t)}
  where {d(t)} is the dimension of {\RN(t)_{+}}. The mapping {enc} from {\RN(t)_{+}}
  to {\RR^d(t)} (the /encoding/) is injective but may not be surjective.
  Depending on the topology of {\RN(t)_{+}} the inverse mapping {dec}
  from {\RR^d(t)} to {\RN(t)_{+}} (the /decoding/) is surjective but may
  be many-to-one. In that case, {dec} is a homomorphic covering of
  {\RN(t)_{+}}. We have {dec(enc(M)) = M} for all {M\in\RN(t)_{+}}, but
  {enc(dec(v))} may be {norm(v)} where {norm} is a branch normalization
  projection of the covering. Anyway, infinitesinal changes in {v} will
  cause infinitesimal changes in {dec(v)}.
  
  These procedures also define an involutive bijection {sym} between the
  sets {\RN(t)_{+}} and {\RN(t)_{-}} of the positive and negative maps
  of each type {t}. The bijection takes {M} to the map {Q M] where {Q} is the
  map that swaps the {x} and {y} homogeneous coordinates while preserving
  the {w} coordinate.
  
  The {encode} and {decode} procedures below implement the {enc}, {dec} and
  {sym} mappings for each specific type.  Namely, {encode} implements 
  {enc(M)} for positive maps {M}, and {enc(sym(M))} for negative ones.
  The {decode} procedure implements {dec(v)} if the {sgn} argument is {+1},
  and {sym(enc(v))} if {sgn} is {-1}. */
 
uint32_t hr2_pmap_encode_num_parameters(hr2_pmap_type_t type);
  /* Returns the dimension {d(type)} of the space
    number of parameters used to encode a projective map of the specified
    {type}. Namely, 0 for {IDENTITY} 2 for {TRANSLATION}, 3 for {CONGRUENCE},
    4 for {SIMILARITY}, and 6 for general {AFFINE} map. */

void hr2_pmap_encode(hr2_pmap_t *M, hr2_pmap_type_t type, uint32_t ny, double y[]);
  /* Converts a positive map {M} into the parameters {enc(M) =
    y[0..ny-1]} where {ny} is
    {d(type)=hr2_pmap_encode_num_parameters(type)}. Assumes that {M} is
    of the given {type}; the result is undefined otherwise. If {M} is a
    negative map, produces (enc(sym(M))} instead. Thus the handedness of
    {M} is lost. */

void hr2_pmap_decode(uint32_t ny, double y[], hr2_pmap_type_t type, sign_t sgn, hr2_pmap_t *M);
  /* The {sgn} must be {-1} or {+1}. If {sgn} is {+1}, converts the
    encoded parameters {y[0..ny-1]} into the positive projective map
    {M=dec(y)} of the given {type}. The number {ny} must be
    {d(type)=hr2_pmap_encode_num_parameters(type)}. If {sgn} is {-1},
    returns {sym(dec(y))} instead. */

#endif
