/* argparser_geo.h -- extends argparser.h for geometric args. */
/* Last edited on 2024-11-20 08:47:52 by stolfi */

#ifndef argparser_geo_H
#define argparser_geo_H

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */

/* This interface provides convenient tools for parsing command
  line arguments whose values are real vectors. */

#define _GNU_SOURCE
#include <stdint.h>
#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r6.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <argparser.h>

r2_t argparser_get_next_r2(argparser_t *pp, double min, double max);
r3_t argparser_get_next_r3(argparser_t *pp, double min, double max);
r4_t argparser_get_next_r4(argparser_t *pp, double min, double max);
r6_t argparser_get_next_r6(argparser_t *pp, double min, double max);
  /* Parses with {argparser_get_next_double} next {N} arguments, which
    must be in the range {[min.. max]}, as the coordinates of a point;
    where {N} is 2,3,4, or 6. */

void argparser_get_next_rn(argparser_t *pp, double p[], uint32_t n, double min, double max);
  /* Parses with {argparser_get_next_double} the next {n} arguments, 
    which must be in the range {[min .. max]}, and stores them in 
    {p[0.n-1]}. */

r3_t argparser_get_next_r3_dir(argparser_t *pp);
  /* Same as {argparser_get_next_r3} but normalizes result to unit length. */

void argparser_get_next_adjust(argparser_t *pp, double *adjP, double min, double max);
  /* Parses the optional modifier "adjust {AMOUNT}" to the 
    most recently parsed command line argument. 
    
    If {adjP} is NULL, does nothing. If {adjP} is not NULL, and the
    next argument is "adjust", parses that keyword and the following
    (required) numeric argument {AMOUNT}, which must be in the range
    {[min _ max]}, and stores that value into {*adjP}. If the "adjust"
    keyword is not present, the procedure parses nothing and leaves
    {*adjP} unchanged. */
        
/* ---------------------------------------------------------------------- */
/* PROJECTIVE MAP TYPE */

hr2_pmap_type_t argparser_get_next_hr2_pmap_type(argparser_t *pp);
  /* Parses the next argument as a {hr2_pmap_type_t} value.
    It should be a string in all-uppercase or all lowercase:
    "IDENTITY" or "identity" for {hr2_pmap_type_IDENTITY},
    "TRANSLATION" or "translation" for {hr2_pmap_type_TRANSLATION},
    and so on.  Fails with error message if the next argument does
    not exist, has been parsed, or is not one of the valid types. */

/* ---------------------------------------------------------------------- */
/* PROJECTIVE MAP AS MATRIX */

hr2_pmap_t argparser_get_next_proj_map_matrix(argparser_t *pp);
  /* Parses the next nine arguments, which must be floating point numbers,
    as the elements of a 3×3 projetive transformation matrix, as 
    per {argparser_proj_map_matrix_HELP} and {argparser_proj_map_matrix_INFO}. */
    
#define argparser_proj_map_matrix_HELP \
  "{D} {TX} {TY}  {DX} {UX} {VX}  {DY} {UY} {VY}"
  
#define argparser_proj_map_matrix_INFO \
  "specifies the elements of the 3×3 projetive transformation" \
  " matrix, in row by row order."
    
/* ---------------------------------------------------------------------- */
/* PROJECTIVE MAP AS FOUR POINT PAIRS */

hr2_pmap_t argparser_get_next_proj_map_from_points(argparser_t *pp);
  /* Parses the next sixteen arguments, which must be floating point
    numbers, as the coordinates of the vertices of two quadrilaters
    {Q1,Q1}, and returns the projective map that takes {Q1} to {Q2}, as
    per {argparser_proj_map_from_points_HELP} and
    {argparser_proj_map_from_points_INFO}. */
    
#define argparser_proj_map_from_points_HELP \
  "        {XIn[1]} {YIn[1]} \\\n" \
  "        {XIn[2]} {YIn[2]} \\\n" \
  "        {XIn[3]} {YIn[3]} \\\n" \
  "        {XIn[4]} {YIn[4]} \\\n" \
  "        \\\n" \
  "        {XOut[1]} {YOut[1]} \\\n" \
  "        {XOut[2]} {YOut[2]} \\\n" \
  "        {XOut[3]} {YOut[3]} \\\n" \
  "        {XOut[4]} {YOut[4]}"
  
#define argparser_proj_map_from_points_INFO \
  "states that the projective transformation must" \
  " take each point {(XIn[i],YIn[i])} of the input image to" \
  " the corresponding point {(XOut[i],YOut[i])} of the output" \
  " image.  No three of the input-side points may be collinear," \
  " and ditto for the output-side points.  The input and output" \
  " coordinates are relative to the user input and output coordinate" \
  " systems, respectively.  The matrix {M} is computed using" \
  " {hr2_pmap_from_four_points} in {hr2.h}"
  
/* ---------------------------------------------------------------------- */
/* PROJECTIVE MAP AS MATRIX OR POINT PAIRS */

hr2_pmap_t argparser_get_proj_map(argparser_t *pp);
  /* Parses an optional description of a 2D projective map from the command line, as
    per {argparser_proj_map_HELP} and
    {argparser_proj_map_INFO}. */
    
#define argparser_proj_map_HELP \
  "[ -matrix " argparser_proj_map_matrix_HELP " \\\n" \
  "    | \\\n" \
  "      -points \\\n" \
  "" argparser_proj_map_from_points_HELP " \\\n" \
  "    ]"
    
#define argparser_proj_map_INFO \
  "The projective map {M} may be specified as a {3×3}" \
  " homogeneous coefficient array (with" \
  " the \"-matrix\" command line argument).  Alternatively, one" \
  " may give four points in the input image and" \
  " four corresponding points in the output image (with \"-points\")"

#define argparser_proj_map_HELP_INFO \
  "  -matrix " argparser_proj_map_matrix_HELP "\n" \
  "    This optional argument " argparser_proj_map_matrix_INFO "This argument is mutually" \
  " exclusive with \"-points\".  If \"-matrix\"" \
  " and \"-points\" are both omitted,  the" \
  " projective transformation is trivial (the identity map).\n" \
  "\n" \
  "  -points \\\n" \
  "" argparser_proj_map_from_points_HELP "\n" \
  "    This optional argument " argparser_proj_map_from_points_INFO "   This argument is mutually" \
  " exclusive with \"-matrix\".   If \"-points\" and" \
  " \"-matrix\" are both omitted, the" \
  " projective transformation is trivial (the identity map)." \
  
/* ---------------------------------------------------------------------- */
/* AFFINE MAP FOR FEATURE EXTRACTION AND MATCHING */

hr2_pmap_t argparser_get_next_feature_map(argparser_t *pp);
  /* Parses the next command line arguments as an affine map that takes a canonical 
    orthogonal sampling grid to a grid of samples on {R^2}, as described by 
    {argparser_feature_map_HELP} and {argparser_feature_map_INFO}. */
    
#define argparser_feature_map_HELP \
  " {CX} {CY}  {A00} {A01} {A10} {A11} / {MAG}"
  
#define argparser_feature_map_INFO \
  "These parameters specify an affine map" \
  " that takes the origin to the point {(CX,CY)} and the unit" \
  " displacement vectors {(1,0)} and {(0,1)} to {(A00,A01)/MAG} and" \
  " {(A10,A11)/MAG}, respectively."

/* Copyright © 2003 by Jorge Stolfi.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/

#endif
