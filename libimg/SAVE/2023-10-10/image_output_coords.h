/* User coordinate system and dimensions for output images. */
/* Last edited on 2023-10-10 08:44:06 by stolfi */

#ifndef image_output_coords_H
#define image_output_coords_H

#include <bool.h>
#include <r2.h>
#include <argparser.h>

/* OUTPUT USER UNIT */

void imgc_parse_output_unit(argparser_t *pp, double *oUnit);
  /* Parses the command line argument "-oUnit" according to the syntax
    specs {imgc_parse_output_unit_HELP}, and returns it in {*oUnit}. If
    \"-oUnit\" is not present, sets it to 1. */

#define imgc_parse_output_unit_HELP \
  "[ -oUnit {OUNIT} ]"

#define imgc_parse_output_unit_INFO_OPTS \
  "  -oUnit {OUNIT}\n" \
  "    This optional argument specifies the size in pixels of " \
  " the user coordinate unit for output images.  If this argument is not specified, the program" \
  " assumes \"-oUnit 1\" (which means that output user units are in pixels)."
  
#define imgc_output_unit_affects_output_org_INFO_OPTS \
  "The option \"-oUnit\" affects the interpretation of \"-oOrg\".  Thus, given \"-oUnit 0.5\", the" \
  " option \"-oOrg 24 30\" will set the user output image origin at 12 pixel columns" \
  " and 15 pixel rows away from the default origin."
  
#define imgc_output_unit_affects_output_size_INFO_OPTS \
  "The option \"-oUnit\" affects the interpretation" \
  " of \"-oSize\".  Thus, given \"-oUnit 2.0\", the" \
  " option \"-oSize 30 40\" will set the actual output" \
  " image size, in pixels, to 60 columns and 80 rows."
  
#define imgc_input_unit_affects_default_output_size_INFO_OPTS \
  "The parameter \"-iUnit\" affects the default output" \
  " image size.  Suppose the input image has 200 by 300 pixels.  Given the" \
  " argument \"-iUnit 0.5\", the default output image size would be 400 by 600" \
  " /user/ units.  Given also \"-oUnit 3.0\", the actual default output" \
  " image size, in pixels, would then be 1200 by 1800."

/* OUTPUT IMAGE ORIGIN */

void imgc_parse_output_center_org(argparser_t *pp, bool_t *oCenter, r2_t *oOrg);
  /* Analogous to {imgc_parse_center_org}, but applying to output
    images, according to {imgc_parse_output_center_org_INFO_OPTS}. 
    Looks for the keywords "-oCenter" and "-oOrg". */

#define imgc_parse_output_center_org_HELP \
  "[ -oCenter | -oOrg {CX_OUT} {CY_OUT} ]"

#define imgc_output_origin_INFO \
  "The default coordinate system origin for output images" \
  " can be overriden by the \"-oOrg\" or \"-oCenter\" argument."

#define imgc_parse_output_center_org_INFO_OPTS(org_default) \
  "  -oCenter\n" \
  "  -oOrg {CX_OUT} {CY_OUT}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the user coordinate system for output" \
  " images.  The \"-oCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-oOrg\" option sets" \
  " the origin displaced by {CX_OUT} and {CY_OUT} from the" \
  " default origin, in the directions specificed by the \"-xAxis\" and \"-yAxis\" options." \
  "  " org_default
  
#define imgc_parse_output_center_org_INFO_OPTS_default_zero \
  "If neither option is specified, the program" \
  " assumes \"-oOrg 0 0\" (which means the default origin)."
  
#define imgc_parse_output_center_org_INFO_OPTS_default_center \
  "If neither option is specified, the program" \
  " assumes \"-oCenter\" (the center of the image)."

/* OUTPUT IMAGE SIZE */

void imgc_parse_output_size(argparser_t *pp, double *oCols, double *oRows);
  /* Parses the output image size option "-oSize", according to
    {imgc_parse_output_size_HELP}. Namely, if "-oSize" is present, sets
    {*oCols} and {*oRows} to the next two arguments. Otherwise leaves
    {*oCols} and {*oRows} unchanged. */
    
#define imgc_parse_output_size_HELP \
  "[ -oSize {NX} {NY} ]"

#define imgc_parse_output_size_INFO_OPTS(osize_default) \
  "  -oSize {NX} {NY}\n" \
  "    Specifies the horizontal and vertical size of the" \
  " output image, in output /user/ units.  These values may be" \
  " fractional. The actual image size in pixels will be rounded up" \
  " to the next integer." \
  "  " osize_default

#define imgc_parse_output_size_INFO_OPTS_default_input \
  "If \"-oSize\" is omitted, assumes that the output image" \
  " size (in output /user/ coords units) is the same as the input" \
  " image size (in input /user/ coords units)."

#endif
