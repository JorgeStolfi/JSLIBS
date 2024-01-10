/* User coordinate system for input images. */
/* Last edited on 2023-08-28 14:36:15 by stolfi */

#ifndef image_input_coords_H
#define image_input_coords_H

#include <bool.h>
#include <r2.h>
#include <argparser.h>
  
/* INPUT USER UNIT */

void imgc_parse_input_unit(argparser_t *pp, double *iUnit);
  /* Parses the command line argument "-iUnit" according to the syntax
    specs {imgc_parse_input_unit_HELP}, and returns it in {*iUnit}. If
    \"-iUnit\" is not present, sets it to 1. */

#define imgc_parse_input_unit_HELP \
  "[ -iUnit {IUNIT} ]"

#define imgc_parse_input_unit_INFO_OPTS \
  "  -iUnit {IUNIT}\n" \
  "    This optional argument specifies the size in pixels of " \
  " the user coordinate unit for input images.  If this argument is not specified, the program" \
  " assumes \"-iUnit 1\" (which means that user units are in pixels)."
  
#define imgc_input_unit_affects_input_org_INFO_OPTS \
  "The option \"-iUnit\" affects the interpretation of \"-iOrg\" as" \
  " well as other coordinate parameters for input images.  Thus, given \"-iUnit 0.5\", the" \
  " option \"-iOrg 24 30\" will set the user input image origin at 12 pixel columns" \
  " and 15 pixel rows away from the default origin."

/* INPUT IMAGE ORIGIN */

void imgc_parse_input_center_org(argparser_t *pp, bool_t *iCenter, r2_t *iOrg);
  /* Parses the mutually exclusive arguments "-iCenter" and "-iOrg",
    that define the coordinate origin for input images,
    according to the syntax specs {imgc_parse_input_center_org_HELP}. If
    "-iCenter"is preent, sets {*iCenter} to TRUE and leaves {*iOrg}
    unchanged. If "-iOrg" is present, sets {*iCenter} to FALSE, and sets
    {*iOrg} to the given point. If neither is present, leaves {*iCenter}
    and {*iOrg} unchanged. */

#define imgc_parse_input_center_org_HELP \
  "[ -iCenter | -iOrg {CX_IN} {CY_IN} ]"

#define imgc_input_origin_INFO \
  "The default coordinate system origin for input images" \
  " can be overriden by the \"-iOrg\" or \"-iCenter\" argument."  

#define imgc_parse_input_center_org_INFO_OPTS(iorg_default) \
  "  -iCenter\n" \
  "  -iOrg {CX_IN} {CY_IN}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the user coordinate system for input" \
  " images.  The \"-iCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-iOrg\" option sets" \
  " the origin displaced {CX_IN} and {CY_IN} /user/ units from the" \
  " default origin, in the directions specified by \"-xAxis\" and \"-yAxis\"." \
  "  " iorg_default
  
#define imgc_parse_input_center_org_INFO_OPTS_default_zero \
  "If neither option is specified, the program" \
  " assumes \"-iOrg 0 0\" (which means the default origin)."
  
#define imgc_parse_input_center_org_INFO_OPTS_default_center \
  "If neither option is specified, the program" \
  " assumes \"-iCenter\" (the center of the image)."

#endif
