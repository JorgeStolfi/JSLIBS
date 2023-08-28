/* Coordinate systems for images. */
/* Last edited on 2023-08-28 02:25:44 by stolfi */

#ifndef image_coords_H
#define image_coords_H

#include <bool.h>
#include <r2.h>
#include <r3x3.h>
#include <hr2.h>
#include <argparser.h>

hr2_pmap_t imgc_coord_sys_map
  ( bool_t xRev, 
    bool_t yRev, 
    double unit,
    bool_t center, 
    r2_t *org, 
    int cols, 
    int rows
  );
  /* Returns a projective map {M} that converts point coordinates from
    some /pixel coordinate system/ (PS), based on pixel indices, to a /user
    coordinate system/ (US).
    
    The PS has the first (X) coordinate increasing with column index,
    the second (Y) coordinate increasing with row index. The image
    domain is the rectangle {[0 _ cols] × [0 _ rows]} where {cols,rows} are
    the image dimensions in pixels. The pixel with indices
    {x} and {y} is assumed to be the unit square with corners {(x,y)} and
    {(x+1,y+1)} in pixel coordinates.
    
    In the US, the X axis increases with the column index if {xRev} is
    FALSE, and opposite to the column index if {xRev} is TRUE. The Y
    axis increases with row index if {yRev} is FALSE, and opposite to
    row index if {yRev} is true. 
    
    The unit of the {US}, on both axes, is assumed to be {unit} times
    the pixel size.
    
    The US origin is located at the center of the image's domain if
    {center} is TRUE, or at the position {org} if {center} is FALSE. In
    the second case, each coordinate of {org} is measured (in user
    units) from the edge of the image's domain that has lowest user
    coordinate (that is, from pixel X coordinate {X=0} if {xRev} is
    FALSE or {X=cols} if {xRev} is true; from pixel Y coordinate {Y=0}
    if {yRev} is FALSE, {Y=rows} if {yRev} is TRUE). */

/* HELP AND INFO TEXTS */

/* The following are standard texts to include in the {PROG_HELP} and {PROG_INFO} strings in 
  programs with let the user define user coordinate systems on both input and output, with 
  variable units. */

#define imgc_input_output_coords_HELP \
  "    " imgc_parse_x_axis_HELP " \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
  "    " imgc_parse_input_unit_HELP " \\\n" \
  "    " imgc_parse_input_center_org_HELP " \\\n" \
  "    " argparser_proj_map_HELP " \\\n" \
  "    " imgc_parse_output_unit_HELP " \\\n" \
  "    " imgc_parse_output_center_org_HELP " \\\n" \
  "    " imgc_parse_output_size_HELP

#define imgc_input_output_coords_intro_INFO \
  "  " imgc_user_axes_INFO "" \
  "  " imgc_pixel_axes_INFO "" \
  "  " imgc_pixel_centers_INFO "" \
  "  " imgc_input_origin_INFO "" \
  "  " imgc_output_origin_INFO

#define imgc_parse_input_output_coords_INFO_OPTS \
  imgc_parse_x_axis_HELP_INFO \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments, including" \
  " {CX_IN}, {CX_OUT}, {X_IN[i]}, {X_OUT[i]}, and the matrix" \
  " coefficients." \
  "  " imgc_parse_x_axis_pbm_default_INFO "\n" \
  "\n" \
  imgc_parse_y_axis_HELP_INFO "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments, including" \
  " {CY_IN}, {CY_OUT}, {Y_IN[i]}, {Y_OUT[i]}, and the matrix" \
  " coefficients." \
  "  " imgc_parse_y_axis_pbm_default_INFO "\n" \
  "\n" \
  imgc_parse_input_unit_HELP_INFO "\n" \
  "\n" \
  imgc_parse_input_center_org_HELP_INFO \
  "  " imgc_input_unit_affects_org_INFO "\n" \
  "\n" \
  imgc_parse_output_unit_HELP_INFO "\n" \
  "\n" \
  imgc_parse_output_size_HELP_INFO \
  "  " imgc_output_unit_affects_output_size_INFO \
  "  " imgc_output_size_defaults_to_input_INFO \
  "  " imgc_input_unit_affects_default_output_size_INFO "\n" \
  "\n" \
  imgc_parse_output_center_org_HELP_INFO \
  "  " imgc_output_unit_org_INFO

/* INDIVIDUAL HELP AND INFO TEXTS */

/* The texts below are elements used to compose the {PROG_HELP}, 
  {PROG_INFO_DESC}, and {PROG_INFO_OPTS} of programs that may 
  give the user moer limited options when defining input and output
  coordinate systems and the output image sizes. */

#define imgc_user_axes_INFO \
  "Points in image domains" \
  " are expressed in terms of a /user coordinate system/" \
  " specific to each image, input or output.  In both systems," \
  " the first (X) axis is horizontal and the second (Y) is vertical.  The" \
  " directions of the axes are determined by" \
  " the \"-xAxis\" and \"-yAxis\" command line arguments.  These" \
  " arguments also affect the default origin of the user coordinate" \
  " system.  In the default user coordinate systems, the origin lies" \
  " at the low corner of the domain (which depends on the axis directions)."

#define imgc_pixel_axes_INFO \
  "Internally, positions on each image are identified by" \
  " a /pixel coordinate system/" \
  " in which the X axis points right, the Y axis points" \
  " down, te origin is at the upper left corner, and" \
  " the unit is the pixel spacing."
  
#define imgc_pixel_centers_INFO \
  "In pixel coordinate systems, pixel *corners* have integer" \
  " coordinates, and pixel *centers* have half-integer coordinates."
  
/* Horizontal axis direction */

void imgc_parse_x_axis(argparser_t *pp, bool_t *xLeft);
  /*  Parses the "-xAxis" command line option, according to
    {imgc_parse_x_axis_HELP}. Sets {*xLeft} to TRUE if "-xAxis left" is
    present, to FALSE if "-xAxis right" is present. If the option is
    not present, leaves {*xLeft} unchanged. */

#define imgc_parse_x_axis_HELP \
  "[ -xAxis { right | r | left | l } ]"

#define imgc_parse_x_axis_HELP_INFO \
  "  -xAxis { right | r | left | l }\n" \
  "    This optional argument specifies the direction of the" \
  " horizontal axis of the user coordinate system in each image.  It also defines the" \
  " default horizontal position of the origin: at the left edge" \
  " of the image if \"right\", or at the right edge" \
  " if \"left\".  The keyword \"-hAxis\" is accepted as synonymous of \"-xAxis\", and" \
  " the directions may be abbreviated as \"r\" and \"l\", respectively."
  
#define imgc_parse_x_axis_pbm_default_INFO \
  "The default is \"-xAxis right\" as in most image processing tools."
  
#define imgc_parse_x_axis_math_default_INFO \
  "The default is \"-xAxis right\" following mathematical tradition."
  
/* Vertical axis direction */

void imgc_parse_y_axis(argparser_t *pp, bool_t *yDown);
  /*  Parses the "-yAxis" command line option, according to
    {imgc_parse_y_axis_HELP}. Sets {*yDown} to TRUE if the "-yAxis
    down" is present, to FALSE if "-yAxis up" is present. If the
    option is not present, leaves {*yDown} unchanged. */

#define imgc_parse_y_axis_HELP \
  "[ -yAxis { down | d | up | u } ]"

#define imgc_parse_y_axis_HELP_INFO \
  "  -yAxis { down | d | up | u }\n" \
  "    This optional argument specifies the direction of" \
  " the vertical axis of the user coordinate system in each image.  It also defines the" \
  " default vertical position of the origin: at the top edge of the" \
  " image if \"down\", or at the bottom edge if \"up\".  The" \
  " keyword \"-vAxis\" is accepted as synonymous of \"-yAxis\", and" \
  " the durections may be abbreviated as \"d\" and \"u\", respectively."
  
#define imgc_parse_y_axis_pbm_default_INFO \
  "The default is \"-yAxis down\" as in most image processing tools." \

#define imgc_parse_y_axis_math_default_INFO \
  "The default is \"-yAxis up\" following mathematical tradition." \
  
/* Input user unit. */

void imgc_parse_input_unit(argparser_t *pp, double *iUnit);
  /* Parses the command line argument "-iUnit" according to the syntax
    specs {imgc_parse_input_unit_HELP}, and returns it in {*iUnit}. If
    \"-iUnit\" is not present, sets it to 1. */

#define imgc_parse_input_unit_HELP \
  "[ -iUnit {IUNIT} ]"

#define imgc_parse_input_unit_HELP_INFO \
  "  -iUnit {IUNIT}\n" \
  "    This optional argument specifies the size of " \
  " the user coordinate unit in pixels.  If this argument is not specified, the program" \
  " assumes \"-iUnit 1\" (which means that user units are in pixels)."
  
#define imgc_input_unit_affects_org_INFO \
  "The option \"-iUnit\" affects the interpretation of \"-iOrg\" as" \
  " well as other coordinate parameters for input images.  Thus, given \"-iUnit 0.5\", the" \
  " option \"-iOrg 24 30\" will set the user input image origin at 12 pixel columns" \
  " and 15 pixel rows away from the default origin."

/* Output user unit. */

void imgc_parse_output_unit(argparser_t *pp, double *oUnit);
  /* Parses the command line argument "-oUnit" according to the syntax
    specs {imgc_parse_output_unit_HELP}, and returns it in {*oUnit}. If
    \"-oUnit\" is not present, sets it to 1. */

#define imgc_parse_output_unit_HELP \
  "[ -oUnit {OUNIT} ]"

#define imgc_parse_output_unit_HELP_INFO \
  "  -oUnit {OUNIT}\n" \
  "    This optional argument specifies the size of " \
  " the user coordinate unit in pixels.  If this argument is not specified, the program" \
  " assumes \"-oUnit 1\" (which means that output user units are in pixels)."
  
#define imgc_output_unit_org_INFO \
  "The option \"-oUnit\" affects the interpretation of \"-oOrg\".  Thus, given \"-oUnit 0.5\", the" \
  " option \"-oOrg 24 30\" will set the user output image origin at 12 pixel columns" \
  " and 15 pixel rows away from the default origin."
  
#define imgc_output_unit_affects_output_size_INFO \
  "The option \"-oUnit\" affects the interpretation" \
  " of \"-oSize\".  Thus, given \"-oUnit 2.0\", the" \
  " option \"-oSize 30 40\" will set the actual user output" \
  " image size to 60 columns and 80 rows."

#define imgc_output_size_defaults_to_input_INFO \
  "If \"-oUnit\" is omitted, assumes that the output image" \
  " size (in output user coords units) is the same as the input" \
  " image size (in input user coords units)."
  
#define imgc_input_unit_affects_default_output_size_INFO \
  "The parameter \"-iUnit\" affects the default output" \
  " image size.  Suppose the input image has 200 by 300 pixels.  Given the" \
  " argument \"-iUnit 0.5\", the default output image size would be 400 by 600" \
  " user units.  Given also \"-oUnit 3.0\", the default output" \
  " image size would then be 1200 by 1800 pixels."

/* Input image origin */

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

#define imgc_parse_input_center_org_HELP_INFO \
  "  -iCenter\n" \
  "  -iOrg {CX_IN} {CY_IN}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the user coordinate system for input" \
  " images.  The \"-iCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-iOrg\" option sets" \
  " the origin displaced {CX_IN} and {CY_IN} /user/ units from the" \
  " default origin, in the directions specified by \"-xAxis\" and \"-yAxis\".  If" \
  " this argument is not specified, the program" \
  " assumes \"-iOrg 0 0\" (which means the default origin)."

/* Output image origin */

void imgc_parse_output_center_org(argparser_t *pp, bool_t *oCenter, r2_t *oOrg);
  /* Analogous to {imgc_parse_center_org}, but applying to output
    images, according to {imgc_parse_output_center_org_HELP_INFO}. 
    Looks for the keywords "-oCenter" and "-oOrg". */

#define imgc_parse_output_center_org_HELP \
  "[ -oCenter | -oOrg {CX_OUT} {CY_OUT} ]"

#define imgc_output_origin_INFO \
  "The default coordinate system origin for output images" \
  " can be overriden by the \"-oOrg\" or \"-oCenter\" argument."

#define imgc_parse_output_center_org_HELP_INFO \
  "  -oCenter\n" \
  "  -oOrg {CX_OUT} {CY_OUT}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the user coordinate system for output" \
  " images.  The \"-oCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-oOrg\" option sets" \
  " the origin displaced by {CX_OUT} and {CY_OUT} from the" \
  " default origin, in the directions specificed by the \"-xAxis\" and \"-yAxis\" options.  If" \
  " this argument is not specified, the program" \
  " assumes \"-oOrg 0 0\" (which means the default origin)."

/* Output image size: */

void imgc_parse_output_size(argparser_t *pp, int *oCols, int *oRows, int max_size);
  /* Parses the output image size option "-oSize", according to {imgc_parse_output_size_HELP}. If
    "-oSize" is present, sets {*oCols} and {*oRows} to the next two
    arguments. Otherwise leaves {*oCols} and {*oRows} unchanged.
    Fails with error if either dimension exceeds {max_size}. */
    
#define imgc_parse_output_size_HELP \
  "[ -oSize {NX} {NY} ]"

#define imgc_parse_output_size_HELP_INFO \
  "  -oSize {NX} {NY}\n" \
  "    Specifies the horizontal and vertical size of the output image, in output /user/ units."

/* Projective maps between images: */

#define imgc_proj_map_INFO \
  "Let's represent a generic point {p} of the plane as" \
  " a three-element homogeneous coordinate row vector [{w},{x},{y}], such that" \
  " ({x/w},{y/w}) are the Cartesian coordinats of {p}.  With these conventions," \
  " the vector-matrix product {p·M} (in that order) will" \
  " give the three homogeneous coordinates [{w'},{x'},{y'}]" \
  " of {M(p)}.  If the computed weight coordinate {w'} is zero or negative, the point" \
  " {M(p)} is invalid, and the image value at {p} is taken to be undefined."

#endif
