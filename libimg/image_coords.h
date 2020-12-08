/* Coordinate systems for images. */
/* Last edited on 2017-03-11 19:05:28 by stolfilocal */

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
    bool_t center, 
    r2_t *org, 
    int cols, 
    int rows
  );
  /* Returns a projective map {M} that converts point coordinates from
    the {float_image.h} native coordinate system (FS) to the user's
    coordinate system (US).
    
    The FS has the first (X) coordinate increasing with column index,
    the second (Y) coordinate increasing with row index. The image
    domain if the rectangle {[0 _ cols] × [0 _ rows]} where {cols,rows} are
    the image's {sz[1]} and {sz[2]} parameters. The pixel with indices
    {x} and {y} is the unit square with corners {(x,y)} and
    {(x+1,y+1)}.
    
    In the US, the X axis increases with the column index if {xRev} is
    FALSE, and opposite to the column index if {xRev} is TRUE. The Y
    axis increases with row index if {yRev} is FALSE, and opposite to
    row index if {yRev} is true. The US origin is located at the
    center of the image's domain if {center} is TRUE, or at the
    position {org} if {center} is FALSE. In the second case, each
    coordinate of {org} is measured from the edge of the image's
    domain that has lowest coordinate (that is, column 0 if {xRev} is
    FALSE, column {cols-1} if {xRev} is true; row 0 if {yRev} is FALSE,
    row {rows-1} if {yRev} is TRUE).
    
    In any case, all coordinates are measured in pixels. The image's
    domain is assumed to be a rectangle of width {cols} and height
    {rows}. This is equivalent to assuming that a pixel with indices
    {[ix,iy]} is a square, displaced {ix} units to the right and {iy}
    units up from the domain's lower left corner. */

/* COMMAND LINE PARSING */

/* General help about image coordinate systems */

#define imgc_axes_INFO \
  "Points in image domains" \
  " are expressed in terms of a coordinate system" \
  " specific to each image, input or output.  In both systems," \
  " the directions of the axes are determined by" \
  " the \"-xAxis\" and \"-yAxis\" command line arguments.  These" \
  " arguments also affect the default origin of the coordinate" \
  " system.  In the default coordinate systems, the origin lies" \
  " at the low corner of the domain (which depends on the axis directions)."
  
#define imgc_pixel_centers_INFO \
  "The unit of length is the pixel size; in the default coordinate systems," \
  " pixel *corners* have integer coordinates, and pixel *centers* have half-integer coordinates."
  
/* Horizontal axis direction */

void imgc_parse_x_axis(argparser_t *pp, bool_t *xLeft);
  /*  Parses the "-xAxis" command line option, according to
    {imgc_parse_x_axis_HELP}. Sets {*xLeft} to TRUE if "-xAxis left" is
    present, to FALSE if "-xAxis right" is present. If the option is
    not present, leaves {*xLeft} unchanged. */

#define imgc_parse_x_axis_HELP \
  "[ -xAxis { right | left } ]"

#define imgc_parse_x_axis_HELP_INFO \
  "  -xAxis { right | left }\n" \
  "    This optional argument specifies the direction of the" \
  " horizontal axis of image coordinate systems.  It also defines the" \
  " default horizontal position of the origin: at the left edge" \
  " of the image if \"right\", or at the right edge" \
  " if \"left\".  The keyword \"-hAxis\" is accepted as synonymous of \"-xAxis\"."
  
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
  "[ -yAxis { down | up } ]"

#define imgc_parse_y_axis_HELP_INFO \
  "  -yAxis { down | up }\n" \
  "    This optional argument specifies the direction of" \
  " the vertical axis of image coordinate systems.  It also defines the" \
  " default vertical position of the origin: at the top edge of the" \
  " image if \"down\", or at the bottom edge if \"up\".  The" \
  " keyword \"-vAxis\" is accepted as synonymous of \"-yAxis\"."
  
#define imgc_parse_y_axis_pbm_default_INFO \
  "The default is \"-yAxis down\" as in most image processing tools." \

#define imgc_parse_y_axis_math_default_INFO \
  "The default is \"-yAxis up\" following mathematical tradition." \
  
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
  "[ -iCenter | -iOrg {CX_IN} {CY_IN} ]" \

#define imgc_input_origin_INFO \
  "The default coordinate system origin for input images" \
  " can be overriden by the \"-iOrg\" or \"-iCenter\" argument."  

#define imgc_parse_input_center_org_HELP_INFO \
  "  -iCenter\n" \
  "  -iOrg {CX_IN} {CY_IN}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the coordinate system for input" \
  " images.  The \"-iCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-iOrg\" option sets" \
  " the origin at the point {(CX_IN,CY_IN)} relative to the" \
  " default origin.  If this argument is not specified, the program" \
  " assumes \"-iOrg 0 0\" (which means the default origin)."

/* Output image origin */

void imgc_parse_output_center_org(argparser_t *pp, bool_t *oCenter, r2_t *oOrg);
  /* Analogous to {imgc_parse_center_org}, but applying to output
    images, according to {imgc_parse_output_center_org_HELP_INFO}. 
    Looks for the keywords "-oCenter" and "-oOrg". */

#define imgc_parse_output_center_org_HELP \
  "[ -oCenter | -oOrg {CX_OUT} {CY_OUT} ]" \

#define imgc_output_origin_INFO \
  "The default coordinate system origin for output images" \
  " can be overriden by the \"-oOrg\" or \"-oCenter\" argument."

#define imgc_parse_output_center_org_HELP_INFO \
  "  -oCenter\n" \
  "  -oOrg {CX_OUT} {CY_OUT}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the coordinate system for output" \
  " images.  The \"-oCenter\" option sets the origin" \
  " at the center of the image domain.  The \"-oOrg\" option sets" \
  " the origin at the point {(CX_OUT,CY_OUT)} relative to the" \
  " default origin.  If this argument is not specified, the program" \
  " assumes \"-oOrg 0 0\" (which means the default origin).\n"

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
  "    Specifies the horizontal and vertical size of the output image."

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
