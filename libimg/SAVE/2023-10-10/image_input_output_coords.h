/* Coordinate systems for images. */
/* Last edited on 2023-09-25 09:43:27 by stolfi */

#ifndef image_coords_H
#define image_coords_H

#include <bool.h>
#include <r2.h>
#include <r3x3.h>
#include <hr2.h>
#include <argparser.h>
#include <argparser_geo.h>
#include <image_input_coords.h>
#include <image_output_coords.h>

hr2_pmap_t imgc_coord_sys_map
  ( bool_t xRev, 
    bool_t yRev, 
    double unit,
    bool_t center, 
    r2_t *org, 
    int32_t cols, 
    int32_t rows
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

void imgc_compute_output_size_in_pixels
  ( int32_t iCols_pix,
    int32_t iRows_pix,
    double iUnit,
    double oCols_usr,
    double oRows_usr,
    double oUnit,
    int32_t *oCols_pix_P,
    int32_t *oRows_pix_P,
    int32_t max_pix
  );
  /* Computes the output image dimensions {*oCols_pix_P,*oRows_pix} in
    pixels given the output image dimensions {oCols_usr,oRows_usr} in
    multiple of the output user coordinate unit, and the size {oUnit} of
    the latter in pixels. Fails with error if either
    computed pixel dimension exceeds {max_pix}.
    
    If any of {oCols_usr,oRows_usr} are negative or zero, assumes that
    it is the same as the corresponding input image dimension in pixels,
    {iCols_pix} or {iRows_pix}, divided by the pixel size {iUnit} of the
    input user coordinate unit.  Otherwise ignores that input image dimension.
    If both {oCols_usr} and {oRows_usr} are positive, ignores {iUnit} altogether.*/

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

#define imgc_parse_input_output_coords_INFO_OPTS(xaxis_def,yaxis_def,iorg_def,oorg_def,osize_def) \
  imgc_parse_x_axis_INFO_OPTS(xaxis_def) \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments, including" \
  " {CX_IN} and {CX_OUT}." "\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(yaxis_def) "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments, including" \
  " {CY_IN} and {CY_OUT}." "\n" \
  "\n" \
  imgc_parse_input_unit_INFO_OPTS "\n" \
  "\n" \
  imgc_parse_input_center_org_INFO_OPTS(iorg_def) \
  "  " imgc_input_unit_affects_input_org_INFO_OPTS "\n" \
  "\n" \
  imgc_parse_output_unit_INFO_OPTS "\n" \
  "\n" \
  imgc_parse_output_size_INFO_OPTS(osize_def) "\n" \
  "\n" \
  "  " imgc_output_unit_affects_output_size_INFO_OPTS \
  "  " imgc_input_unit_affects_default_output_size_INFO_OPTS "\n" \
  "\n" \
  imgc_parse_output_center_org_INFO_OPTS(oorg_def) \
  "  " imgc_output_unit_affects_output_org_INFO_OPTS

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
  
/* COMMAND LINE OPTION PARSING */

/* Complete input and output coord systems options: */

void imgc_parse_input_output_coords_args
  ( argparser_t *pp,
    bool_t *xLeft,
    bool_t *yUp,
    double *iUnit,
    bool_t *iCenter,
    r2_t *iOrg,
    double *oUnit,
    bool_t *oCenter,
    r2_t *oOrg,
    double *oCols,
    double *oRows
  );
  /* Parses the command line options that define the correspondece between
    user and pixel coordinates for input and output images, as well as the 
    output image size.
    
    The boolean {*xLeft} is set to {TRUE} if the option "-xAxis left"
    is present, reset to {FALSE} if the option "-xAxis right" is present,
    and left unchanged if the option "-xAxis" is missing.
    
    Likewise, {*yUp} is set to {TRUE} "-yAxis up" is present, 
    to {FALSE} if "-yAxis down" is present, and left unchanged
    if the option "-yAxis" is missing.
    
    The parameter {*iUnit} is set to the value speficied by "-iUint
    {IUNIT}", or to 1 if that option is not specified. 
    
    The options "-iCenter" and "-iOrg" are mutually exclusive.
    If "-iCenter" is specified, the parameter {*iCenter} is set to {TRUE}, and
    {*iOrg} is left unchanged. If "-iOrg {CX_IN} {CY_IN}" if present, then 
    {*iCenter} is set to   {FALSE}, and {*iOrg} is set to {(CX_IN,CY_IN)} (which are
    in the /user input/ units). If both "-iCenter" and "-iOrg"
    are omitted, leaves {*iCenter} and {*iOrg} unchanged.
    
    The parameters {*oUnit}, {*oCenter}, and {*oOrg} are set in the same
    way from the options "-oUinit {OUNIT}", "-oCenter", and "-oOrg
    {CX_OUT} {CY_OUT}".  Note that the later are /user output/ units.
    
    The parameters {*oCols} and {*oRows} are the output image dimensions
    in /user/ coordinates. They are set to the values specified by the
    option "-oSize {NX} {NY}" if present; otherwise they are left unchanged. */
  
void imgc_parse_x_axis(argparser_t *pp, bool_t *xLeft);
  /*  Parses the "-xAxis" command line option, according to
    {imgc_parse_x_axis_HELP}. Sets {*xLeft} to TRUE if "-xAxis left" is
    present, to FALSE if "-xAxis right" is present. If the option is
    not present, leaves {*xLeft} unchanged. */

#define imgc_parse_x_axis_HELP \
  "[ -xAxis { right | r | left | l } ]"

#define imgc_parse_x_axis_INFO_OPTS(xaxis_default) \
  "  -xAxis { right | r | left | l }\n" \
  "    This optional argument specifies the direction of the" \
  " horizontal axis of the user coordinate system in each image.  It also defines the" \
  " default horizontal position of the origin: at the left edge" \
  " of the image if \"right\", or at the right edge" \
  " if \"left\".  The keyword \"-hAxis\" is accepted as synonymous of \"-xAxis\", and" \
  " the directions may be abbreviated as \"r\" and \"l\", respectively." \
  "  " xaxis_default
  
#define imgc_parse_x_axis_INFO_OPTS_default_pbm \
  "The default is \"-xAxis right\" as in most image processing tools."
  
#define imgc_parse_x_axis_INFO_OPTS_default_math \
  "The default is \"-xAxis right\" following mathematical tradition."
  
/* Vertical axis direction */

void imgc_parse_y_axis(argparser_t *pp, bool_t *yUp);
  /*  Parses the "-yAxis" command line option, according to
    {imgc_parse_y_axis_HELP}. Sets {*yUp} to TRUE if the "-yAxis
    up" is present, to FALSE if "-yAxis down" is present. If the
    option is not present, leaves {*yUp} unchanged. */

#define imgc_parse_y_axis_HELP \
  "[ -yAxis { down | d | up | u } ]"

#define imgc_parse_y_axis_INFO_OPTS(yaxis_default) \
  "  -yAxis { down | d | up | u }\n" \
  "    This optional argument specifies the direction of" \
  " the vertical axis of the user coordinate system in each image.  It also defines the" \
  " default vertical position of the origin: at the top edge of the" \
  " image if \"down\", or at the bottom edge if \"up\".  The" \
  " keyword \"-vAxis\" is accepted as synonymous of \"-yAxis\", and" \
  " the directions may be abbreviated as \"d\" and \"u\", respectively." \
  "  " yaxis_default
  
#define imgc_parse_y_axis_INFO_OPTS_default_pbm \
  "The default is \"-yAxis down\" as in most image processing tools." \

#define imgc_parse_y_axis_INFO_OPTS_default_math \
  "The default is \"-yAxis up\" following mathematical tradition." \

/* Projective maps between images: */

#define imgc_proj_map_INFO \
  "Let's represent a generic point {p} of the plane as" \
  " a three-element homogeneous coordinate row vector [{w},{x},{y}], such that" \
  " ({x/w},{y/w}) are the Cartesian coordinats of {p}.  With these conventions," \
  " the vector-matrix product {p·M} (in that order) will" \
  " give the three homogeneous coordinates [{w'},{x'},{y'}]" \
  " of {M(p)}.  If the computed weight coordinate {w'} is zero or negative, the point" \
  " {M(p)} is invalid, and the image value at {p} is taken to be undefined."

void imgc_print_matrix(FILE *wr, char *name1, char *name2, char *name3, char *tag, r3x3_t *M);
  /* Prints to {wr} the projective matrix {M} to file {wr}, labeled with 
    "{name1} to {name2} coordinates ({name3}{tag})". */

void imgc_print_pmap(FILE *wr, char *name1, char *name2, char *name3, hr2_pmap_t *M);
  /* Prints to {wr} the projective map {M} (direct and inverse matrices)
    to file {wr}, labeled with "{name1} to {name2} coordinates ({name3}{tag})" where
    {tag} is "" and "^{-1}", respectively. */

#endif
