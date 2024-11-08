/* User coordinate system for images. */
/* Last edited on 2024-10-31 11:47:30 by stolfi */

#ifndef image_coords_H
#define image_coords_H

#include <bool.h>
#include <r2.h>
#include <hr2.h>
#include <argparser.h>
#include <argparser_geo.h>

/* COORDINATE AXES */

#define imgc_user_axes_intro_INFO \
  "In user-visible contexts such as the command" \
  " line options, points in the domain of an image" \
  " are expressed in terms of a /user coordinate system/" \
  " specific to each image, input or output.  In every user system," \
  " the first (X) axis is horizontal and the second (Y) is vertical.  The" \
  " directions of the axes are determined by" \
  " the \"-xAxis\" and \"-yAxis\" command line arguments."

#define imgc_pixel_axes_intro_INFO \
  "Internally, positions on each image are identified by" \
  " a /pixel coordinate system/" \
  " in which the X axis points right, the Y axis points" \
  " down, the origin is at the upper left corner, and" \
  " the unit is the pixel spacing."
  
#define imgc_pixel_centers_intro_INFO \
  "In a pixel coordinate system, pixel *corners* have integer" \
  " coordinates, and pixel *centers* have half-integer coordinates."
  
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
  " horizontal axis of the user coordinate system in each" \
  " image.  The keyword \"-hAxis\" is accepted as synonymous of \"-xAxis\", and" \
  " the directions may be abbreviated as \"r\" and \"l\", respectively." \
  "  " xaxis_default
  
#define imgc_parse_x_axis_INFO_OPTS_default_pbm \
  "The default is \"-xAxis right\" as in most image processing tools."
  
#define imgc_parse_x_axis_INFO_OPTS_default_math \
  "The default is \"-xAxis right\" following mathematical tradition."

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
  " the vertical axis of the user coordinate system in each image.  The" \
  " keyword \"-vAxis\" is accepted as synonymous of \"-yAxis\", and" \
  " the directions may be abbreviated as \"d\" and \"u\", respectively." \
  "  " yaxis_default
  
#define imgc_parse_y_axis_INFO_OPTS_default_pbm \
  "The default is \"-yAxis down\" as in most image processing tools." \

#define imgc_parse_y_axis_INFO_OPTS_default_math \
  "The default is \"-yAxis up\" following mathematical tradition." \
  
/* USER UNIT */

void imgc_parse_unit(argparser_t *pp, char *unit_key, double *unit);
  /* Parses the command line argument "{unit_key}" according to the syntax
    specs {imgc_parse_unit_HELP(unit_key,unit_tag)}, and returns it in {*unit}. If
    the argument "{unit_key}" is not given on the command line, sets {*unit} to 1. */

#define imgc_parse_unit_HELP(unit_key,unit_tag) \
  "[ " unit_key " {unit" unit_tag "} ]"

#define imgc_user_unit_intro_INFO(unit_key,which_images) \
  "The unit of length of the user coordinate system for " which_images " can" \
  " be specified by the \"" unit_key "\" argument."  

#define imgc_parse_unit_INFO_OPTS(unit_key,unit_tag,which_images) \
  "  " unit_key " {unit" unit_tag "}\n" \
  "    This optional argument specifies the size in pixels of" \
  " the user coordinate unit for " which_images ".  If this" \
  " argument is not specified, the program" \
  " assumes \"" unit_key " 1\" (that is, user units are pixels)."
  
#define imgc_unit_affects_org_INFO_OPTS(unit_key,org_key,which_images) \
  "The option \"" unit_key "\" affects the interpretation of \"" org_key "\" as" \
  " well as other coordinate parameters for " which_images ".  Thus," \
  " given \"" unit_key " 0.5\", the" \
  " option \"" org_key " 24 30\" will set the user" \
  " origin for " which_images " at 12 pixel columns" \
  " and 15 pixel rows away from the \"low\" corner of the image."
  
#define imgc_unit_affects_size_INFO_OPTS(unit_key,size_key,which_images) \
  "The option \"" unit_key "\" affects the interpretation" \
  " of \"" size_key "\".  Thus, given \"" unit_key " 2.0\", the" \
  " option \"" size_key " 30 40\" will set the actual" \
  " size of " which_images ", in pixels, to 60 columns and 80 rows."
  
#define imgc_unit_affects_default_size_INFO_OPTS(ref_unit_key,size_key,ref_images,out_unit_key,which_images) \
  "The parameter \"" ref_unit_key "\" affects the default" \
  " size for " which_images ".  Suppose the size" \
  " of " ref_images " is 200 by 300 pixels.  Given the" \
  " argument \"" ref_unit_key " 0.5\", the default" \
  " size for " which_images " would be 400 by 600" \
  " /user units/.  Given also \"" out_unit_key " 3.0\", the default" \
  " size for " which_images ", in pixels, would then be 1200 by 1800."

/* IMAGE ORIGIN */

void imgc_parse_center_org(argparser_t *pp, char *ctr_key, bool_t *center, char *org_key, r2_t *org);
  /* Parses the mutually exclusive arguments "{ctr_key}" and "{org_key}",
    that define the coordinate origin for certain image(s),
    according to the syntax specs {imgc_parse_center_org_HELP(ctr_key,org_key,...)}. If
    "{ctr_key}" is preent, sets {*center} to TRUE and leaves {*org}
    unchanged. If "{org_key}" is present, sets {*center} to FALSE, and sets
    {*org} to the given point.  NOTE: if neither is present, leaves {*center}
    and {*org} unchanged.  Therefore, the caller must initialize
    those variables with proper defaults. */

#define imgc_parse_center_org_HELP(ctr_key,org_key,org_tag) \
  "[ " ctr_key " | " org_key " {CX" org_tag "} {CY" org_tag "} ]"

#define imgc_user_origin_intro_INFO(ctr_key,org_key,which_images) \
  "The user coordinate system origin for " which_images " can" \
  " be specified by the \"" ctr_key "\" or \"" org_key "\" argument."  

#define imgc_parse_center_org_INFO_OPTS(ctr_key,org_key,org_tag,org_default) \
  "  " ctr_key "\n" \
  "  " org_key " {CX" org_tag "} {CY" org_tag "}\n" \
  "    These mutually exclusive optional arguments specify" \
  " the origin of the user coordinate system for input" \
  " images.  The \"" ctr_key "\" option sets the origin" \
  " at the center of the image domain.  The \"" org_key "\" option sets" \
  " the origin displaced {CX" org_tag "} and {CY" org_tag "} /user/ units" \
  " in the directions of the X and Y axes, from the" \
  " corner of the image opposite to those directions." \
  "  " org_default
  
#define imgc_parse_center_org_INFO_OPTS_default_zero(org_key) \
  "If neither option is specified, the program" \
  " assumes \"" org_key " 0 0\", which means the" \
  " corner of the image opposite to the directions of the X and Y axes."
  
#define imgc_parse_center_org_INFO_OPTS_default_center(ctr_key) \
  "If neither option is specified, the program" \
  " assumes \"" ctr_key "\" (the center of the image)."

/* IMAGE SIZE */

void imgc_parse_size(argparser_t *pp, char *size_key, double *cols, double *rows);
  /* Parses the image size option "{size_key}", according to
    {imgc_parse_size_HELP(size_key,...)}. Namely, if "{size_key}" is present, sets
    {*cols} and {*rows} to the next two arguments. Otherwise leaves
    {*cols} and {*rows} unchanged. */
    
#define imgc_parse_size_HELP(size_key,size_tag) \
  "[ " size_key " {NX" size_tag "} {NY" size_tag "} ]"

#define imgc_parse_size_INFO_OPTS(size_key,size_tag,which_images,size_default) \
  "  " size_key " {NX" size_tag "} {NY" size_tag "}\n" \
  "    Specifies the horizontal and vertical size of " which_images ", in" \
  " /user units/.  These values may be" \
  " fractional. The actual image size in pixels will be rounded up" \
  " to the next integer." \
  "  " size_default

#define imgc_parse_size_INFO_OPTS_default_input(size_key,ref_images,which_images) \
  "If \"" size_key "\" is omitted, the program assumes that the" \
  " size (in /user units/) of " which_images " is the same as the" \
  " size (in /user units/) of " ref_images "."

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

/* PROJECTIVE MAPS */

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
    
    If {center} is true, the US origin is located at the center of the
    image's domain, and the {org} parameter is ignored. If {center} is
    false, the US origin is located, along each axis, at {org.c[ax]}
    USER units from the edge of the image's domain that has lowest USER
    coordinate on that axis. That is, from the left edge of the image if
    {xRev} is false, or from the right edge if {xRev} is true; and from
    the top edge of the image if {yRev} is false, or from the bottom
    edge {yRev} is true. */

/* COMMAND LINE ARGS FOT PROGS WITH BOTH INPUT IMAGES AND OUTPUT IMAGES */

/* The following are standard texts to include in the {PROG_HELP} and {PROG_INFO} strings in 
  programs with let the user define user coordinate systems on both input and output, with 
  variable units. */

#define imgc_input_output_coords_HELP \
  "    " imgc_parse_x_axis_HELP " \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
  "    " imgc_parse_unit_HELP("-iUnit","In") " \\\n" \
  "    " imgc_parse_center_org_HELP("-iCenter","-iOrg","In") " \\\n" \
  "    " argparser_proj_map_HELP " \\\n" \
  "    " imgc_parse_unit_HELP("-oUnit","Out") " \\\n" \
  "    " imgc_parse_center_org_HELP("-oCenter","-oOrg","Out") " \\\n" \
  "    " imgc_parse_size_HELP("-oSize","")

#define imgc_input_output_coords_intro_INFO(iimages,oimages) \
  "  " imgc_user_axes_intro_INFO "" \
  "  " imgc_user_unit_intro_INFO("-iUnit",iimages) "" \
  "  " imgc_user_origin_intro_INFO("-iCenter","-iOrg",iimages) "\n" \
  "\n" \
  "  " imgc_user_unit_intro_INFO("-oUnit",oimages) "" \
  "  " imgc_user_origin_intro_INFO("-oCenter","-oOrg",oimages) "\n" \
  "\n" \
  "  " imgc_pixel_axes_intro_INFO "" \
  "  " imgc_pixel_centers_intro_INFO

#define imgc_parse_input_output_coords_INFO_OPTS(xaxis_def,yaxis_def,iimages,iorg_def,oimages,oorg_def,osize_def) \
  imgc_parse_x_axis_INFO_OPTS(xaxis_def) \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments, including" \
  " {CXIn} and {CXOut}." "\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(yaxis_def) "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments, including" \
  " {CYIn} and {CYOut}." "\n" \
  "\n" \
  imgc_parse_unit_INFO_OPTS("-iUnit","In",iimages) "\n" \
  "\n" \
  imgc_parse_center_org_INFO_OPTS("-iCenter","-iOrg","In",iorg_def) "\n" \
  "\n" \
  "    " imgc_unit_affects_org_INFO_OPTS("-iUnit","-iOrg",iimages) "\n" \
  "\n" \
  imgc_parse_unit_INFO_OPTS("-oUnit","Out",oimages) "\n" \
  "\n" \
  imgc_parse_size_INFO_OPTS("-oSize","",oimages,osize_def) "\n" \
  "\n" \
  "    " imgc_unit_affects_size_INFO_OPTS("-oUnit","-oSize",oimages) "\n" \
  "\n" \
  "    " imgc_unit_affects_default_size_INFO_OPTS("-iUnit","-oSize",iimages,"-oUnit",oimages) "\n" \
  "\n" \
  imgc_parse_center_org_INFO_OPTS("-oCenter","-oOrg","Out",oorg_def) "\n" \
  "\n" \
  "    " imgc_unit_affects_org_INFO_OPTS("-oUnit","-oOrg",oimages)

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
    {iUnit}", or to 1 if that option is not specified. 
    
    The options "-iCenter" and "-iOrg" are mutually exclusive.
    If "-iCenter" is specified, the parameter {*iCenter} is set to {TRUE}, and
    {*iOrg} is left unchanged. If "-iOrg {CXIn} {CYIn}" if present, then 
    {*iCenter} is set to   {FALSE}, and {*iOrg} is set to {(CXIn,CYIn)} (which are
    in the /user input/ units). If both "-iCenter" and "-iOrg"
    are omitted, leaves {*iCenter} and {*iOrg} unchanged.
    
    The parameters {*oUnit}, {*oCenter}, and {*oOrg} are set in the same
    way from the options "-oUinit {oUnit}", "-oCenter", and "-oOrg
    {CXOut} {CYOut}".  Note that the later are /user output/ units.
    
    The parameters {*oCols} and {*oRows} are the output image dimensions
    in /user/ coordinates. They are set to the values specified by the
    option "-oSize {NX} {NY}" if present; otherwise they are left unchanged. */

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

#endif
