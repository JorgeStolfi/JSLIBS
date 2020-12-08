/* Image file formats. */
/* Last edited on 2020-10-30 22:39:54 by jstolfi */

#ifndef image_file_format_H
#define image_file_format_H

#define _GNU_SOURCE

#include <argparser.h>

/* CODES FOR IMAGE FILE FORMATS */

typedef enum
  { image_file_format_PNG,  /* PNG (Portable Network Graphics) format. */
    image_file_format_PNM,  /* PNM (Portable AnyMap) format; includes PBM, PGM, PPM. */
    image_file_format_JPG,  /* JPEG format. */
    image_file_format_FNI   /* J. Stolfi's float image format. */
  } image_file_format_t;
  /* Codes for different image file formats. */

image_file_format_t image_file_format_from_string(char *str);
  /* Maps the string string {str} to a file format code, 
    as described by {image_file_format_arg_INFO} below.  
    Ignores a leading '.' if present.  Aborts with error 
    if {str} is not among the recognized specs. */

image_file_format_t image_file_format_from_name(char *fname);
  /* Infers the file format from the file name, by looking for 
    the extension and mapping it as described by 
    {image_file_format_arg_INFO} below.  Aborts with error 
    if the extension is not among the recognized ones.  */

image_file_format_t image_file_format_arg_parse(argparser_t *pp, const char *keyword);
  /* Looks for a command line argument equal to {keyword} and parses the next
    argument as a file format spec, as described by {image_file_format_arg_INFO} below.
    Fails with an error message if it is not one of the valid specs. */
  
#define image_file_format_arg_INFO \
  "Allowed values" \
  " are \"pnm\" (meaning the classic Portable" \
  " Bitmap format), \"jpg\" or \"jpeg\" (Joint Photographic Experts Group" \
  " format), or \"png\" (Portable Network Graphics format); or" \
  " their uppercase versions.  The \"pnm\" tag allows any of the" \
  " formats PBM (bilevel), PGM (grayscale), and PPM (RGB color).  The" \
  " tags \"pbm\", \"ppm\", and \"pgm\", in upper or lower case, are" \
  " also accepted as equivalent to \"pnm\", and each of them" \
  " specifies 'any of the three formats'."

#endif
