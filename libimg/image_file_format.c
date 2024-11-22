/* See {image_file_format.h} */
/* Last edited on 2020-10-31 21:18:47 by jstolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bool.h>
#include <argparser.h>

#include <image_file_format.h>

/* INTERNAL PROTOTYPES */

image_file_format_t image_file_format_convert(char *str, bool_t *okP);
  /* Like image_file_format_from_string, but does not abort on error.
    Instead, returns an arbitrary results and sets {*okP} to {FALSE.
    Otherwise sets {*okP} to {TRUE}. */
    
/* IMPLEMENTATIONS */
  
image_file_format_t image_file_format_from_string(char *str)
  {
    bool_t ok;
    image_file_format_t ffmt = image_file_format_convert(str, &ok);  /* Image format code. */
    if (! ok) 
      { fprintf(stderr, "** unsupported file format \"%s\"\n", str);
        assert(FALSE);
      }
    return ffmt;
  }

image_file_format_t image_file_format_from_name(char *fname)
  {
    char *ext = strrchr(fname, '.');
    if (ext == NULL)
      { fprintf(stderr, "** file name \"%s\" has no extension\n", fname);
        assert(FALSE);
      }
    return image_file_format_from_string(ext);
  }

image_file_format_t image_file_format_arg_parse(argparser_t *pp, const char *keyword)
  { argparser_get_keyword(pp, (char *)keyword);
    char *cfmt = argparser_get_next(pp); /* Format argument ("pgm", "JPEG", etc.) */
    bool_t ok;
    image_file_format_t ffmt = image_file_format_convert(cfmt, &ok);  /* Image format code. */
    if (! ok) { argparser_error(pp, "unsupported file format"); }
    return ffmt;
  }
  
image_file_format_t image_file_format_convert(char *str, bool_t *okP)
  {
    (*okP) = TRUE; /* For now. */
    /* Ignore leading '.' if any: */
    if (str[0] == '.') { str++; }
    if ((strcasecmp(str, "jpeg") == 0) || (strcasecmp(str, "jpg") == 0))
      { return image_file_format_JPG; }
    else if (strcasecmp(str, "png") == 0)
      { return image_file_format_PNG; }
    else if
      ( (strcasecmp(str, "pnm") == 0) ||
        (strcasecmp(str, "pgm") == 0) ||
        (strcasecmp(str, "ppm") == 0) ||
        (strcasecmp(str, "pbm") == 0)
      )
      { return image_file_format_PNM; }
    else if (strcasecmp(str, "fni") == 0)
      { return image_file_format_FNI; }
    else
      { (*okP) = FALSE;
        return image_file_format_PNG; /* Arbitrary. */
      }
  }
