/* Parsing a frame name pattern from command line args. */
/* Last edited on 2024-12-05 10:30:16 by stolfi */

#ifndef frame_name_pattern_arg_parse_H
#define frame_name_pattern_arg_parse_H

#include <stdio.h>

#include <argparser.h>

char *frame_name_pattern_arg_parse(argparser_t *pp, const char *keyword);
  /* Looks for a command line argument equal to {keyword}, and
    parses the next argument as a file name pattern for the
    video frame images.  See {frame_name_pattern_arg_INFO}. */

#define frame_name_pattern_arg_parse_INFO \
  "The pattern must contain exactly one embedded '%' formatting" \
  " spec with a 'd' conversion code, as" \
  " in \"frame_%06d.pgm\".  Between the" \
  " two there may be any of the format modifiers accepted by the {sprintf}" \
  " function, except '*' (variable width or precision).  Note that" \
  " a pair of consecutive '%'s is not a formatting spec, but one" \
  " literal '%'.  The pattern must include the file name extension. The file" \
  " name is obtained by inserting the frame number in place of" \
  " the '%...d' spec with {sprintf}."

#endif
