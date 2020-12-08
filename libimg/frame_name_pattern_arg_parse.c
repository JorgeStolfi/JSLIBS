/* See {frame_name_pattern_arg_parse.h} */
/* Last edited on 2017-06-21 22:58:58 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include <argparser.h>

#include <frame_name_pattern_arg_parse.h>

char *frame_name_pattern_arg_parse(argparser_t *pp, const char *keyword)
  { 
    argparser_get_keyword(pp, (char *)keyword);
    char *fpat = argparser_get_next(pp);
    
    /* Check it it has the formatting spec inside: */
    char *p = strchr(fpat, '%');
    char *q = (p == NULL ? NULL : p + strcspn(p, "d"));
    if ((p == NULL) || (q == NULL) || ((*q) == 0))
      { argparser_error(pp, "file name should include a \"%%..d\" format spec"); }
    
    return fpat;
  }
