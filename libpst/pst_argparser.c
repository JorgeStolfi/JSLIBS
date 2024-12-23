/* See pst_argparser.h */
/* Last edited on 2024-12-22 11:38:22 by stolfi */ 

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <float_image.h>
#include <r2.h> 
#include <vec.h> 
#include <affirm.h> 
#include <argparser.h> 

#include <pst_argparser.h>
#include <pst_basic.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

bool_t pst_keyword_present(argparser_t *pp, char *key, bool_t next)
  {
    if (next)
      { return argparser_keyword_present_next(pp, key); }
    else
      { return argparser_keyword_present(pp, key); }
  }

char *pst_parse_next_file_name(argparser_t *pp)
  { /* Peek at next argument: */
    char *nx = argparser_next(pp);
    /* If it is missing or parsed, it is not a file name: */
    if (nx == NULL) { return NULL; }
    /* If it is empty, it is not a file name: */
    if (nx[0] == '\000') { return NULL; }
    /* If it looks like a keyword, it is not a file name: */
    if ((nx[0] == '-') && (nx[1] != '\000')) { return NULL; }
    /* If it begins with a funny character, it is not a file name: */
    if
      ( (nx[0] != '@') &&
        (nx[0] != '/') && 
        (nx[0] != '.') &&
        ((nx[0] < 'A') || (nx[0] > 'Z')) &&
        ((nx[0] < 'a') || (nx[0] > 'z')) &&
        ((nx[0] < '0') || (nx[0] > '9'))
      ) 
      { return NULL; }
    /* Grab the next argument and check if it is a valid filename: */
    nx = argparser_get_next(pp);
    return nx;
  }

name_vec_t pst_parse_file_name_list(argparser_t *pp, uint32_t *NNP)
  { uint32_t NNMax = ((NNP == NULL) || ((*NNP) < 0) ? UINT32_MAX : (*NNP));
    uint32_t NN = 0;
    name_vec_t nvec = name_vec_new(0);
    while (NN < NNMax)
      { char *name = pst_parse_next_file_name(pp);
        if (name == NULL) { break; }
        name_vec_expand(&(nvec), (vec_index_t)NN);
        nvec.e[NN] = name;
        NN++;
      }
    if (NNP != NULL)
      { if ((*NNP) >= 0)
          { if (NN != (*NNP)) { argparser_error(pp, "not enough file names"); } }
        else
          { (*NNP) = NN; }
      }
    name_vec_trim(&(nvec), NN);
    return nvec;
  }
