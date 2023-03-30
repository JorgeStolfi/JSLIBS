/* See nget.h */
/* Last edited on 2023-03-18 11:16:23 by stolfi */

#include <nget.h>
#include <stdint.h>
#include <fget.h>
#include <affirm.h>
#include <stdio.h>
#include <stdlib.h>

void nget_name_eq(FILE *f, char *name)
  { fget_skip_spaces(f);
    fget_match(f, name); 
    fget_skip_spaces(f);
    fget_match(f, "=");
  }

char nget_char(FILE *f, char *name)
  { nget_name_eq(f, name); 
    return fget_char(f);
  }

bool_t nget_bool(FILE *f, char *name)
  { nget_name_eq(f, name); 
    return fget_bool(f);
  }
  
int32_t nget_int32(FILE *f, char *name)
  { nget_name_eq(f, name);
    return fget_int32(f);
  }
  
int64_t nget_int64(FILE *f, char *name)
  { nget_name_eq(f, name);
    return fget_int64(f);
  }

uint32_t nget_uint32(FILE *f, char *name, int32_t base)
  { nget_name_eq(f, name);
    return fget_uint32(f, base);
  }

uint64_t nget_uint64(FILE *f, char *name, int32_t base)
  { nget_name_eq(f, name);
    return fget_uint64(f, base);
  }

double nget_double(FILE *f, char *name)
  { nget_name_eq(f, name); 
    return fget_double(f);
  }

char *nget_string(FILE *f, char *name)
  { nget_name_eq(f, name); 
    return fget_string(f);
  }
