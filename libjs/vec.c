/* See vec.h */
/* Last edited on 2019-01-11 19:55:36 by jstolfi */

#include <vec.h>
#include <stdlib.h>
#include <stdint.h>
#include <affirm.h>
#include <assert.h>

void *vec_alloc(vec_size_t ne, size_t esz)
  { demand(ne <= vec_MAX_SIZE, "too many elements");
    void *e = (ne == 0 ? NULL : malloc(ne*esz));
    affirm((ne == 0) || (e != NULL), "out of mem");
    return e;
  }

void vec_expand(vec_size_t *nep, void **ep, vec_index_t index, size_t esz)
  { demand(index <= vec_MAX_INDEX, "index too large");
    if (index >= (*nep))
      { vec_size_t ne = (*nep);
        assert(ne <= vec_MAX_SIZE);
        if (index + 1 > vec_MAX_SIZE - ne) 
          { ne = vec_MAX_SIZE; } 
        else
          { ne = ne + index + 1; }
        if ((*nep) == 0) { affirm((*ep) == NULL, "bad elem pointer"); } 
        (*ep) = realloc((*ep), ne*esz);
        affirm((*ep) != NULL, "out of mem");
        (*nep) = ne;
      }
  }

void vec_trim(vec_size_t *nep, void **ep, vec_size_t ne, size_t esz)
  { if (ne != (*nep))
      { if (ne == 0) 
          { free((*ep)); (*ep) = NULL; }
        else
          { (*ep) = realloc((*ep), ne*esz);
            affirm((*ep) != NULL, "out of mem");
          }
        (*nep) = ne;
      }
  }

/* SOME USEFUL TYPED VECTORS */

vec_typeimpl(int_vec_t,    int_vec,    int);
vec_typeimpl(uint_vec_t,   uint_vec,   unsigned int);
vec_typeimpl(char_vec_t,   char_vec,   char);
vec_typeimpl(bool_vec_t,   bool_vec,   bool_t);
vec_typeimpl(float_vec_t,  float_vec,  float);
vec_typeimpl(double_vec_t, double_vec, double);
vec_typeimpl(ref_vec_t,    ref_vec,    void*);
vec_typeimpl(string_vec_t, string_vec, char*);


vec_typeimpl(int8_vec_t,    int8_vec,   int8_t);
vec_typeimpl(int16_vec_t,   int16_vec,  int16_t);
vec_typeimpl(int32_vec_t,   int32_vec,  int32_t);
vec_typeimpl(int64_vec_t,   int64_vec,  int64_t);
vec_typeimpl(uint8_vec_t,   uint8_vec,  uint8_t);
vec_typeimpl(uint16_vec_t,  uint16_vec, uint16_t);
vec_typeimpl(uint32_vec_t,  uint32_vec, uint32_t);
vec_typeimpl(uint64_vec_t,  uint64_vec, uint64_t);
