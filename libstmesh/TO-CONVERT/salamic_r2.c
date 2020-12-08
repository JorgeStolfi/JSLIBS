
/* See {salamic_r2.h} */

#include <math.h>
/*
#include <glm/glm.h>
#include <glm/gtc/matrix_transform.h>
*/

#include <salamic_r2.h>
#include <salamic_utils.h>

int salamic_r2_hash (salamic_r2_t *v) {
  int h = 0;
  salamic_utils_hash_combine(&h, salamic_utils_float_hash(v->x));
  salamic_utils_hash_combine(&h, salamic_utils_float_hash(v->y));
  return h;
}

int salamic_r2_Segment_hash (salamic_r2_Segment_t *seg) {
  int h = 0;
  salamic_utils_hash_combine(&h, salamic_r2_hash(&(seg->v[0])));
  salamic_utils_hash_combine(&h, salamic_r2_hash(&(seg->v[1])));
  return h;
}
