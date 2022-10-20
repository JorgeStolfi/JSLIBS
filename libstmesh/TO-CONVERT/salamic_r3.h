/* Cartesian 3D geometry */

#ifndef salamic_r3_H
#define salamic_r3_H

/*
#include <glm/glm.h>
#include <stdint.h>
#include <glm/gtc/matrix_transform.h>
*/

typedef struct salamic_stl_r3_t {
  float x, y, z;
} salamic_stl_r3_t; 
/* A point of {R 3}. */

int32_t salamic_r3_hash(const salamic_stl_r3_t *v);
/* Maps {v} to an integer, pseudorandomly. */

//  float salamic_r3_distTo(const salamic_stl_r3_t *a, const salamic_stl_r3_t *b);
//  /* Euclideand distance. */

//  void transform(salamic_stl_r3_t *v, const glm::mat4 *mat);
//  /* Applies the transformation {mat} to {*v}. */

#endif
