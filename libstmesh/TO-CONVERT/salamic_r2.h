/* Cartesian 2D geometry */

#ifndef salamic_r2_H
#define salamic_r2_H

#include <stdint.h>

using namespace std;

typedef struct salamic_r2_t {
  float x, y;
} salamic_r2_t; 
/* A point of {R 2}. */

int32_t salamic_r2_hash (salamic_r2_t *v);
/* Maps a point {v} to an integer, pseudorandomly. */

typedef struct salamic_r2_Segment_t {
    salamic_r2_t v[2];
} salamic_r2_Segment_t;
/* A segment defined by two points of {R 2}. */

int32_t salamic_r2_Segment_hash (salamic_r2_Segment_t *seg);
/* Maps a segment {seg} to an integer, pseudorandomly. */

typedef vector<salamic_r2_t> salamic_r2_Polygon_t;
/* A polygon is a list of vertices, which are points of {R 2}. */

#endif
