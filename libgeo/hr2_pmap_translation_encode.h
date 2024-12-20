#ifndef hr2_pmap_translation_encode_H
#define hr2_pmap_translation_encode_H

/* Last edited on 2024-12-05 10:27:15 by stolfi */

#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>

void hr2_pmap_translation_encode(hr2_pmap_t *M, double y[]);
  /* Assumes that {M} is a translation. Converts it to two parameters
    {y[0..1]}, which are the Cartesian coordinates of the image of 
    the hither origin {(0,0)}. */

void hr2_pmap_translation_decode(double y[], hr2_pmap_t *M);
  /* Stores into {M} a translaton map with displacement vector {(y[0],y[1])}.
    
    The functions are inverses in the sense that {decode(y',M);
    encode(M,y')} will set {y''} to {y'}, and {encode(M,y); decode(y,M)}
    will have no effect on {M}, if {M} is a translation. */

#endif
