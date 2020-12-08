/* pnmift.h - options, types and prototypes for pnmift.c */
/* Last edited on 2001-09-08 00:04:43 by stolfi */

/* PIXEL POINTERS */

#define MAXROWS 32767
#define MAXCOLS 32767

typedef long pixel_pos;
  /*
    Indices `(h,v)' of a pixel in an image, packed
    as `v*B + h' where `B = 2^16'. */ 
 
typedef long pixel_delta;
  /* 
    Signed index increment `(dh,dv)' between two pixels, packed
    as `(dv+S)*B + (dh+S)' where `B = 2^16' and `S = 2^15'. */
    

