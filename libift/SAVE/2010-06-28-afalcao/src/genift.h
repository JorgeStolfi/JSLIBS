#ifndef _GENIFT_H_
#define _GENIFT_H_

#include "annimg.h"
#include "queue.h"

/* path-cost functions */

typedef int (*PathCost)(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q); /* lambda >= 0 indicates seed pixels */

/* For watershed transform and superior reconstructions */

int Fsuprec(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q);    

/* For boundary tracking with no orientation */

int Fctrack(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q);
  
/* For approximate EDT and related operators */

int Fedt(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q);     

/* For geodesic dilations (chamfer 5x7)*/

int Fgeo(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q); 

/* For regional minima */

int Fini(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q);  

/* For local superior reconstruction */

int Flsuprec(AnnImg *aimg, Image *lambda, Image *handicap, int p, int q); 

/* algorithms 2 and 3 of the paper */ 

AnnImg *IFTFIFO(Image *img, Image *lambda, Image *handicap, AdjRel *A, PathCost Pcost);
AnnImg *IFTLIFO(Image *img, Image *lambda, Image *handicap, AdjRel *A, PathCost Pcost);

#endif
