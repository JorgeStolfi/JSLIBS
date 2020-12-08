/* Polynomial basis for 2D approximation. 
**
** Last edited on 2007-01-04 03:33:52 by stolfi
**
** Copyright © 2003 by Jorge Stolfi and Anamaria Gomide, the
** University of Campinas, Brazil. See the rights and conditions
** notice at the end of this file.
*/

#ifndef FBasPoly2_H
#define FBasPoly2_H

#include <FBas.h>

#include <stdio.h>

#define FBasPoly2_TypeId (FBas_TypeId "P2.")

typedef struct FBasPoly2_Data 
  { FBas_Data pd;  /* Data fields inherited from parent class. */
  } FBasPoly2_Data;

typedef struct FBasPoly2_Mths 
  { FBas_Mths pm;          /* Methods inherited from parent class. */
    FBas_WriteMth *write;  /* Writes the basis to a stream, as a concrete type. */
  } FBasPoly2_Mths;

typedef struct FBasPoly2
  { char *type;          /* Identifier of concrete type. */
    FBasPoly2_Data *d;   /* Internal parameters. */
    FBasPoly2_Mths *m;   /* Methods. */
  } FBasPoly2; 

FBasPoly2 *FBasPoly2_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {FBasPoly2},
    returns {f} cast to that type; otherwise returns NULL. */

FBasPoly2 *FBasPoly2_New(void);
   /* Returns a new polynomial basis object. */

FBasPoly2 *FBasPoly2_Read(FILE *rd);
  /* Reads an {FBasPoly2} from {rd}. It must have been created by
    calling the {write} method of an {FBasPoly2}. */

#endif
