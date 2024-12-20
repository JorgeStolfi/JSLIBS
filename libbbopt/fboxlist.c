/* See fboxlist.h */
/* Last edited on 2024-12-05 10:21:57 by stolfi */

#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>

#include <fboxlist.h>
#include <fbox.h>

FBoxList fbox_cons(FBox *b, FBoxList L)
  { void *v = notnull(malloc(sizeof(FBoxListNode)), "no mem for FBoxListNode");
    FBoxListNode *c = (FBoxListNode *)v;
    c->b = b;
    c->next = L;
    return c;
  }
