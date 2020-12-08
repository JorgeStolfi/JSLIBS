/* See fboxlist.h */
/* Last edited on 2005-06-05 15:23:16 by stolfi */

#include <fbox.h>
#include <fboxlist.h>
#include <affirm.h>
#include <stdlib.h>

FBoxList fbox_cons(FBox *b, FBoxList L)
  { void *v = notnull(malloc(sizeof(FBoxListNode)), "no mem for FBoxListNode");
    FBoxListNode *c = (FBoxListNode *)v;
    c->b = b;
    c->next = L;
    return c;
  }
