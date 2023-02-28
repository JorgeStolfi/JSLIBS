/* See fboxlist.h */
/* Last edited on 2023-02-20 06:38:08 by stolfi */

#define _GNU_SOURCE
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
