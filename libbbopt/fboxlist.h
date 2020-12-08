/* fboxlist.h -- Linked list of f-boxes */
/* Last edited on 2003-09-21 12:01:55 by stolfi */

#ifndef fboxlist_h
#define fboxlist_h

#include <fbox.h>

typedef struct FBoxListNode /* Node of a list of f-boxes. */
  { FBox *b;
    struct FBoxListNode *next;
  } FBoxListNode;
  
typedef FBoxListNode *FBoxList;

FBoxList fbox_cons(FBox *b, FBoxList L);
  /* Prepends {b} to {L}, allocating a new {FBoxListNode}. */

#endif
