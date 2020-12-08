/* ift_queue.h - Priority queue for IFT algorithm */
/* Last edited on 2001-09-30 04:41:58 by stolfi */

#ifndef ift_queue_H
#define ift_queue_H

#include "ift.h"

typedef void ift_queue;

ift_queue *new_ift_queue(ift_tie_breaking tbreak);
  /*
    A new queue with the specified tie-breaking discipline. */

void free_ift_queue(ift_queue *Q);
  /*
    Releases the storage used by `Q'. */

int ift_queue_is_empty(ift_queue *Q);
  /*
    Returns 1 when `Q' is empty, 0 otherwise. */

int ift_queue_is_enqueued(PixelNode *p);
  /*
    Returns 1 when `p' is inserted in some queue, 0 otherwise. */

void ift_queue_insert(ift_queue *Q, PixelNode *p);
  /*
    Inserts a node `p' in the queue, with its current cost. */
    
void ift_queue_delete(ift_queue *Q, PixelNode *p);
  /*
    Removes the node `p' from the queue. */
    
PixelNode *ift_queue_delete_min(ift_queue *Q);
  /*
    Removes from the queue a node of minimum cost, 
    and returns it. Ties are broken according to the
    time when nodes were last inserted in the queue,
    and the queue's `tbreak' attribute. */
    
#endif
