/* Implementation of ift_queue.h */
/* Last edited on 2003-02-04 22:08:37 by stolfi */

#include "ift_queue.h"
#include <stdlib.h>
#include <stdio.h>


#define NUM_BUCKETS (MAX_FINITE_ARC_COST+1)

#define NB NUM_BUCKETS

typedef struct
  { PixelNode *bucket[NB];    /* Node buckets indexed by (finite) cost mod NUM_BUCKETS. */
    PathCost loC, hiC;        /* Cost range of buckets in use is `[loC..hiC]'. */
    int finentries;           /* Number of finite-cost entries. */
    PixelNode *infinite;      /* Bucket for infinite-cost nodes. */
    ift_tie_breaking tbreak;  /* Queue discipline within each bucket. */
  } ift_queue_impl;
  /*
    At any moment, the finite costs of all nodes in the queue are contained
    in the interval `[loC .. hiC]'. All nodes with the same cost
    `C' are stored in a circular doubly linked list that begins with
    node `*(bucket[C % N])' where `N = NUM_BUCKETS', and is connected through
    the `next' and `prev' fields in each node.  All nodes with infinite
    cost are similarly stored in the list that begins with node `*infinite'. 
    
    This representation assumes that the maximum difference between
    the finite costs of any two nodes in the queue is strictly less
    than `N' --- an assumption which should be true in the context of
    the IFT algorithm. */

void ift_bucket_insert(PixelNode **bu, PixelNode *p, ift_tie_breaking tbreak);
  /* 
    Insert pixel `p' into the bucket pointed to by `*bu',
    at the end specified by the `tbreak' parameter
    (FIFO = rear, LIFO = front). */

void ift_bucket_delete(PixelNode **bu, PixelNode *p);
  /* 
    Delete pixel `p' from its current bucket. */

ift_queue * new_ift_queue(ift_tie_breaking tbreak)
  { 
    int i;
    ift_queue_impl *QQ = (ift_queue_impl *)malloc(sizeof(ift_queue_impl));
    if (QQ == NULL) { IFT_ERROR("out of memory"); }
    for (i = 0; i < NB; i++) { QQ->bucket[i] = NULL; }
    QQ->loC = 0;
    QQ->hiC = 0;
    QQ->finentries = 0;
    QQ->infinite = NULL;
    QQ->tbreak = tbreak;
    return QQ;
  }

void free_ift_queue(ift_queue *Q)
  { 
    free(Q);
  }
  
int ift_queue_is_empty(ift_queue *Q)
  { 
    ift_queue_impl *QQ = (ift_queue_impl *)Q;
    return ((QQ->finentries == 0) & (QQ->infinite == NULL));
  }

int ift_queue_is_enqueued(PixelNode *p)
  { 
    if ((p->next == NULL) != (p->prev == NULL)) { IFT_ERROR("inconsistent links"); } 
    return (p->next != NULL);
  }

void ift_queue_insert(ift_queue *Q, PixelNode *p)
  {
    ift_queue_impl *QQ = (ift_queue_impl *)Q;
    PixelNode **bu;
    /* locate the right bucket for `p': */
    if (p->C == INFINITE_PATH_COST) 
      { bu = &(QQ->infinite); }
    else
      { if (p->C > MAX_FINITE_PATH_COST) { IFT_ERROR("cost too big"); }
        if (fmod(p->C, 1.0) != 0.0) { IFT_ERROR("non-integer cost"); }
        bu = &(QQ->bucket[(int)fmod(p->C, (double)NB)]);
        /* Make sure that `p->C' lies in `[QQ->loC .. QQ->hiC-1]': */
        if (QQ->finentries == 0)
          { QQ->loC = p->C;
            QQ->hiC = p->C;
          }
        else if (p->C < QQ->loC)
          { int nb = (int)(QQ->hiC - p->C) + 1;
            int hi = (int)fmod(QQ->hiC, (double)NB);
            /* Reduce upper end of interval, if needed, to make room for `p': */
            /* Actually this shouldn't be needed if we are lazy enough. */
            while ((nb > NB) && (QQ->bucket[hi] == NULL)) 
              { nb--; hi = (hi + NB - 1) % NB; }
            if (nb > NB) { IFT_ERROR("cost range too big"); }
            QQ->loC = p->C;
            QQ->hiC = QQ->loC + (double)(nb - 1);
          }
        else if (p->C > QQ->hiC)
          { int nb = (int)(p->C - QQ->loC) + 1;
            int lo = (int)fmod(QQ->loC, (double)NB);
            /* Increase lower end of interval, if needed, to make room for `p': */
            while ((nb > NB) && (QQ->bucket[lo] == NULL)) 
              { nb--; lo = (lo + 1) % NB; }
            if (nb > NB) { IFT_ERROR("cost range too big"); }
            QQ->hiC = p->C;
            QQ->loC = QQ->hiC - (double)(nb - 1);
          }
        QQ->finentries++;
      }

    /* Now insert `p' in its bucket: */
    ift_bucket_insert(bu, p, QQ->tbreak);
  }
  
void ift_queue_delete(ift_queue *Q, PixelNode *p)
  {
    ift_queue_impl *QQ = (ift_queue_impl *)Q;
    PixelNode **bu;
    if (p->C == INFINITE_PATH_COST) 
      { bu = &(QQ->infinite); }
    else
      { if (QQ->finentries == 0) { IFT_ERROR("bucket not found"); }
        /* Just remove the node from its bucket list: */
        if ((p->C < QQ->loC) || (p->C > QQ->hiC)) 
          { IFT_ERROR("inconsistent queue costs"); }
        bu = &(QQ->bucket[(int)fmod(p->C, (double)NB)]);
        QQ->finentries--;
      }
    ift_bucket_delete(bu, p);
    /* There is no need to fix `QQ->loC' and `QQ->hiC' now. */
  }
    
PixelNode *ift_queue_delete_min(ift_queue *Q)
  { 
    ift_queue_impl *QQ = (ift_queue_impl *)Q;
    PixelNode **bu;
    PixelNode *p;
    if (QQ->finentries == 0) 
      { bu = &(QQ->infinite); }
    else 
      { /* Find first non-empty bucket: */
        int ip = (int)fmod(QQ->loC, (double)NB);
        while (QQ->bucket[ip] == NULL) { ip = (ip + 1) % NB; }
        /* Free the empty ones: */
        QQ->loC = (QQ->bucket[ip])->C;
        if (ip != (int)fmod(QQ->loC, (double)NB)) 
          { IFT_ERROR("wrong bucket"); }
        bu = &(QQ->bucket[ip]);
        QQ->finentries--;
      }
    p = (*bu);
    ift_bucket_delete(bu, p);
    return p;
  }
    
/* BUCKET OPERATIONS */

void ift_bucket_insert(PixelNode **bu, PixelNode *p, ift_tie_breaking tbreak)
  {
    /* Consistency checks: */
    if ((p->next != NULL) && (p->prev != NULL)) 
      { IFT_ERROR("pixel already enqueued"); }
    else if ((p->next == NULL) != (p->prev == NULL)) 
      { IFT_ERROR("inconsistent null links"); }
    else if ((p->next == p) != (p->prev == p)) 
      { IFT_ERROR("inconsistent self links"); }
      
    /* OK, do it: */
    if ((*bu) == NULL)
      { (*bu) = p;
        p->next = p;
        p->prev = p;
      }
    else
      { /* Insert at rear end of list: */
        PixelNode* front = (*bu);
        p->next = front;
        p->prev = front->prev;
        (p->next)->prev = p;
        (p->prev)->next = p;
        /* If LIFO, move `p' to the front: */
        if (tbreak == ift_tbreak_LIFO) { (*bu) = p; }
     }
  }

void ift_bucket_delete(PixelNode **bu, PixelNode *p)
  {
    /* Consistency checks: */
    if ((*bu) == NULL) 
      { IFT_ERROR("bucket is empty"); }
    else if ((p->next == NULL) && (p->prev == NULL)) 
      { IFT_ERROR("pixel not in queue"); }
    else if ((p->next == NULL) != (p->prev == NULL)) 
      { IFT_ERROR("inconsistent null links"); }
    else if ((p->next == p) != (p->prev == p)) 
      { IFT_ERROR("inconsistent self links"); }
      
    /* OK, do it: */
    if ((p->next == p) && (p->prev == p)) 
      { (*bu) = NULL; }
    else 
      { if ((*bu) == p) { (*bu) = p->next; }
        (p->next)->prev = p->prev;
        (p->prev)->next = p->next;
      }
    p->next = NULL; 
    p->prev = NULL;
  }

