/*
  Copyright (C) <2003> <Alexandre Xavier Falcão>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  please see full copyright in COPYING file.
  -------------------------------------------------------------------------
  written by A.X. Falcão <afalcao@ic.unicamp.br>, March 27th 2003

  This program is a collection of functions to create, destroy, and
  manipulate a priority queue.

  A priority queue Q consists of two data structures: a circular
  queue C and a table L that encodes all possible doubly-linked
  lists.
  
  Q requires that the maximum possible increment along the paths be a
  non-negative integer less than the number of buckets in C. An extra
  bucket is created to store infinity costs for the LIFO policy. The
  queue size increases dynamically whenever (maxcost-mincost) >
  (nbuckets-1).
  
  Q->C.first[i] gives the first element that is in bucket i.
  Q->C.last[i]  gives the last  element that is in bucket i.
  Q->C.nbuckets gives the number of buckets in C.
  Q->C.mincost  gives the minimum cost of a node in queue.
  Q->C.maxcost  gives the maximum cost of a node in queue.
  Q->C.tiebreak gives the FIFO or LIFO tie breaking policies
  
  All possible doubly-linked lists are represented in L. Each bucket
  contains a doubly-linked list that is treated as a FIFO.
  
  Q->L.elem[i].next: the next element to i 
  Q->L.elem[i].prev: the previous element to i
  Q->L.elem[i].color: the color of i (WHITE=never inserted, GRAY=inserted,
  BLACK=removed)
  Q->L.nelems: gives the total number of elements that can be
  inserted in Q (It is usually the number of pixels in a given image
  or the number of nodes in a graph)
  Q->L.cost[i]: gives the cost of element i in the graph. 
  
  Insertions and updates are done in O(1).
  Removal may take O(K+1), where K+1 is the number of buckets.     
*/

#include "queue.h"


Queue *CreateQueue(int nbuckets, int nelems, int *cost)
{
  Queue *Q=NULL;

  Q = (Queue *) malloc(1*sizeof(Queue));
  
  if (Q != NULL) {
    Q->C.first = (int *)malloc((nbuckets+1) * sizeof(int));
    Q->C.last  = (int *)malloc((nbuckets+1) * sizeof(int));
    Q->C.nbuckets = nbuckets;
    if ( (Q->C.first != NULL) && (Q->C.last != NULL) ){
      Q->L.elem = (Node *)malloc(nelems*sizeof(Node));
      Q->L.nelems = nelems;
      Q->L.cost   = cost;
      if (Q->L.elem != NULL){
	ResetQueue(Q);
      } else
	Error(MSG1,"CreateQueue");	
    } else
      Error(MSG1,"CreateQueue");
  } else 
    Error(MSG1,"CreateQueue");
  
  return(Q);
}

void ResetQueue(Queue *Q)
{
  int i;

  Q->C.mincost = INT_MAX;
  Q->C.maxcost = INT_MIN;
  SetTieBreak(Q,FIFOBREAK);
  for (i=0; i < Q->C.nbuckets+1; i++)
    Q->C.first[i]=Q->C.last[i]=NIL;
	
  for (i=0; i < Q->L.nelems; i++) {
    Q->L.elem[i].next =  Q->L.elem[i].prev = NIL;
    Q->L.elem[i].color = WHITE;
  }

}

void DestroyQueue(Queue **Q)
{
  Queue *aux;

  aux = *Q;
  if (aux != NULL) {
    if (aux->C.first != NULL) free(aux->C.first);
    if (aux->C.last  != NULL) free(aux->C.last);
    if (aux->L.elem  != NULL) free(aux->L.elem);
    free(aux);
    *Q = NULL;
  }
}

Queue *GrowQueue(Queue **Q, int nbuckets)
{
  Queue *Q1=CreateQueue(nbuckets,(*Q)->L.nelems,(*Q)->L.cost);
  int i,bucket;

  Q1->C.mincost  = (*Q)->C.mincost;
  Q1->C.maxcost  = (*Q)->C.maxcost;
  Q1->C.tiebreak = (*Q)->C.tiebreak;
  for (i=0; i<(*Q)->C.nbuckets; i++) 
    if ((*Q)->C.first[i]!=NIL){
      bucket = (*Q)->L.cost[(*Q)->C.first[i]]%Q1->C.nbuckets;
      Q1->C.first[bucket] = (*Q)->C.first[i];
      Q1->C.last[bucket]  = (*Q)->C.last[i];
    }
  if ((*Q)->C.first[(*Q)->C.nbuckets]!=NIL){
    bucket = Q1->C.nbuckets;
    Q1->C.first[bucket] = (*Q)->C.first[(*Q)->C.nbuckets];
    Q1->C.last[bucket]  = (*Q)->C.last[(*Q)->C.nbuckets];
  }

  for (i=0; i < (*Q)->L.nelems; i++) 
      Q1->L.elem[i]  = (*Q)->L.elem[i];

  DestroyQueue(Q);
  return(Q1);
}


void InsertQueue(Queue **Q, int elem)
{
  int bucket,mincost=(*Q)->C.mincost,maxcost=(*Q)->C.maxcost;

  if ((*Q)->L.cost[elem]==INT_MAX)
    bucket=(*Q)->C.nbuckets;
  else{
    if ((*Q)->L.cost[elem] < mincost) 
      mincost = (*Q)->L.cost[elem];
    if ((*Q)->L.cost[elem] > maxcost) 
      maxcost = (*Q)->L.cost[elem];
    if ((maxcost-mincost) > ((*Q)->C.nbuckets-1)){      
      (*Q) = GrowQueue(Q,2*(maxcost-mincost)+1);
      printf("growing queue\n");
    }
    bucket=(*Q)->L.cost[elem]%(*Q)->C.nbuckets;
    (*Q)->C.mincost = mincost; 
    (*Q)->C.maxcost = maxcost;  
  }
  if ((*Q)->C.first[bucket] == NIL){ 
    (*Q)->C.first[bucket]   = elem;  
    (*Q)->L.elem[elem].prev = NIL;
  }else {
    (*Q)->L.elem[(*Q)->C.last[bucket]].next = elem;
    (*Q)->L.elem[elem].prev = (*Q)->C.last[bucket];
  }
  
  (*Q)->C.last[bucket]     = elem;
  (*Q)->L.elem[elem].next  = NIL;
  (*Q)->L.elem[elem].color = GRAY;
}

int RemoveQueue(Queue *Q)
{
  int elem=NIL, next, prev;
  int last, current=Q->C.mincost%Q->C.nbuckets;

  /** moves to next element **/

  if (Q->C.first[current] == NIL) {
    last = current;
    
    current = (current + 1) % (Q->C.nbuckets);
    
    while ((Q->C.first[current] == NIL) && (current != last)) {
      current = (current + 1) % (Q->C.nbuckets);
    }
    
    if (Q->C.first[current] != NIL){
      Q->C.mincost = Q->L.cost[Q->C.first[current]];
    }else{
      if (Q->C.first[Q->C.nbuckets] != NIL){
	current = Q->C.nbuckets;
	Q->C.mincost = Q->L.cost[Q->C.first[current]];
      }else{
	Error("Queue is empty\n","RemoveQueue");
      }
    }
  }

  if (Q->C.tiebreak == LIFOBREAK) {
    elem = Q->C.last[current];
    prev = Q->L.elem[elem].prev;
    if (prev == NIL) {         /* there was a single element in the list */
      Q->C.last[current] = Q->C.first[current]  = NIL;    
    }
    else {
      Q->C.last[current]   = prev;
      Q->L.elem[prev].next = NIL;
    }
  } else { /* Assume FIFO policy for breaking ties */
    elem = Q->C.first[current];
    next = Q->L.elem[elem].next;
    if (next == NIL) {         /* there was a single element in the list */
      Q->C.first[current] = Q->C.last[current]  = NIL;    
    }
    else {
      Q->C.first[current] = next;
      Q->L.elem[next].prev = NIL;
    }
  }

  Q->L.elem[elem].color = BLACK;

  return elem;
}

void RemoveQueueElem(Queue *Q, int elem)
{
  int prev,next,bucket;

  if (Q->L.cost[elem] == INT_MAX)
    bucket = Q->C.nbuckets;
  else
    bucket = Q->L.cost[elem]%Q->C.nbuckets;

  prev = Q->L.elem[elem].prev;
  next = Q->L.elem[elem].next;
  
  /* if elem is the first element */
  if (Q->C.first[bucket] == elem) {
    Q->C.first[bucket] = next;
    if (next == NIL) /* elem is also the last one */
      Q->C.last[bucket] = NIL;
    else
      Q->L.elem[next].prev = NIL;
  }
  else{   /* elem is in the middle or it is the last */
    Q->L.elem[prev].next = next;
    if (next == NIL) /* if it is the last */
      Q->C.last[bucket] = prev;
    else 
      Q->L.elem[next].prev = prev;
  }

  Q->L.elem[elem].color = BLACK;

}

void UpdateQueue(Queue *Q, int elem, int newcost)
{
  RemoveQueueElem(Q, elem);
  Q->L.cost[elem] = newcost;
  InsertQueue(&Q, elem);
}

int EmptyQueue(Queue *Q)
{
  int last,current=Q->C.mincost%Q->C.nbuckets;

  if (Q->C.first[current] != NIL)
    return 0;
  
  last = current;
  
  current = (current + 1) % (Q->C.nbuckets);
  
  while ((Q->C.first[current] == NIL) && (current != last)) {
    current = (current + 1) % (Q->C.nbuckets); 
  }
  
  if (Q->C.first[current] == NIL){
    if (Q->C.first[Q->C.nbuckets] == NIL){
      return(1);
    }
  }

  return (0);
}








