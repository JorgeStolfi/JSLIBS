/* See {uint16_color_table.h} */
/* Last edited on 2024-12-26 20:08:50 by stolfi */

/* Copied from Jef Poskanzer's libppm3.c - ppm utility library part 3
**
** Copyright (C) 1989, 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <uint16_image.h>

#include <uint16_color_table.h>

#define HASH_SIZE 20023

#define HASH_MASK (((uint32_t)1 << 30) - 1)
  /* A mask of ones, much greater than {HASH_SIZE} but safe against overflow. */

uint32_t uint16_color_table_hash_pixel(uint32_t chns, uint16_t clr[]);
  /* Computes a hash index for the color value {clr[0..chns-1]}. */
  
/* IMPLEMENTATIONS */

uint32_t uint16_color_table_hash_pixel(uint32_t chns, uint16_t clr[])
  { uint64_t sum = 0;
    for (uint32_t c = 0; c < chns; c++)
      { /* Shouldn't overflow, for any reasonable {chns}: */
        uint32_t w = 30013 + 418*c; 
        sum = sum + (uint64_t)(w * (uint32_t)clr[c]);
      }
    uint64_t mask = ((uint64_t)1 << 32) - 1;
    uint32_t hash = (uint32_t)(sum & mask) ^ (uint32_t)(sum >> 32);
    hash = hash % HASH_SIZE;
    return hash;
  }

uint16_color_table_node_t *uint16_color_table_lookup(uint16_color_table_t *cht, uint32_t chns, uint16_t clr[], bool_t add)
  { for (int32_t c = 0; c < chns; c++) { demand(clr[c] <= cht->maxval, "invalid color value {clr}"); }
    uint32_t hash = uint16_color_table_hash_pixel(chns, clr);
    uint16_color_table_node_t *nd = cht->buckets[hash]; /* The node that matches {clr}, or {NULL} if none. */
    uint16_color_table_node_t *nd_prev = NULL; /* The previous node in bucket, or {NULL} if first. */
    uint16_color_table_node_t *nd_prev_prev = NULL; /* Node before {nd_prev}, or {NULL} if none. */
    while (nd != NULL)
      { if (uint16_image_same_color(chns, nd->clr, clr)) { break; }
        nd_prev_prev = nd_prev;
        nd_prev = nd;
        nd = nd->next;
      }
    if ((nd == NULL) && add)
      { /* Append a new node to the list: */
        nd = talloc(1, uint16_color_table_node_t);
        nd->clr = clr;
        nd->count = 0;
        nd->index = -1;
        nd->next = NULL;
        if (nd_prev == NULL) { cht->buckets[hash] = nd; } else { nd_prev->next = nd; }
        cht->NN++;
      }
    if (nd != NULL)
      { /* Increment the count of {nd}, and maybe bubble it up: */
        (nd->count)++;
        if ((nd_prev != NULL) && (nd->count > nd_prev->count))
          { /* Move the node up one notch in the list: */
            if (nd_prev_prev == NULL) { cht->buckets[hash] = nd; } else { nd_prev_prev->next = nd; }
            nd_prev->next = nd->next;
            nd->next = nd_prev;
          }
      }
    return nd;
  }

uint16_color_table_t *uint16_color_table_build
  ( uint16_t **samples,
    uint32_t chns, 
    uint32_t cols, 
    uint32_t rows, 
    uint16_t maxval,
    uint32_t NN_max
  )
  {
    demand(chns > 0, "invalid channel count {chns}");
    demand(cols >= 0, "invalid column count {cols}");
    demand(rows >= 0, "invalid row count {rows}");
    demand(maxval > 0, "invalid {maxval}");
    demand(NN_max > 0, "invalid {NN_max}");

    uint16_color_table_t *cht = talloc(1, uint16_color_table_t);
    cht->buckets = talloc(HASH_SIZE, uint16_color_table_node_t*);
    cht->maxval = maxval;
    cht->NN = 0;

    /* Go through the {samples} arrays, storing the colors in the table. */
    for (int32_t row = 0; row < rows; ++row)
      { uint16_t *sP = samples[row];
        for (int32_t col = 0; col < cols; ++col)
          { /* Grab a pixel {p}, expanding/truncating to 3 samples: */ 
            uint16_t *clr = sP;
            sP += chns;
            bool_t add = cht->NN < NN_max;
            uint16_color_table_node_t *nd = uint16_color_table_lookup(cht, chns, clr, add);
            if (nd == NULL)
              { /* Not found, and could not add because too many colors: */
                assert(cht->NN >= NN_max);
                uint16_color_table_free(cht); 
                return NULL; 
              }
          }
      }
    assert(cht->NN <= NN_max);
    return cht;
  }

uint16_color_table_node_t** uint16_color_table_nodes(uint16_color_table_t *cht)
  { uint32_t NN = cht->NN;
    uint16_color_table_node_t **chv = talloc(NN, uint16_color_table_node_t *);
    uint32_t k = 0;
    for (uint32_t i = 0; i < HASH_SIZE; ++i)
      { uint16_color_table_node_t *nd = cht->buckets[i];
        while (nd != NULL)
          { /* Add the new entry. */
            assert(k < NN);
            chv[k] = nd; k++;
            nd = nd->next;
          }
      }
    assert(k == NN);
    return chv;
  }

uint16_t* uint16_color_table_samples_from_nodes(uint32_t chns, uint32_t NN, uint16_color_table_node_t** nds)
  { uint32_t NS = chns*NN;
    uint16_t *smp = talloc(NS, uint16_t);
    uint16_t *p = smp;
    for (int32_t i = 0; i < NN; i++)
      { uint16_color_table_node_t *nd = nds[i];
        for (int32_t c = 0; c < chns; c++)
          { (*p) = nd->clr[c]; p++; }
      }
    assert(p - smp == NS);
    return smp;
  }

void uint16_color_table_free(uint16_color_table_t *cht)
  { uint32_t NF = 0; /* Total nodes reclaimed. */
    for (int32_t i = 0; i < HASH_SIZE; ++i)
      { uint16_color_table_node_t *nd = cht->buckets[i];
        while (nd != NULL)
          { uint16_color_table_node_t *ndnext = nd->next;
            free(nd); NF++;
            nd = ndnext;
          }
      }
    assert(NF == cht->NN);
    free(cht->buckets);
    free(cht);
  }
