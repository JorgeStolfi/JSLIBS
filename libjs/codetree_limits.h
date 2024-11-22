/* Limits for {codetree.h} parameters. */
/* Last edited on 2024-11-20 01:56:23 by stolfi */

#ifndef codetree_limits_H
#define codetree_limits_H

#include <stdint.h>
#include <stdio.h>

/* 
  These limits are somewhat arbitrary but are expected to avoid overflow
  bugs in parameter passing and the implementation.
  They should not be increased without a careful revision of the code. */

#define codetree_data_MAX_VALUE  (((int32_t)1) << 29)
  /* The set {V} of values that may be encoded or decoded using a code
    tree must be a subset of {0..codetree_data_MAX_VALUE}. The value is quite
    a bit less than {INT32_MAX} to help avoid overflow bugs. */

#define codetree_MAX_LEAVES ((uint64_t)(codetree_data_MAX_VALUE + 1))
  /* The maximum number of leaves in any tree. */

#define codetree_MAX_INTERNALS ((uint64_t)(codetree_MAX_LEAVES - 1))
  /* The maximum number of internal nodes in any tree. */

#define codetree_node_MAX_VALUE  (codetree_data_MAX_VALUE)
#define codetree_node_MIN_VALUE  (-((int32_t)codetree_MAX_INTERNALS - 1))
  /* The {value} field of a {codetree_node_t} can be in
    {codetree_node_MIN_VALUE..codetree_NODE_MAX_VALUE}. 
    negative values denote internal nodes of the tree. */

#define codetree_MAX_NODES ((uint64_t)(codetree_MAX_LEAVES + codetree_MAX_INTERNALS))
  /* The maximum number of internal nodes in any tree. */
  
#define codetree_MAX_SAMPLES ((uint64_t)((1LU << 51) - 1))
  /* The maximum number of samples that can be decoded in one go.
    Somewhat arbitrary. */
  
#define codetree_MAX_BITS ((uint64_t)((1LU << 54) - 1))
  /* Hopefully enough for most uses. */
  
#define codetree_MAX_BYTES ((uint64_t)((1LU << 51) - 1))
  /* Hopefully enough for most uses. */
    
#define codetree_MAX_DELTA ((uint64_t)(2*codetree_MAX_NODES - 1))
  /* Max value stored in a {delta} table. */

#endif
