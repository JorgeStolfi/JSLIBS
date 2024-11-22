/* Binary tree based encoding and decoding. */
/* Last edited on 2024-11-20 01:49:12 by stolfi */

#ifndef codetree_H
#define codetree_H

#include <stdint.h>
#include <stdio.h>

#include <codetree_limits.h>

/* CODE TREES

  A /code tree/ is a data structure that can be used to parse and decode
  a stream of bits into a sequence of non-negative integer values from
  some set {V}, with possibly different number of bits used for
  different values. See {codetree_encode} and {codetree_decode}
  below. */

typedef uint32_t codetree_data_value_t;
  /* Type of an element from the domain set {V}.
    Actually in the range {0..codetree_data_MAX_VALUE}. */

typedef int32_t codetree_node_value_t;
  /* Type of the {value} field of a code tree node, which may be a
  {codetree_data_value_t} or a negative value indicating an internal
  node. Actually in
  {codetree_node_MIN_VALUE..codetree_node_MAX_VALUE} */

typedef struct codetree_node_t
  { codetree_node_value_t value; 
    /* These fields may not exist: */
    struct codetree_node_t *child[2];
  } codetree_node_t;
  /* A {codetree_node_t} {p} is the root of a tree or subtree
    that represents a finite set {p.V} of non-negative integers.
    
    If the field {p.value} of a node {p} is non-negative, the node is a leaf, 
    and the set {p.V} consists of {p.value} only.
    
    If {p.value} is negative, then {p} is an internal node and has two
    non-null children {p.child[0..1]}. The set {p.V} is the union of
    {p.child[0].V} and {p.child[1].V}, with should be disjoint.
    
    The {child} fields of a leaf are irrelevant and should not
    be accessed since they may not have been allocated in memory. */

typedef codetree_node_t codetree_t;
  /* As usual, a pointer to a code tree is a pointer to its root node. */

typedef uint32_t codetree_node_count_t;
  /* A type large enough to store the number of nodes in a code tree. */

codetree_node_t *codetree_new_internal
  ( codetree_node_count_t seq,
    codetree_node_t *child0,
    codetree_node_t *child1
  );
  /* Allocates a new non-leaf node, with the given children. Returns its address. 
    
    The given {child0} and/or {child1} may be {NULL}. Note however that 
    those links must be non-null in a finished valid tree.
    
    The parameter {seq} is the order of creation of the node, starting
    from 0. The {value} field will be set to {codetree_node_MIN_VALUE + seq}. It
    should be at most {codetree_MAX_INTERNALS-1}. */

codetree_node_t *codetree_new_leaf(codetree_data_value_t val);
  /* Allocates a new leaf node with {val} as the {value} field. Returns
    its address. Fails if {val} is {NO_VALUE}. The {child} fields of the
    node should not be accessed as they may not exist. */
  
codetree_node_count_t codetree_free(codetree_t *tree);
  /* Reclaims the storage (with {free}) of all nodes of the tree
    with root node {*tree}. Returns the count of nodes reclaimed.
    If {tree} is {NULL}, does nothing and returns 0.*/
 
codetree_node_count_t codetree_num_leaves(codetree_t *tree);
  /* Returns the number of leaf nodes in the subtree with 
    root node {tree}. */

/* ENCONDING AND DECODING */

typedef uint64_t codetree_sample_count_t;
  /* A type hopefully large enough to store the number of samples
    to be encoded or decoded. */

typedef uint64_t codetree_bit_count_t;
  /* A type hopefully large enough to store the number of bits used 
    to encode any sample sequence. */
  
typedef uint64_t codetree_byte_count_t;
  /* A type hopefully large enough to store the number of bits in any encoded
    value sequence. */
  
typedef uint8_t byte_t;

codetree_bit_count_t codetree_decode
  ( codetree_byte_count_t nb, 
    byte_t buf[], 
    codetree_t *tree, 
    codetree_data_value_t maxval, 
    codetree_sample_count_t ns, 
    codetree_data_value_t smp[]
  );
  /* Decodes the bit sequence stored in the byte array {buf[0..nb-1]}
    into a sequence of {ns} values {smp[0..ns-1]}, using the 
    encoding described by the given {tree}.
    
    DECODING PROCESS
    
    The bytes of {buf} are examined in sequence starting with {buf[0]}. The
    bits of each each byte are examined in sequence from higher (bit 7,
    numeric value {2^7 = 128}) to lower (bit 0, value {2^0 = 1}). If
    only {k} of the 8 bits of byte {buf[nu-1]} are used, they are the
    highest-order ones.
    
    The bit string is parsed into codes for the values in the set
    {tree.V}, possibly with different lengths. The bit code representing
    a value {val} from {tree.V} is the description of the path from the
    root node of {tree} to the leaf {q} with {q.value = val}, where each
    0 means "go to {child[0]}" and each 1 means "go to {child[1]}".
    
    Whenever the procedure reaches a leaf node {q}, it appends the
    corresponding sample {q.value} to {smp}, and parsing starts again
    from the root node of the {tree}; until all the {ns} requested samples have been
    decoded.
    
    The procedure returns the number {nu} of BITS of {buf} that were actually
    scanned.  The number of BYTES used will then be {(nu+7)/8}.
    
    SPECIAL CASES
    
    If {ns} is zero, the procedure will ignore the other arguments
    and return zero.
    
    If {ns} is positive, the set {V} must have at least one value, that
    is, {tree} must not be {NULL}. 
    
    If {tree} is a leaf node, then {V} has a single value {val =
    p->value}, that is encoded with zero bits. In this case the
    procedure will ignore {buf}, will fill {smp[0..ns-1]} with {val},
    and will return zero.
    
    The procedure fails if all the bytes {buf[0.nb-1]} are scanned
    before extracting {ns} complete value codes; or if any of the values
    stored in the tree exceeds {maxval} (which itself must be in
    {0..codetree_data_MAX_VALUE}), or if the tree is malformed in some way. */

typedef uint32_t codetree_delta_t;
  /* A type large enought to contain any element in the {delta}
    encoding table. */

codetree_bit_count_t codetree_encode
  ( codetree_sample_count_t ns,
    codetree_data_value_t smp[], 
    codetree_data_value_t maxval, 
    codetree_node_count_t nd,
    codetree_delta_t delta[], 
    codetree_byte_count_t nb, 
    byte_t buf[]
  );
  /* Encodes the sequence of values {smp[0..ns-1]} as a bit sequence and
    stores it into the byte array {buf[0..nb-1]}, in a way compatible
    with {codetree_decode}. Returns the number {nu} of bits actually used.
    
    ENCODING PROCESS
    
    The bytes of {buf} are filled in sequence starting with {buf[0]}. The
    bits of each each byte are defined in sequence from higher (bit 7,
    numeric value {2^7 = 128}) to lower (bit 0, value {2^0 = 1}). If
    only {k} of the 8 bits of byte {buf[nu-1]} are set, its lowest
    {8-k} bits are set to zero.
    
    The encoding is described by the table {delta[0..nd-1]}. It size
    {nd} must be at least {nd = maxval+|V|}, where {V} is the set of
    allowed sample values (which must be a subset of {0..maxval}). Each
    element of {delta[0..np-1]} must be a number in {0..2*nd-1}.
    
    To encode a sample value {val} from {smp[0..ns-1]}, the procedure
    begins with the index {iv = val} and, while {delta[iv]} is not
    zero, repeately executes {r = delta[iv]%2, iv = iv +
    parent[iv]/2}, until {delta[iv]} is zero. The bits {r} obtained in
    this way, IN REVERSE ORDER, are then appended o the bit stream.
    
    If {delta[val]} is zero for some {val} in {0..maxval}, then {val} is
    not in the representable set {V}. 
     
    The procedure returns the number {nu} of BITS of {buf} that were
    actually used. The number of BYTES used will then be {(nu+7)/8}.
    
    SPECIAL CASES
    
    The procedure fails if the bytes {buf[0..nb-1]} are not enough to
    encode all the values {smp[0..ne-1]}, or any of those values exceeds
    {maxval} (which in turn must not exceed {codetree_data_MAX_VALUE}), or
    any value {smp[k]} is not representable (that is, if delta[smp[k]]
    is zero).
    
    If {ns} is zero, the procedure ignores the other arguments
    and returns zero.
    
    If {ns} is positive, the set {V} must have at least one valid 
    value; that is, at least one of {delta[0..maxval]} must be 
    nonzero. The procedure then will check whether all samples
    {smp[0..ns-1]} are in {V}.
    
    If {V} has only one valid value {val}, its code has zero bits.
    Then the procedure will not store any bits in {buf} and 
    return zero. */

codetree_node_count_t codetree_get_encoding_table
  ( codetree_t *tree, 
    codetree_data_value_t maxval, 
    codetree_node_count_t nd,
    codetree_delta_t delta[]
  );
  /* Stores in {delta[0..nd-1]} the encoding table corresponding to the
    given decoding {tree}.  Returns the number of entries
    of {delta} actually used, which may be less than {nd}. 
    
    The actual table size will be {maxval+|V|} where {|V|} is the 
    number of representable values (the number of leaf nodes in the
    tree).  Thus the table will fit in {delta[0..nd-1]} if
    {nd = 2*(maxval+1)}.
    
    Fails if any leaf node has a {value} field greater than {maxval},
    of if there are two leaf nodes with the same value,
    of if the table needs more than {nd} elements, or the
    tree is malformed. */ 

codetree_t *codetree_get_decoding_tree
  ( codetree_node_count_t nd,
    codetree_delta_t delta[], 
    codetree_data_value_t maxval 
  );
  /* Builds the decoding tree that corresponds to the encoding 
    table {delta[0..nd-1]}, and returns its root.
    
    Assumes that the set {V} of allowable values consists of every
    integer {val} in {0..maxval} such that {delta[val]} is not zero. 
    The actual table size shoudl be {md = maxval + |V|}.
    The procedure faisl if {nd} is less than this value. 
    Otherwise it examines only entries {delta[0..md-1]}. */ 
 
/* VALIDATION */

void codetree_check_tree(codetree_t *tree, codetree_data_value_t maxval);
  /* Verifies the consistency of the given decoding {tree}, which should
    have leaf values in {0..maxval} only. */

void codetree_check_table(codetree_node_count_t nd, codetree_delta_t delta[], codetree_data_value_t maxval);
  /* Verifies the consistency of the encoding table {delta[0..nd-1]},
    which should allow encoding of some values in {0..maxval} only. */

void codetree_check_iso(codetree_t *p, codetree_t *q);
  /* Checks if the subtrees with roots {p} and {q} are 
    isomorphic, meaning that they assign the same codes to the same
    set of values.  The procedure fails with error if they are not. */

codetree_node_count_t codetree_check_codes(codetree_t *tree, codetree_data_value_t maxval, char *code[]);
  /* Every element of {code[0..maxval]} must be either {NULL} or a
    string with characters '0' or '1'. Checks whether {code[val]} is not
    {NULL} if and only if {val} is a member of the set {tree.V}, and
    then checks whether that string is the bit code implied by {tree}
    for the value {val}. If it succeeds, returns the number of
    non-{NULL} codes, that is, the number {|tree.V|} of leaf nodes in
    {tree}. Otherwise it fails with an error message. */

codetree_node_count_t codetree_print_codes(FILE *wr, codetree_t *tree);
  /* Prints to {wr} the set of values {V} stored in the leaves of 
    the given {tree}, and the binary sequences that describe
    the path from {tree}'s root to each leaf.  Returns the number {|V|}
    of valid values. */

void codetree_print_bits(FILE *wr, codetree_byte_count_t nb, byte_t buf[], char *sep);
  /* Prints to {wr} the bits of the bytes {buf[0..nb-1]}, 
    from bit 7 (value 128) to bit 0 (value 1), with the string 
    {sep} between consecutive bytes. */

#endif
