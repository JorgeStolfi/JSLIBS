/* Huffman tree and Huffman decoding for the NASA/JPL/PDS IMQ image format. */
/* Last edited on 2023-02-07 22:07:12 by stolfi */

#ifndef imq_huff_H
#define imq_huff_H

/* 
  The original images produced by the Voyager and Viking missions are
  available in a special IMQ file format defined by NASA/JPL/PDS. In
  these files, each scanline is compressed by using Hufffman encoding
  for the differences between successive pixels. This module provides
  tools for decompressing IMQ images.

  The Huffman encoding is sensitive to the tie-breaking criterion used
  when sorting nodes with the same weight, and when choosing which child
  of an internal node will be the left one. The procedures
  {*imq_huff_build_tree} and {imq_huff_decode} below are designed
  to provide the same bit decoding and scanline integration as the
  {huff_tree} and {DECOMP.C} routines used by NASA/JPL/PDS, e.g. in Mike
  Martin {vdcomp.c} program (1989).
  
  This module considers only the image pixels. It does not provide tools 
  for parsing the IMQ image headers and the VAX-style variable-length records
  that comprise the file.  See the separate program {imqtopgm.c}
  for that. */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <codetree.h>
#include <huff_tree.h>

codetree_t *imq_huff_build_tree(huff_tree_freq_t freq[]);
  /* Builds the Huffman tree {T} for a set of pixel differences
    based on the frequency counts {freq[0..510]}. Returns a pointer 
    to the root node. 
    
    The image pixels are integers in {0.255}, so the difference between
    two pixels {d[i]} is in {-255..+255}. The array {freq} must have 
    exactly 511 elements, and stores the histogram of pixel differences
    observed in the sample, shifted to {0..510}.
    
    More precisely, the value of {freq[ih]} is assumed to be the number
    of times that {pix[i]-pix[i+1]} (in that order) was {ih - 255}. That
    is, {freq[0]} is the number of times a pixel with value 0 is
    followed by one with value 255; {freq[510]} is the number of times
    the reverse happened; and {pix[255]} is the number of times two
    equal pixels occurred next to each other.
    
    This tree can be used with {imq_huff_decode_line} below to
    reconstruct a sequence of pixels given the first pixel and the
    compressed differences between successive pixels.
    
    The numbers stored in the {p->value} fields of the tree leaves are
    actually the complements of the shifted differences above, namely
    {510-ih} intead of {ih}. This changed was needed in order to get the
    proper tie-breaking in the Huffman algortithm, and is taken into
    account by the other procedures in this module. */
  
codetree_bit_count_t imq_huff_decode
  ( uint32_t line_num,
    codetree_byte_count_t nb, 
    byte_t buf[], 
    codetree_t *tree, 
    codetree_byte_count_t np, 
    byte_t pix[]
  );
  /* Decodes the bit sequence stored in the byte array {buf[0..nb-1]}
    into a sequence of {np} 8-bit pixel values {pix[0..np-1]}, 
    using the differences encoding defined by the given {tree}.
    
    The first pixel {pix[0]} will be {buf[0]}. The bits of
    {buf[1..nb-1]} are then parsed as Huffman-encoded pixel differences,
    and each difference is used to compute each successive pixel, until
    {np} pixels are obtained in all.
    
    The procedure fails if the bytes {buf[1..nb-1]} are all parsed
    before completing {np} pixels.
    
    The line number {line_num} (starting from 1, by NASA/JPL convention)
    is used only for diagnostics.
    
    The bits in the IMQ Huffman coding mean '0' go to the right child,
    '1' to the left child. */

void imq_huff_print_codes(FILE *stderr, codetree_t *tree);
  /* Prints the Huffman codes for each pixel-pixel difference
    implied by the given {tree}.  */

#endif
