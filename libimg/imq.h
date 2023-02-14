/* Tools for reading NASA/JPL PDS format images. */
/* Last edited on 2023-02-07 22:01:45 by stolfi */

#ifndef imq_H
#define imq_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

typedef uint8_t byte_t;

void imq_read_header
  ( FILE *rd, 
    char *format_P, 
    uint32_t *record_bytes_P, 
    uint32_t *label_checksum_P,
    uint32_t *NX_P, 
    uint32_t *NY_P, 
    bool_t verbose
  );
  /* Reads and parses the PDS file header (aka "labels") from the input file {rd}.
    
    Extracts the following information from the headers:
    
       the number of bytes per record {record_bytes}
       
       the declared checksum {label_checksum} 
       
       the number of samples per line {NX}
       
       the number of lines {NY}
       
    Returns those fields in {*label_checksum_P}, {*record_bytes_P},
    {*NX_P}, {*NY_P}.
    
    Also sets the character {*format_P} to 'Y' if it is a Voyager image,
    'K' if it is a Viking image.  They have somewhat different formats
    for the image and encoding histograms, among other things. */

uint32_t imq_read_var_length_record(FILE *rd, uint32_t nb, byte_t ibuf[]);
  /* Reads one record from the from input file {rd}.
    
    Assumes that each record begins with a two-byte unsigned integer {length}
    in small-endian format, then {length} bytes, then possibly one more filler byte
    to make the total bytes even.
    
    Stores the {length} bytes read into {ibuf[0..length-1]}, and returns
    the value of {length}. Assumes that {ibuf} has at least {nb}
    characters, and fails if {length > nb}. The extra filler byte, if
    present, is read but not stored.
    
    If end-of-file occurs when reading the first byte of the length
    field, the procedure returns {length = 0} and does not change
    {ibuf}. Otherwise it fails with an error message. */
  
#endif
