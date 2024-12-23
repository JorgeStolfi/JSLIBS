/* See jsaudio_dvf.h */
/* Last edited on 2024-12-21 03:38:48 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <jsaudio_io.h>
#include <jsaudio_dvf.h>

/* INTERNAL PROTOTYPES */

jsaudio_dvf_time_t jsaudio_read_dvf_time(FILE *rd);
  /* Reads a date and time (7 bytes) from file {rd}. */

/* IMPLEMENTATIONS */

jsaudio_dvf_file_header_t jsaudio_dvf_read_file_header(FILE *rd)
  {
    jsaudio_dvf_file_header_t h;
    
    /* Bytes 0..7: */
    char magic[8];                        /* Magic ID of DVF file. */ 
    jsaudio_read_chars(rd, 8, magic);         /* Magic ID of DVF file. */
    assert(strncmp(magic, "MS_VOICE", 8) == 0); 
    
    /* Bytes 8..11: */
    uint32_t unk_008_4 = jsaudio_read_uint32_be(rd); 
    assert(unk_008_4 == 80); 
    
    /* Byte 12 */
    uint8_t unk_012_1 = jsaudio_read_uint8(rd);
    assert(unk_012_1 == 1);
    
    /* Byte 13 */
    uint8_t unk_013_1 = jsaudio_read_uint8(rd);
    assert(unk_013_1 == 2);
    
    /* Bytes 14..15 */
    uint16_t unk_014_2 = jsaudio_read_uint16_be(rd);
    assert(unk_014_2 == 0);
    
    /* Bytes 16..31 */
    char corp1[16];
    jsaudio_read_chars(rd, 16, corp1);
    assert(strncmp(corp1, "SONY CORPORATION", 16) == 0);
    
    /* Byte 32 */
    uint8_t unk_032_1 = jsaudio_read_uint8(rd);
    assert(unk_032_1 == 80);
    
    /* Bytes 33..44 */
    char dll[12];
    jsaudio_read_chars(rd, 12, dll);
    assert(strncmp(dll,"Voice.dll\0\0\0",12) == 0);
    
    /* Bytes 45..48 */
    uint32_t unk_045_4 = jsaudio_read_uint32_be(rd);
    assert(unk_045_4 == 1);
    
    /* Byte 49 */
    uint8_t unk_049_1 = jsaudio_read_uint8(rd);
    assert(unk_049_1 == 0);
    
    /* Byte 50 */
    uint8_t unk_050_1 = jsaudio_read_uint8(rd);
    assert(unk_050_1 == 0);
    
    /* Byte 51 */
    uint8_t unk_051_1 = jsaudio_read_uint8(rd);
    assert(unk_051_1 == 0);
    
    /* Bytes 52..58 */
    h.time_1 = jsaudio_read_dvf_time(rd);

    /* Byte 59 */
    h.unk_059_1 = jsaudio_read_uint8(rd);
    assert(h.unk_059_1 < 4);
    
    /* Byte 60 */
    uint8_t unk_060_1 = jsaudio_read_uint8(rd);
    assert(unk_060_1 == 0);

    /* Byte 61 */
    uint8_t unk_061_1 = jsaudio_read_uint8(rd);
    assert(unk_061_1 == 53);
    
    /* Byte 62 */
    uint8_t unk_062_1 = jsaudio_read_uint8(rd);
    assert(unk_062_1 == 0);
    
    /* Byte 63 */
    uint8_t unk_063_1 = jsaudio_read_uint8(rd);
    assert(unk_063_1 == 1);
    
    /* Byte 64 */
    uint8_t unk_064_1 = jsaudio_read_uint8(rd);
    assert(unk_064_1 == 0);
    
    /* Byte 65 */
    uint8_t unk_065_1 = jsaudio_read_uint8(rd);
    assert(unk_065_1 == 0);
    
    /* Byte 66 */
    uint8_t unk_066_1 = jsaudio_read_uint8(rd);
    assert(unk_066_1 == 28);
    
    /* Byte 67 */
    uint8_t unk_067_1 = jsaudio_read_uint8(rd);
    assert(unk_067_1 == 32);
    
    /* Byte 68 */
    uint8_t unk_068_1 = jsaudio_read_uint8(rd);
    assert(unk_068_1 == 0);
    
    /* Byte 69 */
    uint8_t unk_069_1 = jsaudio_read_uint8(rd);
    assert(unk_069_1 == 0);
    
    /* Byte 70 */
    uint8_t unk_070_1 = jsaudio_read_uint8(rd);
    assert(unk_070_1 == 3);
    
    /* Byte 71 */
    uint8_t unk_071_1 = jsaudio_read_uint8(rd);
    assert(unk_071_1 == 132);
    
    /* Byte 72 */
    uint8_t unk_072_1 = jsaudio_read_uint8(rd);
    assert(unk_072_1 == 4);
    
    /* Byte 73 */
    uint8_t unk_073_1 = jsaudio_read_uint8(rd);
    assert(unk_073_1 == 0);
    
    /* Byte 74 */
    uint8_t unk_074_1 = jsaudio_read_uint8(rd);
    assert(unk_074_1 == 0);
    
    /* Byte 75 */
    uint8_t unk_075_1 = jsaudio_read_uint8(rd);
    assert(unk_075_1 == 1);
    
    /* Byte 76..79 */
    uint32_t unk_076_4 = jsaudio_read_uint32_be(rd);
    assert(unk_076_4 == 0);
    
    /* Byte 80 */
    uint8_t unk_080_1 = jsaudio_read_uint8(rd);
    assert(unk_080_1 == 1);
    
    /* Bytes 81..88 */
    uint64_t unk_081_8 = jsaudio_read_uint64_be(rd);
    assert(unk_081_8 == 20482);
    
    /* Bytes 89..96 */
    uint64_t unk_089_8 = jsaudio_read_uint64_be(rd);
    assert(unk_089_8 == 24581);
    
    /* Bytes 97..104 */
    uint64_t unk_097_8 = jsaudio_read_uint64_be(rd);
    assert(unk_097_8 == 65542);
    
    /* Bytes 105..112 */
    uint64_t unk_105_8 = jsaudio_read_uint64_be(rd);
    assert(unk_105_8 == 4103);
    
    /* Bytes 113..120 */
    uint64_t unk_113_8 = jsaudio_read_uint64_be(rd);
    assert(unk_113_8 == 8201);
    
    /* Bytes 121..128 */
    uint64_t unk_121_8 = jsaudio_read_uint64_be(rd);
    assert(unk_121_8 == 4106);
    
    /* Bytes 129..136 */
    uint64_t unk_129_8 = jsaudio_read_uint64_be(rd);
    assert(unk_129_8 == 4107);
    
    /* Bytes 137..144 */
    uint64_t unk_137_8 = jsaudio_read_uint64_be(rd);
    assert(unk_137_8 == 4109);
    
    /* Bytes 145..152 */
    uint64_t unk_145_8 = jsaudio_read_uint64_be(rd);
    assert(unk_145_8 == 4110);
    
    /* Bytes 153..160 */
    uint64_t unk_153_8 = jsaudio_read_uint64_be(rd);
    assert(unk_153_8 == 262147);
    
    /* Bytes 161..168 */
    uint64_t unk_161_8 = jsaudio_read_uint64_be(rd);
    assert(unk_161_8 == 122884);

    /* Byte 169 */
    uint8_t unk_169_1 = jsaudio_read_uint8(rd);
    assert(unk_169_1 == 0);
    
    /* Byte 170 */
    h.unk_170_1 = jsaudio_read_uint8(rd);
    assert(h.unk_170_1 < 4);
    
    /* Byte 171 */
    h.unk_171_1 = jsaudio_read_uint8(rd);
    /* Arbitrary value? */
    
    /* Byte 172 */
    uint8_t unk_172_1 = jsaudio_read_uint8(rd);
    assert(unk_172_1 == 0);
    
    /* Bytes 173..174 */
    h.unk_173_2 = jsaudio_read_uint16_be(rd);
    /* Arbitrary value? */
    
    /* Bytes 175..176 */
    uint16_t unk_175_2 = jsaudio_read_uint16_be(rd);
    assert(unk_175_2 == 144);
    
    /* Bytes 177..431 */
    for (int32_t i = 177; i <= 431; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Byte 432 */
    uint8_t unk_432_1 = jsaudio_read_uint8(rd);
    assert(unk_432_1 == 1);
    
    /* Bytes 433 */
    uint8_t unk_433_1 = jsaudio_read_uint8(rd);
    assert(unk_433_1 == 0);
    
    /* Bytes 434..447 */
    char corp2[14];
    jsaudio_read_chars(rd, 14, corp2);
    assert(strncmp(corp2, "SONY CORP.    ", 14) == 0);
    
    /* Byte 448 */
    uint8_t unk_448_1 = jsaudio_read_uint8(rd);
    assert(unk_448_1 == 144);
    
    /* Bytes 449..479 */
    for (int32_t i = 449; i <= 479; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Bytes 480..483 */
    h.length_su = jsaudio_read_uint32_be(rd);
    assert(h.length_su % 216 == 0);
    
    /* Bytes 484..495 */
    for (int32_t i = 484; i <= 495; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Bytes 496..499 */
    uint32_t unk_496_4 = jsaudio_read_uint32_be(rd);
    assert(unk_496_4 == 255);
    
    /* Bytes 500..544 */
    for (int32_t i = 500; i <= 544; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Byte 545 */
    h.unk_545_1 = jsaudio_read_uint8(rd);
    /* Arbitrary value? */
    
    /* Byte 546 */
    uint8_t unk_546_1 = jsaudio_read_uint8(rd);
    assert(unk_546_1 == 0);
    
    /* Byte 547 */
    uint8_t unk_547_1 = jsaudio_read_uint8(rd);
    assert(unk_547_1 == 255);
    
    /* Bytes 548..554 */
    h.time_2 = jsaudio_read_dvf_time(rd);

    /* Bytes 555 */
    h.unk_555_1 = jsaudio_read_uint8(rd);
    assert(h.unk_555_1 < 4);
    
    /* Bytes 556..565 */
    jsaudio_read_chars(rd, 24, h.model);
    assert(strncmp(h.model, "SONY ICD-P\0", 11) == 0);
    
    /* Bytes 566..579 */
    for (int32_t i = 566; i <= 579; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Byte 580 */
    uint8_t unk_580_1 = jsaudio_read_uint8(rd);
    assert(unk_580_1 == 34);
    
    /* Byte 581 */
    uint8_t unk_581_1 = jsaudio_read_uint8(rd);
    assert(unk_581_1 == 5);
    
    /* Byte 582 */
    uint8_t unk_582_1 = jsaudio_read_uint8(rd);
    assert(unk_582_1 == 7);
    
    /* Byte 583 */
    uint8_t unk_583_1 = jsaudio_read_uint8(rd);
    assert(unk_583_1 == 27);
    
    /* Byte 584 */
    uint8_t unk_584_1 = jsaudio_read_uint8(rd);
    assert(unk_584_1 == 204);
    
    /* Byte 585 */
    uint8_t unk_585_1 = jsaudio_read_uint8(rd);
    assert(unk_585_1 == 166);
    
    /* Byte 586 */
    uint8_t unk_586_1 = jsaudio_read_uint8(rd);
    assert(unk_586_1 == 172);
    
    /* Byte 587 */
    uint8_t unk_587_1 = jsaudio_read_uint8(rd);
    assert(unk_587_1 == 197);
    
    /* Bytes 588..595 */
    for (int32_t i = 588; i <= 595; i++) { int32_t chr = fgetc(rd); assert(chr == 0); }
    
    /* Bytes 596..599 */
    for (int32_t i = 596; i <= 599; i++) { int32_t chr = fgetc(rd); assert(chr == 255); }
    
    /* Consume stream until first data byte: */
    /* for (int32_t i = 600; i < file_header_size; i++) { int32_t chr = fgetc(rd); assert(chr != EOF); } */
      
    return h;
  }

jsaudio_dvf_time_t jsaudio_read_dvf_time(FILE *rd)
  {
    jsaudio_dvf_time_t t;
    
    t.year = jsaudio_read_uint16_be(rd);
    assert(t.year >= 2000);
    assert(t.year <= 2020);
    
    t.month = jsaudio_read_uint8(rd); 
    assert(t.month >= 1);
    assert(t.month <= 12);
    
    t.day = jsaudio_read_uint8(rd); 
    assert(t.day >= 1);
    assert(t.day <= 31);
    
    t.hour = jsaudio_read_uint8(rd); 
    assert(t.hour <= 23);
    
    t.minute = jsaudio_read_uint8(rd); 
    assert(t.minute <= 59);
    
    t.second = jsaudio_read_uint8(rd); 
    assert(t.second <= 59);
    
    return t;
  }
