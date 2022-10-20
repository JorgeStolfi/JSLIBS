/* Last edited on 2022-10-20 06:03:50 by stolfi */

/*
Original code {tinyxml2.hpp} by Lee Thomason (www.grinninglizard.com)
Modified by J. Stolfi, UNICAMP 2015-09-26

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any
damages arising from the use of this software.

Permission is granted to anyone to use this software for any
purpose, including commercial applications, and to alter it and
redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.


2. Altered source versions must be plainly marked as such, and
must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.
*/

#ifndef salamic_xml_H
#define salamic_xml_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <bool.h>

int32_t salamic_xml_snprintf( char* buffer, size_t size, const char* format, ... );

// namespace salamic_xml
// {
// class XMLDocument;
// class XMLElement;
// class XMLAttribute;
// class XMLComment;
// class XMLText;
// class XMLDeclaration;
// class XMLUnknown;
// class XMLPrinter;

typedef enum
  {
    salamic_xml_NEEDS_ENTITY_PROCESSING         = 0x01,
    salamic_xml_NEEDS_NEWLINE_NORMALIZATION     = 0x02,
    salamic_xml_COLLAPSE_WHITESPACE             = 0x04,

    salamic_xml_NEEDS_FLUSH                     = 0x100,
    salamic_xml_NEEDS_DELETE                    = 0x200,

    salamic_xml_TEXT_ELEMENT                    = salamic_xml_NEEDS_ENTITY_PROCESSING | salamic_xml_NEEDS_NEWLINE_NORMALIZATION,
    salamic_xml_TEXT_ELEMENT_LEAVE_ENTITIES     = salamic_xml_NEEDS_NEWLINE_NORMALIZATION,
    salamic_xml_ATTRIBUTE_NAME                  = 0,
    salamic_xml_ATTRIBUTE_VALUE                 = salamic_xml_NEEDS_ENTITY_PROCESSING | salamic_xml_NEEDS_NEWLINE_NORMALIZATION,
    salamic_xml_ATTRIBUTE_VALUE_LEAVE_ENTITIES  = salamic_xml_NEEDS_NEWLINE_NORMALIZATION,
    salamic_xml_COMMENT                         = salamic_xml_NEEDS_NEWLINE_NORMALIZATION
  } salamic_xml_str_type_t;



typedef struct salamic_xml_str_pair_t
  {
    char *start;
    char *end;
    int32_t flags;
  } salamic_xml_str_pair_t;
  /*
    A strig wrapper. Normally stores the start and end
    pointers into the XML file itself, and will apply normalization
    and entity translation if actually read. Can also store (and memory
    manage) a traditional {char[]}. */

void salamic_xml_str_pair_set(salamic_xml_str_pair_t *str, char* start, char* end, int32_t flags);
// StrPair.Set

char *salamic_xml_str_pair_get_chars (salamic_xml_str_pair_t *str);
// StrPair.GetChars

bool_t salamic_xml_str_pair_is_empty(salamic_xml_str_pair_t *str);
// StrPair.Empty

void salamic_xml_str_pair_from_string(salamic_xml_str_pair_t *str, const char* ch, int32_t flags);
// StrPair.SetChars

char* salamic_xml_str_pair_parse_text(salamic_xml_str_pair_t *str, char* in, const char* endTag, int32_t flags);
// StrPair.ParseText

char* salamic_xml_str_pair_parse_name(salamic_xml_str_pair_t *str, char* in);
// StrPair.ParseName

void salamic_xml_str_pair_reset(salamic_xml_str_pair_t *str);
// StrPair.Reset

void salamic_xml_str_pair_collapse_white_space(salamic_xml_str_pair_t *str);
// StrPair.CollapseWhitespace

#endif
