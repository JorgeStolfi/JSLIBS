/* {jspng.h} - definitions and tools for PNG image files.  */
/* Last edited on 2024-12-05 10:30:51 by stolfi */

#ifndef jspng_H
#define jspng_H

#include <png.h>

#define jspng_MAX_CHANS 4
  /* Max channels expected in a PNG image file. */

#define jspng_MAGIC_BYTES 8
  /* Number of bytes in the PNG image file signature. */

#define jspng_MAX_SIZE (2147483647u)
  /* Max width and height expected in a PNG image file. */

void jspng_dump_info(FILE *wr, const char *func, char *label, png_structp pr, png_infop pi);
  /* Writes to {wr} the main attributes of an image as stored in the
    records {pr} and {pi}, in a human-readable format. The strings
    {func} and {label} are printed in the header of the printout. */

#endif
