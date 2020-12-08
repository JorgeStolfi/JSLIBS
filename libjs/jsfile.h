#ifndef jsfile_H
#define jsfile_H

/* File open with auto bomb-out; read line as new string. */
/* Last edited on 2016-04-01 00:45:56 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

FILE *open_read(const char *name, bool_t verbose);
FILE *open_write(const char *name, bool_t verbose);
  /* Opens file for input or output, respectively. If {name} is "-",
    returns {stdin} or {stdout}, respectively. If {verbose} is TRUE,
    prints a message to {stderr}. Bombs out in case of failure
    (including {name} is NULL or an empty string). */

FILE *open_read_tag_ext(char *name, char *tag, char *ext, bool_t verbose);
FILE *open_write_tag_ext(char *name, char *tag, char *ext, bool_t verbose);
  /* These procedure open a disk file called "{name}{tag}{ext}",
    for reading or writing, respectively.  If that file name is "-",
    returns {stdin} or {stdout}, respectively.  A message is printed
    to {stderr} if {verbose} is true.  The procedures abort on
    any errors. */
     
char *read_line(FILE *f);
  /* Reads the next line from file {f}, until '\n', NUL, or EOF.
    Returns it as a NUL-terminated string, without the final '\n'.
    Replaces TABs by single spaces. Allocates space for the result with 
    {malloc()}. Returns NULL iff the file is already at EOF. 
    
    !!! Should receive a max-length argument for safety !!! */
 
#endif
