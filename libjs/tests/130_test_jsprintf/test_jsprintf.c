#define PROG_NAME "test_jsprintf"
#define PROG_DESC "test of {jsprintf.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 04:02:22 by stolfi */
/* Created on 2024-11-18 by J. Stolfi, UNICAMP */

#define test_jsprintf_COPYRIGHT \
  "Copyright © 2024  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>

#define SBUF_SIZE 1000*1000
  /* Buffer size for {snprintf} */

#define do_test(FMT,...) \
  do \
    { \
      fprintf(stderr, "fmt = |%s|\n", FMT); \
      sres_len = snprintf(sres, SBUF_SIZE, FMT, ##__VA_ARGS__); \
      fprintf(stderr, "slen = %d\n", sres_len); \
      demand(sres_len >= 0, "{snprintf} failed"); \
      demand(sres_len < SBUF_SIZE, "formatted result is way too long"); \
      fprintf(stderr, "sres = |%s|\n", sres); \
      jres = jsprintf(FMT, ##__VA_ARGS__); \
      jres_len = (int32_t)strlen(jres); \
      fprintf(stderr, "jres = |%s|\n", jres); \
      demand(sres_len == jres_len, "lengths differ"); \
      demand(strcmp(sres, jres) == 0, "formatted results differ");  \
      free(jres);  \
      fprintf(stderr, "\n"); \
    } \
  while(FALSE)
  /* Prints the given {...} values with format {FMT}
    using {snprintf} and {jsprintf}, and compares the results.
    Does NOT check whether the allocated area was indeed
    big enough to contain the result string. */

int32_t main(int32_t argc, char **argv);

/* MAIN */

int32_t main(int32_t argc, char **argv)
  {
    char *sres = talloc(SBUF_SIZE, char);
    int32_t sres_len;
    
    char *jres = NULL;
    int32_t jres_len;

    do_test("");
    do_test("!");
    do_test("hello world");
    do_test("¿ sabão × açúcar ?");
    do_test("char %c", '?');
    do_test("string %12s", "¿¿¿¿¿¿¿¿¿¿¿¿");
    do_test("integer %d", -12345);
    do_test("integer %ld", -12345L);
    do_test("integer %lld", -12345LL);
    do_test("integer %*d", 30, -12345);
    do_test("unsigned %u", 54321);
    do_test("unsigned %lu", 54321LU);
    do_test("unsigned %llu", 54321LLU);
    do_test("double %.15f", M_PI);
    do_test("double %*.*f", 40, 15, M_PI);
    do_test("double %f", DBL_MAX);
    fprintf(stderr, "\n=== done ===\n");
    return 0;
  }
