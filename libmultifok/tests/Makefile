# Last edited on 2023-04-28 18:36:40 by stolfi

GOOD_TESTS := \
  0100_make_image \
  0150_show_basis \
  0180_diff_ops
  
JUNK_TESTS :=
  
LATER_TESTS := \
  0300_merge_images \
  0400_sort_stack \

IGNOREDIRS := ${LATER_TESTS}  
  
include ${STOLFIHOME}/programs/c/GENERIC-ROOT-DIR.make

all: check

check:  build-tests check
