# Last edited on 2020-10-11 15:23:35 by jstolfi

PROG := test_nnclassif

TEST_LIB := libclassif.a
TEST_LIB_DIR := ../..

JS_LIBS = \
  libimg.a \
  libgeo.a \
  libjs.a

VERIFY := YES
SHOW_SIZE := 512x

# PROBLEM:NA_R:NA_I:NC:NS_MODEL:NS_REFINE:NS_EVAL:GRID:SEED:NOISE:OPF:IMGSIZE

# Tests of throwing and classifying with each problem.
# No refinement or OPF handicaps.
PROBLEM_TESTS := \
  saturn:02:00:02:000100:000000:000000:Y:003173:0.00:N:512 \
  petals:02:00:04:000100:000000:000000:Y:003173:0.00:N:512 \
  vessel:02:00:03:000100:000000:000000:Y:003173:0.00:N:512 \
  mballs:02:00:04:000100:000000:000000:Y:003173:0.00:N:512 \
  shells:02:00:03:000100:000000:000000:Y:003173:0.00:N:512 \
  staoli:02:00:03:000100:000000:000000:Y:003173:0.00:N:512 \

# Basic tests of all program stages with all problems.
# Include refinement but not OPF handicaps.
BASIC_PROG_TESTS := \
  saturn:02:00:02:000100:000300:001000:Y:003173:0.00:N:256 \
  petals:02:00:04:000100:000300:001000:Y:003173:0.00:N:256 \
  vessel:02:00:03:000100:000300:001000:Y:003173:0.00:N:256 \
  mballs:02:00:04:000100:000300:001000:Y:003173:0.00:N:256 \
  shells:02:00:03:000100:000300:001000:Y:003173:0.00:N:256 \
  staoli:02:00:03:000100:000300:001000:Y:003173:0.00:N:256 \

# Tests to compare the accuracy of the basic 1NN classifier
# with varying model sizes.
NN_ACCURACY_TESTS := \
  shells:02:00:02:000010:000000:001000:Y:003173:0.00:N:256 \
  shells:02:00:02:000100:000000:001000:Y:003173:0.00:N:256 \
  shells:02:00:02:001000:000000:001000:Y:003173:0.00:N:256 \
  \
  mballs:02:00:03:000010:000000:001000:Y:003173:0.00:N:256 \
  mballs:02:00:03:000100:000000:001000:Y:003173:0.00:N:256 \
  mballs:02:00:03:001000:000000:001000:Y:003173:0.00:N:256
  
OPF_HANDICAP_TESTS_1 := \
  mballs:02:00:04:000200:000300:010000:Y:003173:0.00:N:256 \
  mballs:02:00:04:000200:000300:010000:Y:003173:0.00:Y:256
  
OPF_HANDICAP_TESTS_2 := \
  shells:02:00:02:000100:000300:010000:Y:003173:0.00:Y:256 \
  shells:02:00:02:000100:000300:010000:Y:003173:0.00:N:256
  
OPF_HANDICAP_TESTS_3 := \
  staoli:02:00:03:000100:000300:010000:Y:003173:0.00:N:256 \
  staoli:02:00:03:000100:000300:010000:Y:003173:0.00:Y:256
   
OPF_HANDICAP_TESTS_4_NOISY := \
  mballs:02:00:02:000100:000300:010000:Y:003173:0.02:N:256 \
  mballs:02:00:02:000100:000300:010000:Y:003173:0.02:Y:256
   
# Attempt to reproduce tables IV and V of tht
# IJIST paper (Papa, Falcao, Suzuki 2009).
# The |E| sets are random rather than grid because
# the domains are very small.  We use 100 times 
# larger E-sets to get more accurate figures.
OPF_HANDICAP_TESTS_IJIST := \
  saturn:02:00:02:000067:000033:010000:N:419417:0.00:Y:256 \
  petals:02:00:04:000033:000017:005000:N:419417:0.00:Y:256 \
  vessel:02:00:03:000033:000017:005000:N:419417:0.00:Y:256 \
  \
  saturn:02:00:02:000067:000033:010000:N:419417:0.00:N:256 \
  petals:02:00:04:000033:000017:005000:N:419417:0.00:N:256 \
  vessel:02:00:03:000033:000017:005000:N:419417:0.00:N:256 \
 
# Tests of effects of irrelevant attributes.
IRRELEVANT_ATTR_TESTS := \
  mballs:02:00:02:000100:000300:010000:Y:003173:0.00:N:256 \
  mballs:02:01:02:000100:000300:010000:Y:003173:0.00:N:256 \
  mballs:02:01:02:000100:000300:010000:Y:003173:0.00:Y:256
 
OTHER_TESTS := \
  vessel:02:00:03:000100:000300:010000:Y:003173:0.02:N:256 \
  \
  mballs:02:00:03:000100:000300:010000:Y:003173:0.00:N:256 \
  \
  mballs:02:00:03:000100:000300:010000:Y:003173:0.00:N:256 \
  \
  mballs:02:00:03:001000:003000:010000:Y:003173:0.00:N:256 \
  mballs:02:00:03:001000:003000:010000:Y:003173:0.00:N:256 \
  shells:02:00:04:000100:000300:010000:Y:003173:0.00:N:256 \
  shells:02:00:04:000100:000300:010000:Y:003173:0.00:N:256
 
# TESTS := ${PROBLEM_TESTS}
# TESTS := ${BASIC_PROG_TESTS}
# TESTS := ${OPF_HANDICAP_TESTS_1}
# TESTS := ${OPF_HANDICAP_TESTS_2}
# TESTS := ${OPF_HANDICAP_TESTS_3}
# TESTS := ${OPF_HANDICAP_TESTS_4_NOISY}
TESTS := ${OPF_HANDICAP_TESTS_IJIST}
# TESTS := ${NN_ACCURACY_TESTS}
# TESTS := ${IRRELEVANT_ATTR_TESTS}

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

check: 
	for t in ${TESTS}; do \
          ${MAKE} TEST=$$t single ; \
        done
        
clean::
	( cd out && rm -f {*/,}*.{dat,ppm,eps} )
        
ifneq "/${TEST}" "/"        
########################################################################
# Make with ${TEST}

TEST_FIELDS := ${subst :, ,${TEST}}
PROBLEM :=     ${word  1, ${TEST_FIELDS}}
NA_R :=        ${word  2, ${TEST_FIELDS}}
NA_I :=        ${word  3, ${TEST_FIELDS}}
NC :=          ${word  4, ${TEST_FIELDS}}
NS_MODEL :=    ${word  5, ${TEST_FIELDS}}
NS_REFINE :=   ${word  6, ${TEST_FIELDS}}
NS_EVAL :=     ${word  7, ${TEST_FIELDS}}
GRID :=        ${word  8, ${TEST_FIELDS}}
SEED :=        ${word  9, ${TEST_FIELDS}}
NOISE :=       ${word 10, ${TEST_FIELDS}}
OPF :=         ${word 11, ${TEST_FIELDS}}
IMAGE_SIZE :=  ${word 12, ${TEST_FIELDS}}

NAME :=  ${subst .,v,${subst :,-,${TEST}}}
OUTDIR := out/${NAME}
PREFIX := ${OUTDIR}/t

single: ${PROG}
	mkdir -p ${OUTDIR}
	${PROG} \
            -problem ${PROBLEM} \
            -samples ${NS_MODEL} ${NS_REFINE} ${NS_EVAL} \
            -seed ${SEED} \
            -prefix ${PREFIX} \
            -attributes ${NA_R} ${NA_I} \
            -classes ${NC} \
            -noise ${NOISE} \
            -opf ${OPF} \
            -grid ${GRID} \
            -verify ${VERIFY} \
            -image ${IMAGE_SIZE} 2
	-display -title '%f' -filter point -resize '${SHOW_SIZE}' ${PREFIX}-*.ppm
        
        
# end ${TEST}
########################################################################endif
endif
