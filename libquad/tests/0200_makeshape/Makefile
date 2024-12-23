# Last edited on 2024-12-22 11:34:13 by stolfi

PROG := test_makeshape
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

# Format: ${NAME}

TESTS := \
  torus \
  bitorus \
  tritorus \
  klein \
  klein2 \
  klein3 \
  projective \
  tetra \
  stick \
  ring \
  cube \
  sausage \
  orange \
  fork \
  star \

JS_LIBS := \
  libjs.a

include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make
 
all: check

check:  
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t do-test ; \
        done
        
######################################################################
# Section for recursive make
# Client must define ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst :, ,${TEST}}
NAME :=    ${word 1, ${TEST_ARGS}}

do-test: ${PROG} out
	${PROG} -shape ${NAME} > out/${NAME}.oct

endif
# End of ${TEST} section
######################################################################
