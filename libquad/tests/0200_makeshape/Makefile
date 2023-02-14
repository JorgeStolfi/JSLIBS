# Last edited on 2009-03-08 00:09:16 by stolfi

PROG := makeshape
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

# Format: ${REFINE}:${NAME}

TESTS := \
  01:torus \
  01:bitorus \
  01:tritorus \
  01:klein \
  01:klein2 \
  01:klein3 \
  01:projective \
  01:tetra \
  01:stick \
  01:ring \
  01:cube \
  01:sausage \
  01:orange \
  01:fork \
  01:star \

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
REFINE :=  ${word 1, ${TEST_ARGS}}
NAME :=    ${word 2, ${TEST_ARGS}}

do-test: ${PROG} out
	${PROG} -shape ${NAME} -refine ${REFINE} > out/${NAME}.oct

endif
# End of ${TEST} section
######################################################################
