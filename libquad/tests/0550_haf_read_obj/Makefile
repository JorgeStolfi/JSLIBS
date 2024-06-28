# Last edited on 2024-06-27 10:45:47 by stolfi

PROG := test_haf_read_obj
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

# Format: ${OBJNAME}

TESTS := \
  holey
  
OLTHER_TESTS := \
  holey \
  pinch \
  cubic \
  tetra

JS_LIBS := \
  libgeo.a \
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
OBJNAME :=   ${word 1, ${TEST_ARGS}}

IN_OBJ_FILE := in/${OBJNAME}.obj
OT_OBJ_FILE := out/${OBJNAME}-wr.obj

do-test: ${PROG} out
	${PROG} ${OBJNAME}

endif
# End of ${TEST} section
######################################################################
