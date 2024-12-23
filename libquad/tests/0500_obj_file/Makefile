# Last edited on 2024-12-22 11:24:57 by stolfi

PROG := test_obj_file
 
TEST_LIB := libquad.a

TEST_LIB_DIR := ../..

# Format: ${OBJNAME}

TESTS := \
  holey
  
OLTHER_TESTS := \
  tetra \
  cubic

JS_LIBS := \
  libgeo.a \
  libgausol.a \
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
