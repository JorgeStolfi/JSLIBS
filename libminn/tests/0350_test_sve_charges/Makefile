# Last edited on 2024-11-07 18:59:14 by stolfi

PROG = test_sve_charges
 
TEST_LIB := libminn.a
TEST_LIB_DIR := ../..

JS_LIBS := \
  libimg.a \
  libgeo.a \
  libjs.a
  
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

# Test run arguments: {NQ}:{DIM}:{SYM}:{LMAX}
# 
# where {SYM} is the symmetry class, {LMAX} is the relative maximum 
# electron pair connector to be added to the 3D model.
 
TESTS := \
  060:3:0:0.500 \
  060:3:1:0.500

NON_TESTS := \
  001:2:0:2.100 \
   \
  002:1:0:2.100 \
   \
  003:2:0:1.500 \
  003:2:1:1.500 \
  003:2:2:1.500 \
   \
  004:3:0:1.500 \
  004:3:2:1.500 \
   \
  005:3:0:1.300 \
  005:3:1:1.300 \
  005:3:3:1.300 \
   \
  006:3:0:1.300 \
  006:3:1:1.300 \
  006:3:2:1.300 \
   \
  007:3:0:1.200 \
  007:3:1:1.200 \
  007:3:3:1.200 \
   \
  008:3:0:1.150 \
  008:3:1:1.150 \
  008:3:2:1.150 \
   \
  009:3:0:1.150 \
  009:3:1:1.150 \
  009:3:3:1.150 \
   \
  010:3:0:1.150 \
  010:3:1:1.100 \
   \
  011:3:0:1.100 \
  011:3:1:1.100 \
   \
  012:3:0:1.150 \
  012:3:1:1.150 \
  012:3:2:1.150 \
   \
  013:3:0:1.150 \
  013:3:1:1.150 \
  013:3:3:1.150 \
   \
  014:3:0:1.150 \
  014:3:1:1.000 \
   \
  015:3:0:0.900 \
  015:3:1:0.900 \
   \
  016:3:0:0.900 \
  016:3:1:0.900 \
   \
  017:3:0:0.900 \
  017:3:1:0.900 \
   \
  018:3:0:0.900 \
  018:3:1:0.900 \
   \
  019:3:0:0.900 \
  019:3:1:0.900 \
   \
  020:3:0:0.800 \
  020:3:1:0.800 \
  020:3:2:0.800 \
   \
  021:3:0:0.800 \
  021:3:1:0.800 \
  021:3:3:0.800 \
   \
  022:3:0:0.800 \
  022:3:1:0.800 \
   \
  023:3:0:0.800 \
  023:3:1:0.800 \
   \
  024:3:0:0.800 \
  024:3:1:0.800 \
   \
  030:3:0:0.700 \
  030:3:1:0.700 \
   \
  031:3:0:0.700 \
  031:3:1:0.700 \
   \
  032:3:0:0.700 \
  032:3:1:0.700 \
   \
  033:3:0:0.700 \
  033:3:1:0.700 \
   \
  036:3:0:0.650 \
  036:3:1:0.650 \
   \
  060:3:0:0.500 \
  060:3:1:0.500

.PHONY:: all-tests single-test

all: check

check:  all-tests

all-tests: ${PROG}
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt single-test ; \
        done
	
########################################################################
# Recursive make -- caller must define ${TEST}

TEST = TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

TESTARGS := ${subst :, ,${TEST}}
NQ :=      ${word 1, ${TESTARGS}}
DIM :=     ${word 2, ${TESTARGS}}
SYM :=     ${word 3, ${TESTARGS}}
LMAX := ${word 4, ${TESTARGS}}

NAME := test-${NQ}-${DIM}-${SYM}

SOLFILE := out/${NAME}-sol.dat
POTFILE := out/${NAME}-pot.dat

single-test: ${PROG} 
	time ${PROG} ${NQ} ${DIM} ${SYM}
	if [[ -r ${POTFILE} ]] ; then \
          show-potential-2d.sh ${NQ} ${POTFILE} & \
        fi
	if [[ -r ${SOLFILE} ]] ; then \
          show-solution-3d.sh ${NQ} ${LMAX} ${SOLFILE} & \
        fi
	@-xmessage "OK"

endif
# end ${TEST} make
########################################################################

clean::
	-/bin/rm -f out/*.dat
