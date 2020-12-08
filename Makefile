# Makefile for C libraries
# Last edited on 2020-12-08 13:50:22 by jstolfi

# ----------------------------------------------------------------------
# Good libraries, in import-compatible order:

LIBS_DIR := 00-LIBS-STATUS.dir

LIBS_OK :=         ${shell gawk '/^[+] /   { print $$2; }' ${LIBS_DIR}}
LIBS_UNCERTAIN :=  ${shell gawk '/^[?] /   { print $$2; }' ${LIBS_DIR}}
LIBS_BUGGY :=      ${shell gawk '/^[!] /   { print $$2; }' ${LIBS_DIR}}
LIBS_JUNK :=       ${shell gawk '/^[J] /   { print $$2; }' ${LIBS_DIR}}

LIBS_ALL :=        ${LIBS_OK} ${LIBS_UNCERTAIN} ${LIBS_BUGGY}
LIBS_NOT_OK :=     ${LIBS_UNCERTAIN} ${LIBS_BUGGY}
LIBS_NOT_BUGGY :=  ${LIBS_OK} ${LIBS_UNCERTAIN} 

# ----------------------------------------------------------------------
# Library tests to run:

LIB_TESTS_DIR := 00-LIB-TESTS-STATUS.dir

LIB_TESTS_OK :=         ${shell gawk '/^[+] /   { print $$2; }' ${LIB_TESTS_DIR}}
LIB_TESTS_UNCERTAIN :=  ${shell gawk '/^[?] /   { print $$2; }' ${LIB_TESTS_DIR}}
LIB_TESTS_BUGGY :=      ${shell gawk '/^[!] /   { print $$2; }' ${LIB_TESTS_DIR}}
LIB_TESTS_JUNK :=       ${shell gawk '/^[J] /   { print $$2; }' ${LIB_TESTS_DIR}}

LIB_TESTS_ALL :=        ${LIB_TESTS_OK} ${LIB_TESTS_UNCERTAIN} ${LIB_TESTS_BUGGY}
LIB_TESTS_NOT_OK :=     ${LIB_TESTS_UNCERTAIN} ${LIB_TESTS_BUGGY}
LIB_TESTS_NOT_BUGGY :=  ${LIB_TESTS_OK} ${LIB_TESTS_UNCERTAIN} 

# ----------------------------------------------------------------------
# Libraries to process:

# LIBS := ${LIBS_UNCERTAIN}
# LIBS := ${LIBS_NOT_OK}
# LIBS := ${LIBS_ALL}
LIBS := ${LIBS_NOT_BUGGY}
 
# ----------------------------------------------------------------------
# Library tests to run:

# LIB_TESTS := ${LIB_TESTS_UNCERTAIN}
# LIB_TESTS := ${LIB_TESTS_NOT_OK}
LIB_TESTS := ${LIB_TESTS_ALL}

# ----------------------------------------------------------------------
# Actions to perform:

# all:  all-clean all-build all-check
all:  all-clean all-build
# all:  all-check
 
all-clean:
	@for dir in ${LIBS_ALL}; do \
          ( echo '= = all-clean = = = = = = = = = = = = = = = = = = = = = = = = = = = ='; \
            echo "$$dir"; \
            cd $$dir/. && \
            make uninstall clean ; \
          ) ; \
        done

all-build:
	@for dir in ${LIBS}; do \
          ( echo '= = all-build = = = = = = = = = = = = = = = = = = = = = = = = = = = = ='; \
            echo "$$dir"; \
            cd $$dir/. && \
            make depend build-libs install ; \
          ) ; \
        done

all-check:
	@for dir in ${LIB_TESTS}; do \
          ( echo '= = all-check = = = = = = = = = = = = = = = = = = = = = = = = = = = ='; \
            echo "$$dir"; \
            cd $$dir/. && \
            make depend build-tests check ; \
          ) ; \
        done
