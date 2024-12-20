# Last edited on 2024-12-04 19:34:51 by stolfi

.PHONY:: run-f run-c diffs 

all: check diffs

build: build-f build-c

check: run-f run-c 

actions: ${ACTIONS}

# Output files of the C and Fortran versions of the code:

build-f:
	cd f && make build

run-f:
	cd f && rm -fv out/*.txt
	cd f && make check

build-c:
	cd c && make build

run-c:
	cd c && rm -fv out/*.txt
	cd c && make check

# DAMN MINUS ZEROS!!!!
diffs: f/.log c/.log f/out/test1.txt c/out/test1.txt f/out/test2.txt c/out/test2.txt
	@echo "=== f/.log c/.log ===" 1>&2
	prdiff -bB f/.log c/.log
	@for fn in test1 test2 ; do \
          echo "=== f/out/$${fn}.txt c/out/$${fn}.txt ===" ; \
          cat f/out/$${fn}.txt | sed -e 's:[-]\(0[.]000*\b\): \1:g' > .tf ; \
          cat c/out/$${fn}.txt | sed -e 's:[-]\(0[.]000*\b\): \1:g' > .tc ; \
          prdiff -bB .tf .tc ; \
        done

clean:: 
	for dir in f c ; do \
          ( cd $$dir && make clean ) ; \
        done
