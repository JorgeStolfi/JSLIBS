# Last edited on 2023-02-15 12:41:49 by stolfi

.PHONY:: diffs

all: check 

check:  diffs

actions: ${ACTIONS}

# Output files of the C and Fortran versions of the code:
C_OUTPUT := out/testc.txt
F_OUTPUT := out/testf.txt

diffs: c/${C_OUTPUT} f/${F_OUTPUT}
	diff -b f/${F_OUTPUT} c/${C_OUTPUT} \
          | prettify-diff-output

f/${F_OUTPUT}: 
	cd f && make ${F_OUTPUT}

c/${C_OUTPUT}: 
	cd c && make ${C_OUTPUT}

clean:: 
	for dir in f c ; do \
          ( cd $$dir && make clean ) ; \
        done
