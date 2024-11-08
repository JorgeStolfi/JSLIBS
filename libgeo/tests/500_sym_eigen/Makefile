# Last edited on 2024-11-06 04:32:32 by stolfi

.PHONY:: diffs f-out c-out

all: check 

check:  diffs

actions: ${ACTIONS}

# Output files of the C and Fortran versions of the code:
C_OUTPUT := out/testc.txt
F_OUTPUT := out/testf.txt

diffs: c/${C_OUTPUT} f/${F_OUTPUT}
	prdiff -bB f/${F_OUTPUT} c/${C_OUTPUT}

f/${F_OUTPUT}: f-out

f-out:
	cd f && make ${F_OUTPUT}

c/${C_OUTPUT}: c-out

c-out:
	cd c && make ${C_OUTPUT}

clean:: 
	for dir in f c ; do \
          ( cd $$dir && make clean ) ; \
        done
