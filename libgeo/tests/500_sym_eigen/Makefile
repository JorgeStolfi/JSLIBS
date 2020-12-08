# Last edited on 2019-04-09 19:41:05 by jstolfi

.PHONY:: diffs

all: check

check:  diffs

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
