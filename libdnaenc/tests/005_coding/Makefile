# Last edited on 2020-10-11 15:24:04 by jstolfi

TEST_LIB  := libdnaenc.a
TEST_LIB_DIR  := ../..

PROG := test_dnae_coding

JS_LIBS := \
  libjs.a
 
include ${STOLFIHOME}/programs/c/GENERIC-LIB-TEST.make

all: check

# check:  clean run show
check:  clean run

# Average samples per histogram bin:
NSAMPLES := 200

OUTDIR := out
HISTNAME := ${OUTDIR}/hist

clean:: 
	-/bin/rm -f ${HISTNAME}-*.pgm
  
run: ${HISTNAME}.txt ${HISTNAME}-e.png ${HISTNAME}-p.png ${HISTNAME}-v.png

${HISTNAME}.txt: ${PROG} 
	${PROG} ${NSAMPLES} > ${HISTNAME}.txt
        
${HISTNAME}-e.png: ${HISTNAME}.txt plot_errors.sh     
	plot_errors.sh ${HISTNAME}
        
${HISTNAME}-p.png: ${HISTNAME}.txt plot_probs.sh     
	plot_probs.sh ${HISTNAME}
        
${HISTNAME}-v.png: ${HISTNAME}.txt plot_values.sh     
	plot_values.sh ${HISTNAME}

show: 
	@echo "** show what?"
