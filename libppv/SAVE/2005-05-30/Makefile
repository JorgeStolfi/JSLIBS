# Last edited on 2007-04-16 21:09:02 by stolfi

LIBNAME = libppv

HFILES := \
  ppv_types.h \
  ppv_array.h \
  ppv_io.h
  
NON_HFILES := \
  ppv_packing.h \
  ppv_convert.h
  
CFILES := \
  ppv_array.c \
  ppv_io.c
  
NON_CFILES := \
  ppv.c
  
LIBFILE = \
  ${LIBNAME}.a
  
.PHONY: all clean uninstall install check

# all: clean uninstall install
all: uninstall install
# all: check

HOFILES := ${subst .h,.ho,${HFILES}}
OFILES := ${subst .c,.o,${CFILES}}

uninstall:
	cd ${INC} && /bin/rm -f ${HFILES} ${HOFILES}
	cd ${LIB} && /bin/rm -f ${LIBFILE}

install: ${HOFILES} ${LIBFILE}
	cp -p ${HFILES} ${INC}
	cp ${HOFILES} ${INC}
	cp -p ${LIBFILE} ${LIB}
  
check:
	cd tests && make all

clean:
	-/bin/rm -f *.o *.ho *.a core
	cd tests && make clean

# ======================================================================

INC = ${HOME}/PUB/include
LIB = ${HOME}/PUB/${PLATFORM}/lib
  
IFLAGS := \
  -I. -I${INC}
  
CFLAGS := \
  -g -ffloat-store -frounding-math \
  -Wall -Wpointer-arith -Wmissing-prototypes
  
%.o: %.c ;\
  /usr/bin/gcc -c ${CFLAGS} ${IFLAGS} $*.c
  
%.ho: %.h ;\
  /usr/bin/gcc -o $*.ho -c ${CFLAGS} ${IFLAGS} -x c $*.h \
  || /bin/rm -f $*.ho
  
${LIBFILE}: ${OFILES} ;\
  rm -f $*.a && ar crv $*.a ${OFILES} && ranlib $*.a

# ====================================================================
# Dependencies:

DEPFILE := Deps.make

depend: 
	/bin/rm -f ${DEPFILE}
	extract-ho-deps ${IFLAGS} -I/usr/include ${HFILES} ${CFILES} \
          | egrep -v ': /usr/include' \
          > ${DEPFILE}

# Include specific dependencies extracted by "make depend"
include ${DEPFILE}
