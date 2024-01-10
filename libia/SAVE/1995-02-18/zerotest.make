# Last edited on 2007-04-16 21:10:26 by stolfi

PROG = zerotest
  
INC = ${HOME}/PUB/include
LIB = ${HOME}/PUB/${PLATFORM}/lib

HFILES =

HOFILES =

OFILES =
  
MAINOBJ = \
  zerotest.o
  
LIBS = \
  ${LIB}/libia.a \
  ${LIB}/libflt.a \
  ${LIB}/libjs.a

OPTFLAGS := -g

IFLAGS := \
  -I${INC}

CFLAGS := \
  ${OPTFLAGS} \
  -ansi -ffloat-store -frounding-math \
  -Wall  -Wpointer-arith -Wmissing-prototypes \
  

all: cleanup ${HOFILES} ${PROG} run

cleanup: 

install: 

%.o: %.c ;\
  ${CC} -c ${CFLAGS} ${IFLAGS} $*.c
  
%.ho: %.h ;\
  ${CC} -o $*.ho -c ${CFLAGS} ${IFLAGS} -x c $*.h \
  || /bin/rm -f $*.ho
  
${PROG}: ${PROG}.o ${OFILES} ${LIBS} ;\
  rm -f ${PROG} && \
  ${CC} ${LDFLAGS} -o ${PROG} ${PROG}.o ${OFILES} -L${LIB} ${LIBS} -lrt -lm
  
run: ${PROG}
	${PROG} 

# ====================================================================
# Dependencies:

DEPFILE := Deps.make

depend: 
	/bin/rm -f ${DEPFILE}
	extract-ho-deps ${IFLAGS} -I/usr/include ${HFILES} ${CFILES} ${PROG}.c \
          | egrep -v ': /usr/include' \
          > ${DEPFILE}

# Include specific dependencies extracted by "make depend"
include ${DEPFILE}
