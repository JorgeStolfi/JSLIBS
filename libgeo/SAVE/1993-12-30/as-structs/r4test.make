INC = ${HOME}/PUB/gcc/include
LIB = ${HOME}/PUB/gcc/lib

HFILES =

HOFILES =

OFILES =
  
PROG = \
  r4test
  
MAINOBJ = \
  $(PROG).o

LIBS = \
  $(LIB)/libgeo.a \
  $(LIB)/libjs.a

GCCFLAGS = \
  -I$(INC) \
  -pg -g \
  -ansi \
  -Wall -Wtraditional -Wpointer-arith -Wmissing-prototypes \
  -mfpu -ffloat-store -frounding-math

all: cleanup $(HOFILES) $(PROG) install

cleanup: ;\
  rm -f $(PROG)

install:

%.o: %.c ;\
  gcc -c $(GCCFLAGS) $*.c
  
%.ho: %.h ;\
  gcc -o $*.ho -c $(GCCFLAGS) -x c $*.h \
  || /bin/rm -f $*.ho
  
$(PROG): $(OFILES) $(MAINOBJ) $(LIBS) ;\
  gcc -pg -o $(PROG) $(MAINOBJ) $(OFILES) $(LIBS) -lm &&\
  $(PROG)

# Dependencies of .h files:
  
# Dependencies of .c files:

r4test.o:: r4test.c $(INC)/r4.ho $(INC)/r4x4.ho \
           $(INC)/js.ho $(INC)/ioprotos.ho

