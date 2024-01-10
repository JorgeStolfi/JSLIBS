
OBJS =  \
  foifloat.o \
  foimisc.o \
  interval.o \
  foi.o \
  plt0.o \
  plte.o \
  pstools.o
  
GCCFLAGS = \
  -g \
  -ansi -ffloat-store -frounding-math \
  -Wall -Wtraditional -Wpointer-arith -Wmissing-prototypes \
  -mfpu

all: libfoi.a

libfoi.a: $(OBJS) ;\
  gcc -c $(OBJS) && ar crv $*.a $(OBJS) && ranlib $*.a
  
%.o: %.c ;\
  gcc -c $(GCCFLAGS) $*.c
  
%.ho: %.h ;\
  gcc -o $*.ho -c $(GCCFLAGS) -x c $*.h || /bin/rm -f $*.ho
  
# Dependencies of .h files: 
  
foi.ho::      foifloat.ho foimisc.ho interval.ho

interval.ho:: foifloat.ho

# Dependencies for .c files:

foi.o::      foi.ho foimisc.ho foifloat.ho interval.ho iomisc.ho

foifloat.o:: foifloat.ho foimisc.o iomisc.ho

foimisc.o::  foimisc.ho iomisc.ho

interval.o:: interval.ho foifloat.ho foimisc.ho iomisc.ho

plt0.o::     plt0.ho pstools.ho foimisc.ho iomisc.ho

plte.o::     plte.ho pstools.ho foimisc.ho iomisc.ho

pstools.o::  pstools.ho foimisc.ho iomisc.ho

