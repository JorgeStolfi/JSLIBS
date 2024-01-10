
OBJS =  \
  foi1g.o \
  foi2z.o \
  foi2q.o \
  foi2qe.o \
  fgraph.o \
  zeros2.o
  
GCCFLAGS = \
  -g \
  -ansi \
  -Wall -Wtraditional -Wpointer-arith -Wmissing-prototypes \
  -mfpu -ffloat-store -frounding-math

all: libfoi.a foitest

foitest: $(OBJS) libfoi.a foitest.o ;\
  gcc -o foitest foitest.o $(OBJS) -lfoi -lm
  
libfoi.a:: ;\
  $(MAKE) -f libfoi.make

%.o: %.c ;\
  gcc -c $(GCCFLAGS) $*.c
  
%.ho: %.h ;\
  gcc -o $*.ho -c $(GCCFLAGS) -x c $*.h || /bin/rm -f $*.ho
  
/* Dependencies of .h files: */
  
fgraph.ho::   foifloat.ho interval.ho

foi1g.ho::    foi.ho foifloat.ho interval.ho

foi2q.ho::    foi.ho foifloat.ho interval.ho

foi2qe.ho::   foi.ho foifloat.ho interval.ho

foi2z.ho::    foi.ho foifloat.ho interval.ho

zeros2.ho::   foifloat.ho interval.ho

# Dependencies for .c files:

fgraph.o::   fgraph.ho foifloat.ho foimisc.ho iomisc.ho interval.ho plt0.ho

foi1g.o::    foi1g.ho foi.ho foifloat.ho foimisc.ho interval.ho iomisc.ho \
             fgraph.ho plt0.ho 

foi2q.o::    foi2q.ho foi.ho foifloat.ho foimisc.ho interval.ho iomisc.ho \
             zeros2.ho plt0.ho 

foi2qe.o::   foi2qe.ho foi.ho foifloat.ho foimisc.ho interval.ho iomisc.ho \
             zeros2.ho plte.ho 

foi2z.o::    foi2z.ho foi.ho foifloat.ho foimisc.ho interval.ho iomisc.ho \
             zeros2.ho plt0.ho 

foitest.o::  foi.ho foifloat.ho foimisc.ho interval.ho iomisc.ho \
             foi1g.ho foi2z.ho foi2q.ho foi2qe.o \
             f0.c f1.c f2.c f3.c f4.c f5.c f6.c f7.c \
             g0.c g1.c g2.c g3.c g4.c g5.c g6.c g7.c g8.c g9.c

zeros2.o::   zeros2.ho foifloat.ho foimisc.ho iomisc.ho interval.ho plt0.ho
