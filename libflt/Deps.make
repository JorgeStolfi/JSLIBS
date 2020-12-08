flt.ho: flt.h

flteval.ho: flteval.h
flteval.ho: flt.ho
flteval.ho: pcode.ho

pcode.ho: pcode.h

flt.o: flt.c
flt.o: /home/jstolfi/include/affirm.ho
flt.o: /home/jstolfi/include/jsrandom.ho
flt.o: flt.ho

fltasm.o: fltasm.c
fltasm.o: flt.ho

flteval.o: flteval.c
flteval.o: flteval.ho
flteval.o: flt.ho
flteval.o: pcode.ho
flteval.o: /home/jstolfi/include/affirm.ho

pcode.o: pcode.c
pcode.o: pcode.ho
pcode.o: /home/jstolfi/include/affirm.ho
pcode.o: /home/jstolfi/include/jsfile.ho

