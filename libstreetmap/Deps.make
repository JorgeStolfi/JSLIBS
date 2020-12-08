stheap.ho: stheap.h

stimage.ho: stimage.h
stimage.ho: stmap.ho

stmap.ho: stmap.h
stmap.ho: /home/jstolfi/include/quad.ho
stmap.ho: /home/jstolfi/include/r2.ho
stmap.ho: /home/jstolfi/include/vec.ho
stmap.ho: /home/jstolfi/include/pswr.ho

stheap.o: stheap.c
stheap.o: stheap.ho
stheap.o: /home/jstolfi/include/affirm.ho

stimage.o: stimage.c
stimage.o: /home/jstolfi/include/fget.ho
stimage.o: /home/jstolfi/include/nget.ho
stimage.o: /home/jstolfi/include/affirm.ho
