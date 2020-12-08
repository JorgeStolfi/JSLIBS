stpoly.ho: stpoly.h
stpoly.ho: /home/jstolfi/include/vec.ho
stpoly.ho: /home/jstolfi/include/bool.ho
stpoly.ho: /home/jstolfi/include/i2.ho
stpoly.ho: /home/jstolfi/include/r2.ho
stpoly.ho: stpoly_STP.ho

stpoly_STP.ho: stpoly_STP.h
stpoly_STP.ho: /home/jstolfi/include/vec.ho
stpoly_STP.ho: /home/jstolfi/include/bool.ho
stpoly_STP.ho: /home/jstolfi/include/r2.ho
stpoly_STP.ho: /home/jstolfi/include/i2.ho

stpoly_rep.ho: stpoly_rep.h
stpoly_rep.ho: /home/jstolfi/include/vec.ho
stpoly_rep.ho: /home/jstolfi/include/bool.ho
stpoly_rep.ho: /home/jstolfi/include/i2.ho
stpoly_rep.ho: /home/jstolfi/include/r2.ho
stpoly_rep.ho: stpoly_STP.ho
stpoly_rep.ho: stpoly.ho

stpoly.o: stpoly.c
stpoly.o: /home/jstolfi/include/bool.ho
stpoly.o: /home/jstolfi/include/affirm.ho
stpoly.o: /home/jstolfi/include/i2.ho
stpoly.o: /home/jstolfi/include/r2.ho
stpoly.o: /home/jstolfi/include/bvtable.ho
stpoly.o: /home/jstolfi/include/bvhash.ho
stpoly.o: stpoly_STP.ho
stpoly.o: stpoly.ho

stpoly_STP.o: stpoly_STP.c
stpoly_STP.o: /home/jstolfi/include/bool.ho
stpoly_STP.o: /home/jstolfi/include/argparser.ho
stpoly_STP.o: /home/jstolfi/include/fget.ho
stpoly_STP.o: /home/jstolfi/include/jsmath.ho
stpoly_STP.o: stpoly_STP.ho

stpoly_rep.o: stpoly_rep.c
stpoly_rep.o: /home/jstolfi/include/bool.ho
stpoly_rep.o: /home/jstolfi/include/affirm.ho
stpoly_rep.o: /home/jstolfi/include/i2.ho
stpoly_rep.o: /home/jstolfi/include/r2.ho
stpoly_rep.o: stpoly.ho
stpoly_rep.o: stpoly_rep.ho

