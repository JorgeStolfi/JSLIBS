oct.ho: oct.h
oct.ho: /home/jstolfi/include/vec.ho
oct.ho: /home/jstolfi/include/bool.ho

oct_enum.ho: oct_enum.h
oct_enum.ho: oct.ho

oct_shapes.ho: oct_shapes.h
oct_shapes.ho: oct.ho
oct_shapes.ho: /home/jstolfi/include/bool.ho

quad.ho: quad.h
quad.ho: /home/jstolfi/include/vec.ho

oct.o: oct.c
oct.o: /home/jstolfi/include/vec.ho
oct.o: /home/jstolfi/include/bool.ho
oct.o: /home/jstolfi/include/jswsize.ho
oct.o: /home/jstolfi/include/affirm.ho
oct.o: /home/jstolfi/include/fget.ho
oct.o: /home/jstolfi/include/nget.ho
oct.o: /home/jstolfi/include/filefmt.ho
oct.o: /home/jstolfi/include/jsmath.ho
oct.o: oct.ho
oct.o: oct_enum.ho

oct_enum.o: oct_enum.c
oct_enum.o: oct.ho
oct_enum.o: oct_enum.ho
oct_enum.o: /home/jstolfi/include/enum_orbits.ho

oct_shapes.o: oct_shapes.c
oct_shapes.o: oct_shapes.ho
oct_shapes.o: oct.ho
oct_shapes.o: /home/jstolfi/include/bool.ho

quad.o: quad.c
quad.o: /home/jstolfi/include/vec.ho
quad.o: /home/jstolfi/include/bool.ho
quad.o: /home/jstolfi/include/jswsize.ho
quad.o: /home/jstolfi/include/affirm.ho
quad.o: /home/jstolfi/include/jsmath.ho
quad.o: /home/jstolfi/include/fget.ho
quad.o: /home/jstolfi/include/nget.ho
quad.o: /home/jstolfi/include/filefmt.ho
quad.o: quad.ho

