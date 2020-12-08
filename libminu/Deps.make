minu_brent.ho: minu_brent.h
minu_brent.ho: minu_gen.ho

minu_gen.ho: minu_gen.h
minu_gen.ho: /home/jstolfi/include/bool.ho

minu_herm.ho: minu_herm.h
minu_herm.ho: minu_gen.ho

minu_js.ho: minu_js.h
minu_js.ho: minu_gen.ho

minu_brent.o: minu_brent.c
minu_brent.o: /home/jstolfi/include/affirm.ho
minu_brent.o: minu_gen.ho
minu_brent.o: minu_brent.ho

minu_gen.o: minu_gen.c
minu_gen.o: minu_gen.ho

minu_herm.o: minu_herm.c
minu_herm.o: minu_gen.ho
minu_herm.o: minu_herm.ho
minu_herm.o: /home/jstolfi/include/affirm.ho

minu_js.o: minu_js.c
minu_js.o: minu_gen.ho
minu_js.o: minu_js.ho
minu_js.o: /home/jstolfi/include/affirm.ho

