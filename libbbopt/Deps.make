bbgoal.ho: bbgoal.h
bbgoal.ho: /home/jstolfi/include/flt.ho
bbgoal.ho: /home/jstolfi/include/ia.ho

bbopt.ho: bbopt.h
bbopt.ho: fbox.ho
bbopt.ho: fboxlist.ho
bbopt.ho: /home/jstolfi/include/ia.ho
bbopt.ho: /home/jstolfi/include/bool.ho

fbb_f1_ia.ho: fbb_f1_ia.h
fbb_f1_ia.ho: bbgoal.ho

fbb_f2_ia.ho: fbb_f2_ia.h
fbb_f2_ia.ho: bbgoal.ho

fbb_f3_ia.ho: fbb_f3_ia.h
fbb_f3_ia.ho: bbgoal.ho

fbox.ho: fbox.h
fbox.ho: /home/jstolfi/include/ia.ho

fboxheap.ho: fboxheap.h
fboxheap.ho: fbox.ho

fboxlist.ho: fboxlist.h
fboxlist.ho: fbox.ho

bbgoal.o: bbgoal.c
bbgoal.o: bbgoal.ho
bbgoal.o: fbb_f1_ia.ho
bbgoal.o: fbb_f2_ia.ho
bbgoal.o: fbb_f3_ia.ho
bbgoal.o: /home/jstolfi/include/affirm.ho
bbgoal.o: /home/jstolfi/include/bool.ho

bbopt.o: bbopt.c
bbopt.o: fbox.ho
bbopt.o: fboxheap.ho
bbopt.o: fboxlist.ho
bbopt.o: /home/jstolfi/include/aa.ho
bbopt.o: /home/jstolfi/include/ia.ho
bbopt.o: /home/jstolfi/include/affirm.ho
bbopt.o: /home/jstolfi/include/ps.ho
bbopt.o: bbopt.ho

fbb_f1_ia.o: fbb_f1_ia.c
fbb_f1_ia.o: fbb_f1_ia.ho
fbb_f1_ia.o: bbgoal.ho
fbb_f1_ia.o: /home/jstolfi/include/ia.ho

fbb_f2_ia.o: fbb_f2_ia.c
fbb_f2_ia.o: fbb_f2_ia.ho
fbb_f2_ia.o: bbgoal.ho
fbb_f2_ia.o: /home/jstolfi/include/ia.ho
fbb_f2_ia.o: /home/jstolfi/include/bool.ho
fbb_f2_ia.o: /home/jstolfi/include/affirm.ho

fbb_f3_ia.o: fbb_f3_ia.c
fbb_f3_ia.o: fbb_f3_ia.ho
fbb_f3_ia.o: bbgoal.ho
fbb_f3_ia.o: /home/jstolfi/include/ia.ho

fbox.o: fbox.c
fbox.o: /home/jstolfi/include/ia.ho
fbox.o: fbox.ho

fboxheap.o: fboxheap.c
fboxheap.o: fbox.ho
fboxheap.o: fboxheap.ho
fboxheap.o: /home/jstolfi/include/affirm.ho

fboxlist.o: fboxlist.c
fboxlist.o: fbox.ho
fboxlist.o: fboxlist.ho
fboxlist.o: /home/jstolfi/include/affirm.ho

