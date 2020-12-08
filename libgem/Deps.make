gem.ho: gem.h
gem.ho: /home/jstolfi/include/vec.ho

gem_bary.ho: gem_bary.h
gem_bary.ho: gem.ho

gem_print_face.ho: gem_print_face.h
gem_print_face.ho: gem.ho

gem_print_graph.ho: gem_print_graph.h
gem_print_graph.ho: gem.ho

gem.o: gem.c
gem.o: /home/jstolfi/include/bool.ho
gem.o: /home/jstolfi/include/affirm.ho
gem.o: gem.ho

gem_bary.o: gem_bary.c
gem_bary.o: /home/jstolfi/include/affirm.ho
gem_bary.o: gem.ho
gem_bary.o: gem_bary.ho

gem_print_face.o: gem_print_face.c
gem_print_face.o: /home/jstolfi/include/jsfile.ho
gem_print_face.o: gem.ho
gem_print_face.o: gem_print_face.ho

gem_print_graph.o: gem_print_graph.c
gem_print_graph.o: /home/jstolfi/include/jsfile.ho
gem_print_graph.o: gem.ho
gem_print_graph.o: gem_print_graph.ho

