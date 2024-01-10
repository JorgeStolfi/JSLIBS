
LIB=./lib
INCLUDE=./include
BIN=./bin
SRC=./src
OBJ=./obj
DEMO=./demo
#FLAGS= -g -Wall 
FLAGS= -O3 -Wall

libgift: $(LIB)/libgift.a
	echo "libgift.a built..."

$(LIB)/libgift.a: \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/queue.o \
$(OBJ)/image.o \
$(OBJ)/adjacency.o \
$(OBJ)/seeds.o \
$(OBJ)/annimg.o \
$(OBJ)/genift.o \

	ar csr $(LIB)/libgift.a \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/queue.o \
$(OBJ)/image.o \
$(OBJ)/adjacency.o \
$(OBJ)/seeds.o \
$(OBJ)/annimg.o \
$(OBJ)/genift.o \

$(OBJ)/common.o: $(SRC)/common.c
	gcc $(FLAGS) -c $(SRC)/common.c -I$(INCLUDE) \
	-o $(OBJ)/common.o 

$(OBJ)/comptime.o: $(SRC)/comptime.c
	gcc $(FLAGS) -c $(SRC)/comptime.c -I$(INCLUDE) \
	-o $(OBJ)/comptime.o 

$(OBJ)/queue.o: $(SRC)/queue.c
	gcc $(FLAGS) -c $(SRC)/queue.c -I$(INCLUDE) \
	-o $(OBJ)/queue.o 

$(OBJ)/image.o: $(SRC)/image.c
	gcc $(FLAGS) -c $(SRC)/image.c -I$(INCLUDE) \
	-o $(OBJ)/image.o 

$(OBJ)/adjacency.o: $(SRC)/adjacency.c
	gcc $(FLAGS) -c $(SRC)/adjacency.c -I$(INCLUDE) \
	-o $(OBJ)/adjacency.o 

$(OBJ)/seeds.o: $(SRC)/seeds.c
	gcc $(FLAGS) -c $(SRC)/seeds.c -I$(INCLUDE) \
	-o $(OBJ)/seeds.o 

$(OBJ)/annimg.o: $(SRC)/annimg.c
	gcc $(FLAGS) -c $(SRC)/annimg.c -I$(INCLUDE) \
	-o $(OBJ)/annimg.o 

$(OBJ)/genift.o: $(SRC)/genift.c
	gcc $(FLAGS) -c $(SRC)/genift.c -I$(INCLUDE) \
	-o $(OBJ)/genift.o 

clean: 
	rm $(LIB)/libgift.a; rm $(OBJ)/*.o; rm $(BIN)/*; rm ./examples/libps/libps.a; rm ./examples/libps/*.o; rm ./examples/libps/*.ho	





