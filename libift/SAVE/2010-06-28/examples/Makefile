INC=../include
LIB=../lib
BIN=../bin

#FLAGS= -g -Wall 
FLAGS= -O3 -Wall


$@.c: $(LIB)/libgift.a libps/libps.a $@.c 
	cd ../; make; cd examples
	cd libps; make; cd ../
	gcc $(FLAGS) $@.c -o $(BIN)/$@ -I$(INC) -L$(LIB) -lgift -lm libps/libps.a

clean:
	rm *.pgm; rm *~; rm time.txt



