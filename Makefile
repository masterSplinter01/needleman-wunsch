CC=g++
CFLAGS = -Wall -Werror

all: nw

nw: main.o needleman_wunsch.o
	$(CC) $(CFLAGS) main.o needleman_wunsch.o -o nw
main.o: ./src/main.cpp
	$(CC) $(CFLAGS) -c ./src/main.cpp

needleman_wunsch.o: ./src/needleman_wunsch.cpp ./src/needleman_wunsch.h
	$(CC) $(CFLAGS) -c ./src/needleman_wunsch.cpp

clean:
	rm -rf *.o bot

valgrind:
	 valgrind --leak-check=yes --show-reachable=yes ./nw

