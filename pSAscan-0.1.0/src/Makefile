LIB_DIR = ${HOME}/lib                                                                                                                                   
INC_DIR = ${HOME}/include
SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++0x -DNDEBUG -O3 -march=native -I$(INC_DIR) -L$(LIB_DIR)

all: psascan
psascan: main.cpp
	$(CC) $(CFLAGS) -o psascan ./psascan_src/utils.cpp main.cpp -ldivsufsort -ldivsufsort64 -fopenmp 
clean:
	/bin/rm -f *.o
nuclear:
	/bin/rm -f psascan *.o
