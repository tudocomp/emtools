SHELL = /bin/sh
CC = $(CXX)
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -pthread -march=native
AUX_PAR_FLAGS = -fopenmp
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++0x -g2 -O3 -pthread

all: construct_lcp_sequential \
     construct_lcp_parallel \
     construct_plcp_sequential \
     construct_plcp_parallel \
     construct_lcp_from_plcp_sequential \
     construct_lcp_from_plcp_parallel

construct_lcp_sequential:
	$(CC) $(CFLAGS) -o construct_lcp_sequential em_succinct_irreducible_src/utils.cpp main_construct_lcp_array.cpp
construct_lcp_parallel:
	$(CC) $(CFLAGS) -o construct_lcp_parallel em_succinct_irreducible_src/utils.cpp main_construct_lcp_array.cpp $(AUX_PAR_FLAGS)
construct_plcp_sequential:
	$(CC) $(CFLAGS) -o construct_plcp_sequential em_succinct_irreducible_src/utils.cpp main_construct_plcp_bitvector.cpp
construct_plcp_parallel:
	$(CC) $(CFLAGS) -o construct_plcp_parallel em_succinct_irreducible_src/utils.cpp main_construct_plcp_bitvector.cpp $(AUX_PAR_FLAGS)
construct_lcp_from_plcp_sequential:
	$(CC) $(CFLAGS) -o construct_lcp_from_plcp_sequential em_succinct_irreducible_src/utils.cpp main_construct_lcp_from_plcp.cpp
construct_lcp_from_plcp_parallel:
	$(CC) $(CFLAGS) -o construct_lcp_from_plcp_parallel em_succinct_irreducible_src/utils.cpp main_construct_lcp_from_plcp.cpp $(AUX_PAR_FLAGS)

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f construct_lcp_sequential \
        construct_lcp_parallel \
        construct_plcp_sequential \
        construct_plcp_parallel \
        construct_lcp_from_plcp_sequential \
        construct_lcp_from_plcp_parallel \
        *.o
