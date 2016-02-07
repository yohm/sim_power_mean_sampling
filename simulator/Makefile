CPP=$(shell (which icpc 2>/dev/null || which g++-4.9 || which g++) | tail -n 1)
# OPT=-O3 -Wall -DNDEBUG -Wno-c++11-extensions
OPT=-O2 -Wall -std=c++11
# OPT=-O0 -Wall -pg -Wno-c++11-extensions
INCLUDE=

all: parametrized_sampling.out parametrized_sampling_truncation.out

HEADERS=network.hpp sampling.hpp
SRCS=network.cpp sampling.cpp
OBJS=network.o sampling.o

network.o: network.cpp network.hpp
	$(CPP) $(OPT) -c $(INCLUDE) network.cpp

sampling.o: sampling.cpp $(HEADERS)
	$(CPP) $(OPT) -c $(INCLUDE) sampling.cpp

parametrized_sampling.out: $(OBJS) $(HEADERS) parametrized_sampling.cpp Makefile
	$(CPP) $(OPT) $(INCLUDE) ${OBJS} parametrized_sampling.cpp -o $@

parametrized_sampling_truncation.out: $(OBJS) $(HEADERS) parametrized_sampling_truncation.cpp Makefile
	$(CPP) $(OPT) $(INCLUDE) ${OBJS} parametrized_sampling_truncation.cpp -o $@

clean:
	rm -f *.out *~ *.bak *.o
