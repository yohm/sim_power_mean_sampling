CPP=$(shell (which icpc 2>/dev/null || which g++-4.9 || which g++) | tail -n 1)
# OPT=-O3 -Wall -DNDEBUG -Wno-c++11-extensions
OPT=-O2 -Wall -std=c++11

all: analyzer.out

HEADERS1=network.hpp
SRCS1=analysis_main.cpp network.cpp
analyzer.out: $(SRCS1) $(HEADERS1) Makefile
	$(CPP) $(OPT) $(INCLUDE) ${SRCS1} -o $@

clean:
	rm -f *.out *~ *.bak
