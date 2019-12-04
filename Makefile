CXX=clang++
RM=rm -f
CPPFLAGS=-g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs)
INC=-I$${HOME}/.local/igraph-0.7.1/include/igraph/
LIB=-L$${HOME}/.local/igraph-0.7.1/lib/

SRCS=./src/toycities.cpp
OBJS=./run

all:
	clang++ src/toycities.cpp $(INC) $(LIB) -ligraph -o run

run:
	clang++ src/toycities.cpp $(INC) $(LIB) -ligraph -o run; ./run

clean:
	$(RM) run
