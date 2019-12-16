CXX := clang++
SRCDIR := src
BUILDDIR := build
TARGET := bin/run

RM := rm -f

SRCEXT := cpp
SOURCES := $(shell find ${SRCDIR} -type f -name *.${SRCEXT})
OBJECTS := $(patsubst ${SRCDIR}/%,${BUILDDIR}/%,${SOURCES:.${SRCEXT}=.o})
CXXFLAGS := -g -Wall -O3

INC := -I$${HOME}/.local/igraph-0.7.1/include/igraph/ \
	-I$${HOME}/projects/Catch2/single_include/ -I./include/
LIB := -L$${HOME}/.local/igraph-0.7.1/lib/ -ligraph -L./lib/

all: ${TARGET}

${TARGET}: ${OBJECTS}
	${CXX} $^ -o ${TARGET} ${LIB}

${BUILDDIR}/%.o: ${SRCDIR}/%.${SRCEXT}
	$(info Building libs...)
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} -c -o $@ $<

clean:
	$(info Cleaning...)
	${RM} -r ${BUILDDIR} ${TARGET}

test:
	$(info Running tests...)
	${CXX} ${CXXFLAGS} test/tester.cpp ${INC} ${LIB} -o bin/tester

.PHONY: all clean
