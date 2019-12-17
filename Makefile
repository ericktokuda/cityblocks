CXX := clang++

SRCDIR := src
BUILDDIR := build
TESTDIR := test
BINDIR := bin

RUN := ${BINDIR}/run
TESTRUN := ${BINDIR}/test

RM := rm -f
SRCEXT := cpp

SRC := $(shell find ${SRCDIR} -type f -name *.${SRCEXT})
OBJ := $(patsubst ${SRCDIR}/%,${BUILDDIR}/%,${SRC:.${SRCEXT}=.o})

TESTSRC := $(shell find ${TESTDIR} -type f -name *.${SRCEXT})
TESTOBJ := $(patsubst ${TESTDIR}/%,${BUILDDIR}/%,${TESTSRC:.${SRCEXT}=.o})

CXXFLAGS := -g -Wall -O3 -std=c++11

INC := -I$${HOME}/.local/igraph-0.7.1/include/igraph/ -I./include/
INCTEST := -I$${HOME}/projects/Catch2/single_include/
LIB := -L$${HOME}/.local/igraph-0.7.1/lib/ -ligraph -L./lib/ -L./build/

all: ${RUN}

${RUN}: ${OBJ}
	${CXX} $^ -o ${RUN} ${LIB}

${BUILDDIR}/%.o: ${SRCDIR}/%.${SRCEXT}
	$(info Building libs...)
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} -c -o $@ $<

test: ${TESTRUN}

${TESTRUN}: ${TESTOBJ} ${BUILDDIR}/toycities.o
	${CXX} $^ -o ${TESTRUN} ${LIB}

${BUILDDIR}/tests.o: ${TESTDIR}/tests.cpp
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} ${INCTEST} -c -o $@ $<

${BUILDDIR}/tests-main.o: ${TESTDIR}/tests-main.cpp
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} ${INCTEST} -c -o $@ $<

clean:
	$(info Cleaning...)
	${RM} -r ${BUILDDIR} ${RUN} ${TESTRUN}

.PHONY: all clean test
