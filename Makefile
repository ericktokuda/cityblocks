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

CXXFLAGS := -g -Wall -O3

INC := -I$${HOME}/.local/igraph-0.7.1/include/igraph/ \
	-I$${HOME}/projects/Catch2/single_include/ -I./include/
LIB := -L$${HOME}/.local/igraph-0.7.1/lib/ -ligraph -L./lib/ -L./build/

all: ${RUN}

${RUN}: ${OBJ}
	${CXX} $^ -o ${RUN} ${LIB}

${BUILDDIR}/%.o: ${SRCDIR}/%.${SRCEXT}
	$(info Building libs...)
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} -c -o $@ $<

test: ${TESTRUN}

${TESTRUN}: ${OBJ} ${TESTOBJ} ${OBJ}
	${CXX} $^ -o ${TESTRUN} ${LIB}

${BUILDDIR}/%.o: ${TESTDIR}/%.${SRCEXT}
	$(info Building libs...)
	@mkdir -p ${BUILDDIR}
	${CXX} ${CXXFLAGS} ${INC} -c -o $@ $<

clean:
	$(info Cleaning...)
	${RM} -r ${BUILDDIR} ${RUN} ${TESTRUN}

.PHONY: all clean test
