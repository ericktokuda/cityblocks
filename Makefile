CXX := clang++
SRCDIR := src
BUILDDIR := build
TARGET := bin/run

RM := rm -f

SRCEXT := cpp
SOURCES := $(shell find ${SRCDIR} -type f -name *.${SRCEXT})
OBJECTS := $(patsubst ${SRCDIR}/%,${BUILDDIR}/%,${SOURCES:.${SRCEXT}=.o})
CFLAGS := -g # -Wall

INC := -I$${HOME}/.local/igraph-0.7.1/include/igraph/
LIB := -L$${HOME}/.local/igraph-0.7.1/lib/ -ligraph

${TARGET}: ${OBJECTS}
	@echo " Linking..."
	@echo " ${CXX} $^ -o ${TARGET} ${LIB}"
	${CXX} $^ -o ${TARGET} ${LIB}

${BUILDDIR}/%.o: ${SRCDIR}/%.${SRCEXT}
	@mkdir -p ${BUILDDIR}
	@echo " ${CXX} ${CFLAGS} ${INC} -c -o $@ $<"; ${CXX} ${CFLAGS} ${INC} -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " ${RM} -r ${BUILDDIR} ${TARGET}"; ${RM} -r ${BUILDDIR} ${TARGET}

test:
	${CXX} ${CFLAGS} test/tester.cpp ${INC} ${LIB} -o bin/tester

.PHONY: all clean
