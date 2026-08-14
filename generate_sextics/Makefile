###############################################
# C
###############################################


CC :=  gcc -fopenmp 
BASE_CFLAGS := -Wall -DGCC_INLINE -DEXCLUDE_MAIN 
CFLAGS := -O3 -DPRINT $(BASE_CFLAGS) -DDEBUG=0 #-g -fsanitize=address

###############################################
# Pari configs
###############################################
#PARIHOME=$(GP_INSTALL_PREFIX)
PARIHOME=/usr
PARIDIR=$(PARIHOME)/include/pari
PARILIB=$(PARIHOME)/lib

LDFLAGS := -L$(PARILIB) -Wl,-rpath,$(PARILIB)
LDLIBS := -lgmp -lpari -lm

###############################################
# Detect files and set paths accordingly
###############################################
 
LIBDIR := lib
SRCDIR := src
BUILDDIR := build
TARGET := sextic

SRCFILES := $(wildcard ${SRCDIR}/*.c)
OBJFILES := $(patsubst ${SRCDIR}/%.c, ${BUILDDIR}/%.o, ${SRCFILES})

# excluding main 
CORE_OBJFILES := $(filter-out ${BUILDDIR}/sextic_generator.o, ${OBJFILES})

INCLUDES := -I$(LIBDIR) -I$(PARIDIR)

MKDIR_P = mkdir -p

###############################################
# Compile files
###############################################
all: ${BUILDDIR} ${OBJFILES} ${TARGET}


${BUILDDIR}/%.o : ${SRCDIR}/%.c
	${CC} ${CFLAGS} $(INCLUDES) -c -o $@ $<


${TARGET} : ${OBJFILES}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $^ ${LDLIBS}

.PHONY: all clean debug profile noprint


debug: CFLAGS := -g -DPRINT -DDEBUG=1 $(BASE_CFLAGS)
debug: all

profile: CFLAGS := -O3 -g -pg -DPRINT $(BASE_CFLAGS) -DDEBUG=0
profile: LDFLAGS += -pg
profile: all

noprint: CFLAGS := -O3 $(BASE_CFLAGS) -DDEBUG=0
noprint: all

##################################################
# TEST BUILDER
###################################################
TEST_DIR := test
TEST_SRC_DIR := ${TEST_DIR}/src
TEST_LIB_DIR := ${TEST_DIR}/lib
TEST_BUILD_DIR := ${BUILDDIR}/tests

TEST_LIBS := -I${TEST_LIB_DIR} ${INCLUDES}

TEST_SRCFILES := $(wildcard ${TEST_SRC_DIR}/*.c)
TEST_OBJFILES := $(patsubst ${TEST_SRC_DIR}/%.c, ${TEST_BUILD_DIR}/%.o, ${TEST_SRCFILES})

TEST_CASES := $(patsubst ${TEST_BUILD_DIR}/%.o, ${TEST_DIR}/%, ${TEST_OBJFILES})

#################
# Tests Compiler
#################
.PHONY: tests
# Set noprint
tests: CFLAGS := -O3 $(BASE_CFLAGS) -DDEBUG=0
tests: ${BUILDDIR} ${TEST_BUILD_DIR} ${TEST_CASES}

${TEST_DIR}/% : ${TEST_BUILD_DIR}/%.o ${CORE_OBJFILES} 
	@echo "[Building Test Case]" $@
	${CC} ${CFLAGS} ${LIBS} ${TEST_LIBS} -o $@ $^ ${LDLIBS}

${TEST_BUILD_DIR}/%.o: ${TEST_SRC_DIR}/%.c
	@echo "[Building Test]" $@
	${CC} ${CFLAGS} ${LIBS} ${TEST_LIBS} -c -o $@ $^

#################
# Creates the build directories we need
#################

${BUILDDIR}:
	${MKDIR_P} ${BUILDDIR}

${TEST_BUILD_DIR}:
	${MKDIR_P} ${TEST_BUILD_DIR}
#################

clean:
	rm -f ${BUILDDIR}/*.o *.out
	rm -f ${TEST_BUILD_DIR}/*.o *.out
	rm -f ${TEST_CASES}
	rm -f ${TARGET}
	rmdir  ${TEST_BUILD_DIR}
	rmdir  ${BUILDDIR}
