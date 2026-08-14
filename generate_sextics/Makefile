###############################################
# C
###############################################

# Added -no-pie to prevent the modern GCC security feature from blocking libpari
CC :=  gcc -fopenmp -no-pie
BASE_CFLAGS := -Wall -DGCC_INLINE -DEXCLUDE_MAIN -fgnu89-inline
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
 
LIBDIR := include
SRCDIR := src
BUILDDIR := build
TARGET := sextic

# WARNING: Make sure there are no backup files like 'sextic-tmp.c' in the src/ folder!
SRCFILES := $(wildcard ${SRCDIR}/*.c)
OBJFILES := $(patsubst ${SRCDIR}/%.c, ${BUILDDIR}/%.o, ${SRCFILES})

INCLUDES := -I$(LIBDIR) -I$(PARIDIR)

MKDIR_P = mkdir -p

###############################################
# Compile files
###############################################
all: ${BUILDDIR} ${OBJFILES} ${TARGET}

# Compiling object files (No libraries needed here, just includes)
${BUILDDIR}/%.o : ${SRCDIR}/%.c
	${CC} ${CFLAGS} $(INCLUDES) -c -o $@ $<

# Linking the final executable (LDFLAGS before objects, LDLIBS strictly at the end!)
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

#################
# Creates the build directories we need
#################

${BUILDDIR}:
	${MKDIR_P} ${BUILDDIR}

#################

clean:
	rm -f ${BUILDDIR}/*.o *.out
	rmdir ${BUILDDIR} 2>/dev/null || true
	rm -f ${TARGET}
