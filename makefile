_processor = $(shell uname -p)
EXECUTABLE = main

MOL_VERSION = 0.0.6
MOL_INCLUDE ="\"mol.$(MOL_VERSION).h\""
CPPFLAGS := -D _MOL_VERSION_="\"$(MOL_VERSION)\"" -D _MOL_INCLUDE_=$(MOL_INCLUDE) $(CPPFLAGS)
CPPFLAGS := -D ATOM_PRM="\"atom.$(MOL_VERSION).prm\"" $(CPPFLAGS)


LDFLAGS = -L$(HOME)/lib
CPPFLAGS := -I$(HOME)/include $(CPPFLAGS)

CFLAGS =  -ffast-math  -O3 -Wall -W -Wextra -Wpointer-arith -Wcast-qual -Winline -std=gnu99 -g -O0 -fopenmp -pedantic
LIBS = -lmol.$(MOL_VERSION) -ljansson -lm

OBJS = main.o utils.o rotamers.o graph.o mwis.o rotamer_placement.o pack.o

LIBS := $(LIBS)

$(EXECUTABLE):	$(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS)   $(LIBS) -o $@

%.o: %.c %.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

all:
	$(EXECUTABLE)
$(EXECUTABLE).mpi: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS)  $(LIBS) -o $@
$(EXECUTABLE).bgl: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@
clean:
	$(RM) $(OBJS)
	$(RM) $(EXECUTABLE) $(EXECUTABLE).mpi $(EXECUTABLE).bgl
TAGS:
	ctags $(EXECUTABLE).c ~/src/libmol/trunk/mol.$(MOL_VERSION)/*.c
tags: TAGS
