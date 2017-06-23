PROG= bxcheck

all:$(PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LDFLAGS = $(HTSLIB_static_LDFLAGS)
HTSLIB_LIBS = $(HTSLIB_static_LIBS)

CC       = gcc
CPPFLAGS =
CFLAGS   = -g -Wall -Wc++-compat -O2
LDFLAGS  =
LIBS     =

# TODO Use configure or htslib.pc to add -rdynamic/-ldl conditionally
ALL_CPPFLAGS = -I. $(HTSLIB_CPPFLAGS) $(CPPFLAGS)
ALL_LDFLAGS  = $(DYNAMIC_FLAGS) $(HTSLIB_LDFLAGS) $(LDFLAGS)
ALL_LIBS     = -lm -lz -ldl $(LIBS)

EXTRA_CPPFLAGS = -I. -I$(HTSDIR)

OBJS = main.o dist.o cov.o

.SUFFIXES:.c .o
.PHONY:all

force:

.c.o:
	$(CC) $(CFLAGS) $(EXTRA_CPPFLAGS) $(ALL_CPPFLAGS) -c -o $@ $<

main.o: main.c $(htslib_sam_h) $(htslib_khash_h)
dist.o: dist.c dist.h
cov.o: cov.c cov.h

bxcheck: $(HTSLIB) $(OBJS)
	$(CC) $(ALL_LDFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread $(HTSLIB_LIBS) $(GSL_LIBS) $(ALL_LIBS)

clean:
	-rm -f gmon.out *.o *~ $(PROG)

