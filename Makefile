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

OBJS = main.o stats.o dist.o cov.o trim.o version.o bxhash.o



PACKAGE_VERSION = 0.1

# If building from a Git repository, replace $(PACKAGE_VERSION) with the Git
# description of the working tree: either a release tag with the same value
# as $(PACKAGE_VERSION) above, or an exact description likely based on a tag.
# $(shell), :=, etc are GNU Make-specific.  If you don't have GNU Make,
# comment out this conditional.
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

# If you don't have GNU Make but are building from a Git repository, you may
# wish to replace this with a rule that always rebuilds version.h:
# version.h: force
#   echo '#define BXCHECK_VERSION "`git describe --always --dirty`"' > $@
version.h:
	echo '#define BXCHECK_VERSION "$(PACKAGE_VERSION)"' > $@

print-version:
	@echo $(PACKAGE_VERSION)


.SUFFIXES:.c .o
.PHONY:all clean print-version

force:

.c.o:
	$(CC) $(CFLAGS) $(EXTRA_CPPFLAGS) $(ALL_CPPFLAGS) -c -o $@ $<

version.o: version.h version.c
main.o: main.c version.h
stats.o: stats.c cov.h dist.h bxhash.h $(htslib_sam_h) $(htslib_khash_h)
dist.o: dist.c dist.h
cov.o: cov.c cov.h rbuf.h dist.h
trim.o: trim.c bxhash.h
bxhash.o: bxhash.h bxhash.c $(htslib_khash_h)

bxcheck: $(HTSLIB) $(OBJS)
	$(CC) $(ALL_LDFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread $(HTSLIB_LIBS) $(GSL_LIBS) $(ALL_LIBS)

bxreads: $(HTSLIB) bxreads.c
	$(CC) $(CFLAGS) $(EXTRA_CPPFLAGS) $(ALL_CPPFLAGS) $(ALL_LDFLAGS) -o $@ bxreads.c $(HTSLIB) -lpthread $(HTSLIB_LIBS) $(GSL_LIBS) $(ALL_LIBS)

clean:
	-rm -f gmon.out version.h *.o *~ $(PROG)

