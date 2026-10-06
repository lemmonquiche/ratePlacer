# Makefile for ratePlacer and gfl
#
#   make                  build bin/gfl and bin/ratePlacer
#   make example          build, then run the bundled example (example/)
#   make install          copy the binaries to $(PREFIX)/bin (default /usr/local/bin)
#   make debug            build bin/gfl_debug and bin/ratePlacer_debug with sanitizers
#   make clean            remove build products and example output
#   make help             list these targets
#
# Options (set on the command line, e.g. `make OPENMP=0 CC=gcc-14`):
#   CC=...                C compiler (default: cc)
#   OPENMP=0              build ratePlacer without OpenMP (single-threaded)
#   NATIVE=0              do not pass -march=native (for binaries run on other machines)
#   LIBOMP_PREFIX=...     macOS only: location of libomp if not found through Homebrew
#   PREFIX=...            install location (default: /usr/local)

CC       ?= cc
CFLAGS   ?= -O3 -g
OPENMP   ?= 1
NATIVE   ?= 1
PREFIX   ?= /usr/local
BINDIR   := bin
LDLIBS   := -lm

UNAME_S  := $(shell uname -s)

ifeq ($(NATIVE),1)
  CFLAGS += -march=native
endif

# OpenMP: gcc and Linux clang take -fopenmp. Apple clang needs Homebrew's libomp.
OMP_CFLAGS :=
OMP_LIBS   :=
OMP_NOTE   :=
ifeq ($(OPENMP),1)
  CC_IS_APPLE_CLANG := $(shell $(CC) --version 2>/dev/null | grep -q 'Apple clang' && echo 1)
  ifeq ($(UNAME_S)$(CC_IS_APPLE_CLANG),Darwin1)
    LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(wildcard $(LIBOMP_PREFIX)/include/omp.h),)
      OMP_CFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
      OMP_LIBS   := -L$(LIBOMP_PREFIX)/lib -lomp
    else
      OMP_NOTE := libomp not found, so ratePlacer was built without OpenMP (single-threaded). For multithreading, run 'brew install libomp' and rebuild.
    endif
  else
    OMP_CFLAGS := -fopenmp
    OMP_LIBS   := -fopenmp
  endif
endif

DEBUG_CFLAGS := -O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer

GFL_SRC  := src/get_frac_like.c
GFL_DEPS := $(GFL_SRC) src/tools.h
RP_SRC   := src/RatePlacer.c src/minfunc.c
RP_DEPS  := $(RP_SRC) src/tools.h src/opt.h src/minfunc.h src/jph.h

.PHONY: all gfl ratePlacer example check install uninstall debug clean help

all: $(BINDIR)/gfl $(BINDIR)/ratePlacer

gfl: $(BINDIR)/gfl
ratePlacer: $(BINDIR)/ratePlacer

$(BINDIR):
	mkdir -p $@

$(BINDIR)/gfl: $(GFL_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $(GFL_SRC) $(LDLIBS)

$(BINDIR)/ratePlacer: $(RP_DEPS) | $(BINDIR)
	$(CC) $(CFLAGS) $(OMP_CFLAGS) -o $@ $(RP_SRC) $(OMP_LIBS) $(LDLIBS)
	$(if $(OMP_NOTE),@echo "NOTE: $(OMP_NOTE)")

example: all
	sh example/run_example.sh

check: example

install: all
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(BINDIR)/gfl $(BINDIR)/ratePlacer $(DESTDIR)$(PREFIX)/bin

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/gfl $(DESTDIR)$(PREFIX)/bin/ratePlacer

debug: $(BINDIR)/gfl_debug $(BINDIR)/ratePlacer_debug

$(BINDIR)/gfl_debug: $(GFL_DEPS) | $(BINDIR)
	$(CC) $(DEBUG_CFLAGS) -o $@ $(GFL_SRC) $(LDLIBS)

$(BINDIR)/ratePlacer_debug: $(RP_DEPS) | $(BINDIR)
	$(CC) $(DEBUG_CFLAGS) $(OMP_CFLAGS) -o $@ $(RP_SRC) $(OMP_LIBS) $(LDLIBS)

clean:
	rm -rf $(BINDIR) example/output

help:
	@sed -n '2,15p' Makefile | sed 's/^# \{0,1\}//'
