LIBS_PATH=libs

ifdef DEBUG
	OPT = -O0 -g -ggdb
else
	OPT = -O3
endif

# Build with sanitizers: make SANITIZE=1 test
ifdef SANITIZE
	OPT = -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined \
	      -fno-sanitize-recover=all
	TGTFLAGS += -fsanitize=address,undefined -fno-sanitize-recover=all
endif

# Tune for this machine's CPU (e.g. 8-lane AVX2 for the SIMD fill):
# make NATIVE=1. The binaries may not run on other CPUs
ifdef NATIVE
	OPT += -march=native
endif

# Force the 32-bit lanes per SIMD vector (4, 8 or 16), to test the widths
# this CPU would not pick: make LANES=16 test. Not with clang on x86-64
# Linux, whose dispatch builds take only 4
ifdef LANES
	OPT += -DALN_VL=$(LANES)
endif

CFLAGS ?= -Wall -Wextra
# -Wno-psabi: GCC notes that 8-lane vectors are passed differently without
# AVX; the SIMD helpers are always inlined, so no call passes one
CFLAGS += -std=c11 $(OPT) -Wno-psabi
OBJFLAGS = -fPIC
LINKFLAGS = -lalign -lstrbuf -lpthread -lz

INCS=-I $(LIBS_PATH) -I src
LIBS=-L $(LIBS_PATH)/string_buffer -L src
LINK=-lalign -lstrbuf -lpthread -lz

# Compile and bundle all non-main files into library
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)
HDRS=$(wildcard src/*.h)

all: bin/needleman_wunsch bin/smith_waterman bin/lcs src/libalign.a examples

# Rebuild everything when the compiler or flags change (e.g. SANITIZE=1
# after a plain build), rather than mixing objects built different ways
FLAGS_STAMP = src/.build-flags
BUILD_FLAGS = $(CC) $(CFLAGS) $(OBJFLAGS) $(TGTFLAGS)
# Updated as the Makefile is read (not by a dry run, make -n), and written
# again by the rule below if a target such as clean removes it
ifeq (,$(findstring n,$(firstword -$(MAKEFLAGS))))
ifneq ($(BUILD_FLAGS),$(shell cat $(FLAGS_STAMP) 2>/dev/null))
$(shell echo '$(BUILD_FLAGS)' > $(FLAGS_STAMP))
endif
endif
$(FLAGS_STAMP):
	echo '$(BUILD_FLAGS)' > $@

# Build libraries only if they're downloaded
src/libalign.a: $(OBJS)
	[ -d libs/string_buffer ] && cd libs && $(MAKE)
	rm -f $@
	ar -csr $@ $(OBJS)

# -MMD -MP: each object also depends on the headers it includes (src/*.d)
%.o: %.c $(FLAGS_STAMP)
	$(CC) $(CFLAGS) $(OBJFLAGS) $(INCS) -MMD -MP -c $< -o $@

-include $(OBJS:.o=.d)

# The tools and tests link the same library the tests check
bin/needleman_wunsch: src/tools/nw_cmdline.c src/libalign.a $(HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) $(LINKFLAGS)

bin/smith_waterman: src/tools/sw_cmdline.c src/libalign.a $(HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) $(LINKFLAGS)

bin/lcs: src/tools/lcs_cmdline.c src/libalign.a $(HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) $(LINKFLAGS)

bin/nucleics_align_tests: src/tools/tests.c src/libalign.a $(HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) $(LINK)

bin/test_v130: tests/test_v130.c src/libalign.a $(HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(INCS) $(LIBS) $(LINK)

examples: src/libalign.a
	cd examples; $(MAKE) LIBS_PATH=$(abspath $(LIBS_PATH)) SANFLAGS="$(TGTFLAGS)"

bin:
	mkdir -p bin

clean:
	rm -rf bin src/*.o src/*.d src/libalign.a src/.build-flags
	cd examples && $(MAKE) clean

test: bin/nucleics_align_tests bin/test_v130 bin/needleman_wunsch bin/smith_waterman bin/lcs
	./bin/nucleics_align_tests
	./bin/test_v130
	sh tests/cli/run.sh

.PHONY: all clean examples test
