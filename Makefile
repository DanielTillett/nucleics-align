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
# AVX; the SIMD helpers are always inlined, so no call passes one.
# override: appended even to a CFLAGS given on the make command line
override CFLAGS += -std=c11 $(OPT) -Wno-psabi
OBJFLAGS = -fPIC

# The C++ examples use the compiler matching CC unless CXX is given, so a
# sanitizer build links one runtime (clang's objects with clang's)
ifeq ($(origin CXX),default)
  ifneq (,$(findstring clang,$(CC)))
    CXX = clang++
  endif
endif

LIB_INCS = -I include -I src/lib
CLI_INCS = -I include -I src/lib -I src/cli
LIBS = -L src
LIB_LINK = -lalign -lpthread
CLI_LINK = -lalign -lpthread -lz

# Core library sources (zero CLI, zero zlib, zero libs/ dependencies)
LIB_SRCS = $(wildcard src/lib/*.c)
LIB_OBJS = $(LIB_SRCS:.c=.o)
LIB_HDRS = $(wildcard include/*.h) $(wildcard src/lib/*.h)

# Shared CLI support objects (linked only into CLI tools & CLI-scoring tests)
CLI_COMMON_SRCS = src/cli/alignment_cmdline.c src/cli/alignment_scoring_load.c src/cli/cli_seqio.c
CLI_COMMON_OBJS = $(CLI_COMMON_SRCS:.c=.o)
CLI_HDRS = $(wildcard src/cli/*.h)

ALL_OBJS = $(LIB_OBJS) $(CLI_COMMON_OBJS)

all: bin/needleman_wunsch bin/smith_waterman bin/lcs src/libalign.a examples

# Rebuild everything when the compiler or flags change (e.g. SANITIZE=1
# after a plain build), rather than mixing objects built different ways
FLAGS_STAMP = src/.build-flags
BUILD_FLAGS = $(CC) $(CFLAGS) $(OBJFLAGS) $(TGTFLAGS)
ifeq (,$(findstring n,$(firstword -$(MAKEFLAGS))))
ifneq ($(BUILD_FLAGS),$(shell cat $(FLAGS_STAMP) 2>/dev/null))
$(shell echo '$(BUILD_FLAGS)' > $(FLAGS_STAMP))
endif
endif
$(FLAGS_STAMP):
	echo '$(BUILD_FLAGS)' > $@

src/libalign.a: $(LIB_OBJS)
	rm -f $@
	ar -csr $@ $(LIB_OBJS)

# -MMD -MP: each object also depends on the headers it includes (*.d)
src/lib/%.o: src/lib/%.c $(FLAGS_STAMP)
	$(CC) $(CFLAGS) $(OBJFLAGS) $(LIB_INCS) -MMD -MP -c $< -o $@

src/cli/%.o: src/cli/%.c $(FLAGS_STAMP)
	$(CC) $(CFLAGS) $(OBJFLAGS) $(CLI_INCS) -MMD -MP -c $< -o $@

-include $(ALL_OBJS:.o=.d)

bin/needleman_wunsch: src/cli/nw_cmdline.c $(CLI_COMMON_OBJS) src/libalign.a $(LIB_HDRS) $(CLI_HDRS) | bin
	$(CC) -o $@ $< $(CLI_COMMON_OBJS) $(CFLAGS) $(TGTFLAGS) $(CLI_INCS) $(LIBS) $(CLI_LINK)

bin/smith_waterman: src/cli/sw_cmdline.c $(CLI_COMMON_OBJS) src/libalign.a $(LIB_HDRS) $(CLI_HDRS) | bin
	$(CC) -o $@ $< $(CLI_COMMON_OBJS) $(CFLAGS) $(TGTFLAGS) $(CLI_INCS) $(LIBS) $(CLI_LINK)

bin/lcs: src/cli/lcs_cmdline.c $(CLI_COMMON_OBJS) src/libalign.a $(LIB_HDRS) $(CLI_HDRS) | bin
	$(CC) -o $@ $< $(CLI_COMMON_OBJS) $(CFLAGS) $(TGTFLAGS) $(CLI_INCS) $(LIBS) $(CLI_LINK)

bin/nucleics_align_tests: tests/tests.c $(CLI_COMMON_OBJS) src/libalign.a $(LIB_HDRS) $(CLI_HDRS) | bin
	$(CC) -o $@ $< $(CLI_COMMON_OBJS) $(CFLAGS) $(TGTFLAGS) $(CLI_INCS) $(LIBS) $(CLI_LINK)

bin/test_v130: tests/test_v130.c src/libalign.a $(LIB_HDRS) | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) $(LIB_INCS) $(LIBS) $(LIB_LINK)

bin/test_v2_api: tests/test_v2_api.c src/libalign.a include/nucleics_align.h | bin
	$(CC) -o $@ $< $(CFLAGS) $(TGTFLAGS) -I include $(LIBS) $(LIB_LINK)

examples: src/libalign.a
	$(MAKE) -C examples SANFLAGS="$(TGTFLAGS)" CC="$(CC)" CXX="$(CXX)"

bin:
	mkdir -p bin

clean:
	rm -rf bin src/*.o src/*.d src/lib/*.o src/lib/*.d src/cli/*.o src/cli/*.d src/libalign.a src/.build-flags
	cd examples && $(MAKE) clean

test: bin/nucleics_align_tests bin/test_v130 bin/test_v2_api bin/needleman_wunsch bin/smith_waterman bin/lcs
	./bin/nucleics_align_tests
	./bin/test_v130
	./bin/test_v2_api
	sh tests/cli/run.sh

.PHONY: all clean examples test
