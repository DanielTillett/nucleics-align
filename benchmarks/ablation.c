/*
 * benchmarks/ablation.c
 * Ten-configuration component ablation benchmark for Nucleics Align (C11 API).
 * Accepts a single configuration name so each configuration can be timed in a
 * fresh process with exact peak RSS (MiB) and zero cross-run allocator effects.
 *
 * Build (from repo root after `make`):
 *   cc -O3 benchmarks/ablation.c -I src -I libs -L src -L libs/string_buffer \
 *      -lalign -lstrbuf -lpthread -lz -o /tmp/ablation
 * Usage:
 *   /tmp/ablation <pair.fa> <mode|all> [threads]
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "needleman_wunsch.h"
#include "hirschberg.h"

static double now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

static char *read_fasta_seq(FILE *f, size_t *out_len) {
  size_t cap = 4096, len = 0;
  char *buf = malloc(cap);
  if (!buf) exit(1);
  char line[8192];
  int started = 0;
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '>') {
      if (started) {
        fseek(f, -(long)strlen(line), SEEK_CUR);
        break;
      }
      started = 1;
      continue;
    }
    started = 1;
    for (char *p = line; *p; p++) {
      if (*p == '\n' || *p == '\r' || *p == ' ' || *p == '\t') continue;
      if (len + 2 > cap) {
        cap *= 2;
        buf = realloc(buf, cap);
        if (!buf) exit(1);
      }
      buf[len++] = *p;
    }
  }
  buf[len] = '\0';
  *out_len = len;
  return buf;
}

static double time_nw(const char *a, size_t la, const char *b, size_t lb,
                      const scoring_t *sc, bool compact,
                      const alignment_options_t *opts, score_t *out_score) {
  nw_aligner_t *nw = needleman_wunsch_new();
  nw->compact = compact;
  nw->options = opts;
  alignment_t *aln = alignment_create(la + lb + 16);
  double t0 = now_ms();
  aln_error_t err;
  if (needleman_wunsch_align_checked(a, b, la, lb, sc, nw, aln, &err) != ALN_OK) {
    alignment_fprint_error(stderr, &err);
    exit(2);
  }
  double dt = now_ms() - t0;
  *out_score = aln->score;
  alignment_free(aln);
  needleman_wunsch_free(nw);
  return dt;
}

static double time_hb(const char *a, size_t la, const char *b, size_t lb,
                      const scoring_t *sc, const alignment_options_t *opts,
                      score_t *out_score) {
  alignment_t *aln = alignment_create(la + lb + 16);
  double t0 = now_ms();
  aln_error_t err;
  if (hirschberg_align_checked(a, b, la, lb, sc, opts, aln, &err) != ALN_OK) {
    alignment_fprint_error(stderr, &err);
    exit(2);
  }
  double dt = now_ms() - t0;
  *out_score = aln->score;
  alignment_free(aln);
  return dt;
}

static void run_one(const char *mode, const char *a, size_t la,
                    const char *b, size_t lb, const scoring_t *sc, int threads) {
  alignment_options_t o;
  alignment_options_default(&o);
  score_t s = 0;
  double ms = 0.0;

  if (strcmp(mode, "scalar_3mat_1t") == 0) {
    o.striped = false; o.band = false; o.diff = false; o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, false, &o, &s);
  } else if (strcmp(mode, "striped_32b_1t") == 0) {
    o.striped = true; o.diff = false; o.band = false; o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else if (strcmp(mode, "striped_diff8_1t") == 0) {
    o.striped = true; o.diff = true; o.band = false; o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else if (strcmp(mode, "banded_nowfa_1t") == 0) {
    o.striped = true; o.diff = true; o.band = true; o.band_wfa_budget = 0;
    o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else if (strcmp(mode, "banded_wfa_1t") == 0) {
    o.striped = true; o.diff = true; o.band = true; o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else if (strcmp(mode, "banded_nt") == 0) {
    o.striped = true; o.diff = true; o.band = true; o.max_threads = threads;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else if (strcmp(mode, "linear_1t") == 0) {
    o.hb_max_threads = 1;
    ms = time_hb(a, la, b, lb, sc, &o, &s);
  } else if (strcmp(mode, "linear_nt") == 0) {
    o.hb_max_threads = threads;
    ms = time_hb(a, la, b, lb, sc, &o, &s);
  } else if (strcmp(mode, "wfa_1t") == 0) {
    o.wfa_align = true; o.max_threads = 1;
    ms = time_nw(a, la, b, lb, sc, true, &o, &s);
  } else {
    fprintf(stderr, "Unknown mode: %s\n", mode);
    exit(1);
  }
  printf("%s\t%.3f\t%d\n", mode, ms, s);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <pair.fa> <mode|all> [threads]\n", argv[0]);
    return 1;
  }
  const char *mode = argv[2];
  int threads = argc >= 4 ? atoi(argv[3]) : 8;
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  size_t la = 0, lb = 0;
  char *a = read_fasta_seq(f, &la);
  char *b = read_fasta_seq(f, &lb);
  fclose(f);

  scoring_t sc;
  scoring_system_default(&sc);

  if (strcmp(mode, "all") == 0) {
    const char *modes[] = {
      "scalar_3mat_1t", "striped_32b_1t", "striped_diff8_1t",
      "banded_nowfa_1t", "banded_wfa_1t", "banded_nt",
      "linear_1t", "linear_nt", "wfa_1t"
    };
    for (size_t i = 0; i < sizeof(modes)/sizeof(modes[0]); i++)
      run_one(modes[i], a, la, b, lb, &sc, threads);
  } else {
    run_one(mode, a, la, b, lb, &sc, threads);
  }

  free(a);
  free(b);
  return 0;
}
