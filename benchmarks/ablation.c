/*
 * benchmarks/ablation.c
 * Component ablation benchmark for Nucleics Align (C11 library API).
 * Times each optimisation layer on a pair of sequences from a FASTA file
 * (default DNA scoring +1 / -2 / -4 / -1) and asserts identical optimal scores.
 *
 * Build (from repo root after `make`):
 *   cc -O3 benchmarks/ablation.c -I src -I libs -L src -L libs/string_buffer \
 *      -lalign -lstrbuf -lpthread -lz -o /tmp/ablation
 * Usage:
 *   /tmp/ablation <pair.fa> [threads]
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

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <pair.fa> [threads]\n", argv[0]);
    return 1;
  }
  int threads = argc >= 3 ? atoi(argv[2]) : 8;
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror(argv[1]); return 1; }
  size_t la = 0, lb = 0;
  char *a = read_fasta_seq(f, &la);
  char *b = read_fasta_seq(f, &lb);
  fclose(f);

  scoring_t sc;
  scoring_system_default(&sc);

  alignment_options_t o;
  score_t ref_score = 0, s = 0;
  double ms;

  /* 1. Baseline scalar 3-matrix (1 thread) */
  alignment_options_default(&o);
  o.striped = false; o.band = false; o.diff = false; o.max_threads = 1;
  ms = time_nw(a, la, b, lb, &sc, false, &o, &ref_score);
  printf("scalar_3mat_1t\t%.3f\t%d\n", ms, ref_score);

  /* 2. + Compact 1-byte traceback (scalar, 1 thread) */
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("compact_scalar_1t\t%.3f\t%d\n", ms, s);

  /* 3. + Striped SIMD (no diff, no band, 1 thread) */
  o.striped = true; o.diff = false; o.band = false; o.max_threads = 1;
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("striped_32b_1t\t%.3f\t%d\n", ms, s);

  /* 4. + Difference-encoded 8-bit lanes (no band, 1 thread) */
  o.diff = true; o.band = false; o.max_threads = 1;
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("striped_diff8_1t\t%.3f\t%d\n", ms, s);

  /* 5. + Certified exact banding (1 thread) */
  o.band = true; o.max_threads = 1;
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("banded_1t\t%.3f\t%d\n", ms, s);

  /* 6. + Multi-threading (N threads) */
  o.max_threads = threads;
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("banded_%dt\t%.3f\t%d\n", threads, ms, s);

  /* 7. Linear space (1 thread and N threads) */
  alignment_options_default(&o);
  o.hb_max_threads = 1;
  ms = time_hb(a, la, b, lb, &sc, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("linear_1t\t%.3f\t%d\n", ms, s);

  o.hb_max_threads = threads;
  ms = time_hb(a, la, b, lb, &sc, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("linear_%dt\t%.3f\t%d\n", threads, ms, s);

  /* 8. Wavefront (--wfa, 1 thread) */
  alignment_options_default(&o);
  o.wfa_align = true; o.max_threads = 1;
  ms = time_nw(a, la, b, lb, &sc, true, &o, &s);
  if (s != ref_score) { fprintf(stderr, "Score mismatch!\n"); return 3; }
  printf("wfa_1t\t%.3f\t%d\n", ms, s);

  free(a);
  free(b);
  return 0;
}
