/*
 percore.c
 url: https://github.com/DanielTillett/nucleics-align
 license: BSD 3-Clause, see LICENSE
 date: Nov 2026
 */

// Library speed on one thread and on the default threads, for one pair of
// sequences (the two lines of a file from make_pairs.sh): the full-matrix
// global alignment, the linear-space score, the Smith-Waterman score scan
// and the best local hit. Build and run from the repository root:
//   make src/libalign.a
//   cc -O2 -I include -I src/lib benchmarks/percore.c -L src -lalign
//      -lpthread -o percore && ./percore pair.txt   (one command)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include "hirschberg.h"

static char* read_line(FILE *fh)
{
  size_t cap = 1024, len = 0;
  char *s = malloc(cap);
  int c;
  if(s == NULL) exit(EXIT_FAILURE);
  while((c = fgetc(fh)) != EOF && c != '\n') {
    if(len + 1 == cap && (s = realloc(s, cap *= 2)) == NULL) exit(EXIT_FAILURE);
    s[len++] = (char)c;
  }
  s[len] = '\0';
  return s;
}

static double seconds(clockid_t clock)
{
  struct timespec t;
  clock_gettime(clock, &t);
  return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
  FILE *fh = argc == 2 ? fopen(argv[1], "r") : NULL;
  if(fh == NULL) {
    fprintf(stderr, "usage: %s <pair file>\n", argv[0]);
    return EXIT_FAILURE;
  }
  char *a = read_line(fh), *b = read_line(fh);
  fclose(fh);
  const size_t la = strlen(a), lb = strlen(b);
  const double cells = (double)la * (double)lb;

  scoring_t scoring;
  scoring_init(&scoring, 1, -2, -4, -1, false, false, false, false, false,
               false);
  nw_aligner_t *nw = needleman_wunsch_new();
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *aln = alignment_create(256);
  nw->compact = true;
  smith_waterman_set_compact(sw, true);
  smith_waterman_set_min_score(sw, 1);

  alignment_options_t opts;
  const char *names[] = {"global, full matrix", "global, score only",
                         "local, score scan", "local, best hit"};
  int threads[] = {1, 0}, t, task, rep;

  printf("%zu x %zu cells\n", la, lb);
  for(t = 0; t < 2; t++) {
    alignment_options_default(&opts);
    opts.max_threads = opts.hb_max_threads = threads[t];
    nw->options = &opts;
    smith_waterman_set_options(sw, &opts);
    for(task = 0; task < 4; task++) {
      double best_wall = 1e30, best_cpu = 0;
      score_t score = 0;
      size_t x, y;
      for(rep = 0; rep < 3; rep++) {
        double w0 = seconds(CLOCK_MONOTONIC);
        double c0 = seconds(CLOCK_PROCESS_CPUTIME_ID);
        switch(task) {
          case 0:
            needleman_wunsch_align2(a, b, la, lb, &scoring, nw, aln);
            score = aln->score;
            break;
          case 1:
            hirschberg_score_opts(a, b, la, lb, &scoring, &opts, &score);
            break;
          case 2:
            alignment_sw_scan(a, la, b, lb, &scoring, 0, &opts, &score, &x, &y);
            break;
          case 3:
            smith_waterman_align2(a, b, la, lb, &scoring, sw);
            score = smith_waterman_fetch(sw, aln) ? aln->score : 0;
            break;
        }
        double wall = seconds(CLOCK_MONOTONIC) - w0;
        double cpu = seconds(CLOCK_PROCESS_CPUTIME_ID) - c0;
        if(wall < best_wall) { best_wall = wall; best_cpu = cpu; }
      }
      printf("%-8s %-22s %8.3f s wall %8.3f s CPU %6.2f Gcells/CPU-s  "
             "score %i\n", t == 0 ? "1 thread" : "default", names[task],
             best_wall, best_cpu, cells / best_cpu / 1e9, score);
    }
  }

  needleman_wunsch_free(nw);
  smith_waterman_free(sw);
  alignment_free(aln);
  free(a);
  free(b);
  return EXIT_SUCCESS;
}
