/*
 examples/nw_example.c
 Nucleics Align v2.0 C global alignment example
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nucleics_align.h"

int main(int argc, char *argv[])
{
  if(argc != 3) {
    printf("usage: ./nw_example <seq1> <seq2>\n");
    return EXIT_FAILURE;
  }

  na_scoring_t *sc = na_scoring_new(1, -2, -4, -1);
  na_aligner_t *al = na_aligner_new();
  na_alignment_t *res = na_alignment_new();
  if(!sc || !al || !res) {
    fprintf(stderr, "Error: out of memory\n");
    na_scoring_free(sc); na_aligner_free(al); na_alignment_free(res);
    return EXIT_FAILURE;
  }

  na_scoring_set_flags(sc, NA_FREE_START_GAP | NA_FREE_END_GAP);
  na_scoring_set_pair(sc, 'a', 'c', -2);
  na_scoring_set_pair(sc, 'c', 'a', -1);

  na_status_t st = na_global_align(al, argv[1], strlen(argv[1]),
                                   argv[2], strlen(argv[2]), sc, res);
  if(st != NA_OK) {
    char msg[128];
    na_error_format(na_aligner_last_error(al), msg, sizeof(msg));
    fprintf(stderr, "%s\n", msg);
    na_scoring_free(sc); na_aligner_free(al); na_alignment_free(res);
    return EXIT_FAILURE;
  }

  printf("seqA: %s\n", na_alignment_a(res));
  printf("seqB: %s\n", na_alignment_b(res));
  printf("alignment score: %d\n", na_alignment_score(res));

  na_alignment_free(res);
  na_aligner_free(al);
  na_scoring_free(sc);
  return EXIT_SUCCESS;
}
