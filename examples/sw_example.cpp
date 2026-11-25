/*
 examples/sw_example.cpp
 Nucleics Align v2.0 C++ local alignment example
 */

#include <cstring>
#include <iostream>
#include "nucleics_align.h"

int main(int argc, char *argv[])
{
  if(argc != 3) {
    std::cout << "Usage: " << argv[0] << " <seq1> <seq2>\n";
    return 1;
  }

  na_scoring_t *sc = na_scoring_new(1, -1, -4, -1);
  na_scoring_set_flags(sc, NA_CASE_SENSITIVE);
  na_aligner_t *al = na_aligner_new();
  na_alignment_t *res = na_alignment_new();

  na_status_t st = na_local_align(al, argv[1], std::strlen(argv[1]),
                                  argv[2], std::strlen(argv[2]), sc);
  if(st != NA_OK) {
    std::cerr << "Error: " << na_strerror(st) << "\n";
    na_alignment_free(res); na_aligner_free(al); na_scoring_free(sc);
    return 1;
  }

  while((st = na_local_next(al, res)) == NA_OK) {
    size_t pos_a = 0, pos_b = 0;
    na_alignment_span(res, &pos_a, NULL, &pos_b, NULL);
    std::cout << "seqA: " << na_alignment_a(res) << " [start:" << pos_a << "]\n";
    std::cout << "seqB: " << na_alignment_b(res) << " [start:" << pos_b << "]\n";
    std::cout << "alignment score: " << na_alignment_score(res) << "\n\n";
  }

  na_alignment_free(res);
  na_aligner_free(al);
  na_scoring_free(sc);
  return st == NA_DONE ? 0 : 1;
}
