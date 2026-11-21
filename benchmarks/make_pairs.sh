#!/bin/sh
# Benchmark inputs: pairs of random DNA sequences, the same on every run.
#   benchmarks/make_pairs.sh <length> <divergence> <seed> > pair.txt
# divergence: 1 (1% substitutions, 0.2% insertions and deletions), 10 (10%,
# 2%) or rand (two unrelated sequences). Writes the two sequences, one per
# line, as the tools read them with --file
length=${1:?length} div=${2:?divergence} seed=${3:-11}
awk -v n="$length" -v div="$div" -v seed="$seed" '
function base() { return substr("acgt", int(rand() * 4) + 1, 1) }
BEGIN {
  srand(seed)
  a = ""
  for(i = 0; i < n; i++) a = a base()
  if(div == "rand") {
    b = ""
    for(i = 0; i < n; i++) b = b base()
  } else {
    sub_rate = div / 100; indel_rate = sub_rate / 5
    b = ""
    for(i = 1; i <= n; i++) {
      r = rand()
      if(r < indel_rate / 2) continue                       # deletion
      if(r < indel_rate) b = b base()                       # insertion
      if(rand() < sub_rate) b = b base(); else b = b substr(a, i, 1)
    }
  }
  print a; print b
}'
