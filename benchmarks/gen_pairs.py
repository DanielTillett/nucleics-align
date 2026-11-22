#!/usr/bin/env python3
"""Deterministic, platform-independent DNA pair generator for benchmarks.

Generates a random DNA sequence A of length N and a mutated sequence B where:
  - if div == 'rand': B is an independent random sequence of length N
  - otherwise: div is the nominal substitution percentage (e.g. 1, 10, 0.1).
    Each base is substituted to one of the OTHER three bases with probability
    div/100, or deleted/inserted with probability (div/100)/5 each.

Uses Python's MT19937 PRNG seeded with an integer seed, producing identical
sequences across Linux, macOS, and Windows.

Usage: python3 benchmarks/gen_pairs.py <length> <div|rand> <seed> [out.fa]
"""
import random
import sys

BASES = "ACGT"
OTHERS = {b: [x for x in BASES if x != b] for b in BASES}


def gen_pair(length: int, div: str, seed: int):
    rng = random.Random(seed)
    seq_a = "".join(rng.choice(BASES) for _ in range(length))
    if div == "rand":
        seq_b = "".join(rng.choice(BASES) for _ in range(length))
        return seq_a, seq_b
    p_sub = float(div) / 100.0
    p_del = p_sub / 5.0
    p_ins = p_sub / 5.0
    out = []
    for ch in seq_a:
        r = rng.random()
        if r < p_del:
            continue
        elif r < p_del + p_ins:
            out.append(rng.choice(BASES))
            out.append(ch)
        elif r < p_del + p_ins + p_sub:
            out.append(rng.choice(OTHERS[ch]))
        else:
            out.append(ch)
    return seq_a, "".join(out)


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: gen_pairs.py <length> <div|rand> <seed> [out.fa]")
    length = int(sys.argv[1])
    div = sys.argv[2]
    seed = int(sys.argv[3])
    a, b = gen_pair(length, div, seed)
    text = f">seqA_len{length}_div{div}_seed{seed}\n{a}\n>seqB_len{length}_div{div}_seed{seed}\n{b}\n"
    if len(sys.argv) >= 5:
        with open(sys.argv[4], "w") as f:
            f.write(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
