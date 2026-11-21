#!/bin/sh
# Time the command line tools on the benchmark pairs (see README.md here).
#   benchmarks/run.sh [lengths]        e.g. benchmarks/run.sh "5000 20000"
# Run from the repository root after `make`. Needs GNU time (/usr/bin/time
# with -f) for CPU time and peak memory. Prints one tab-separated row per
# run: task, length, divergence, mode, wall s, CPU s, peak MB, score
cd "$(dirname "$0")/.." || exit 1
lengths=${1:-"5000 20000 50000"}
tmp=$(mktemp -d) || exit 1
trap 'rm -rf "$tmp"' EXIT

best_of_3() {
  best=""
  for rep in 1 2 3; do
    if ! /usr/bin/time -f "%e %U %S %M" -o "$tmp/time" "$@" > "$tmp/out" \
                                                           2> "$tmp/err"; then
      echo "Failed: $*" >&2
      cat "$tmp/err" >&2
      return 1
    fi
    read -r wall user sys kb < "$tmp/time"
    if [ -z "$best" ] || awk -v a="$wall" -v b="$best" 'BEGIN { exit !(a < b) }'; then
      best=$wall
      cpu=$(awk -v u="$user" -v s="$sys" 'BEGIN { printf "%.2f", u + s }')
      mb=$((kb / 1024))
    fi
  done
  printf '%s\t%s\t%s\n' "$best" "$cpu" "$mb"
}

printf 'task\tlength\tdivergence\tmode\twall_s\tcpu_s\tpeak_mb\tscore\n'
for n in $lengths; do
  for d in 1 10 rand; do
    sh benchmarks/make_pairs.sh "$n" "$d" 11 > "$tmp/pair.txt"
    for mode in full linear scoreonly local; do
      case $mode in
        full)      task=global; cmd="./bin/needleman_wunsch --printscores" ;;
        linear)    task=global; cmd="./bin/needleman_wunsch --linear --printscores" ;;
        scoreonly) task=score;  cmd="./bin/needleman_wunsch --scoreonly" ;;
        local)     task=local;  cmd="./bin/smith_waterman --maxhits 1" ;;
      esac
      times=$(best_of_3 $cmd --file "$tmp/pair.txt") || exit 1
      score=$(grep -o 'score: -\{0,1\}[0-9]*' "$tmp/out" | head -1 | cut -d' ' -f2)
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$n" "$d" "$mode" "$times" "$score"
    done
  done
done
