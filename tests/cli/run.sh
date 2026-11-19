#!/bin/sh
# Compare the command line tools' stdout, stderr and exit status with the
# expected files in tests/cli/expected (see tests/cli/cases).
#   tests/cli/run.sh            check (run from the repository root)
#   tests/cli/run.sh --update   write the expected files from this build
# The build time in the usage text is replaced by <date>.

cd "$(dirname "$0")/../.." || exit 1
dir=tests/cli
update=0
[ "$1" = "--update" ] && update=1
[ $update = 1 ] && mkdir -p $dir/expected
tmp=$(mktemp -d) || exit 1
trap 'rm -rf "$tmp"' EXIT

pass=0
fail=0

# run <command> <file prefix>: output to <prefix>.out, .err and .status
run() {
  sh -c "$1" > "$2.out" 2> "$2.err.raw" < /dev/null
  echo $? > "$2.status"
  sed 's/(compiled: [^)]*)/(compiled: <date>)/' "$2.err.raw" > "$2.err"
  rm -f "$2.err.raw"
}

same() {
  cmp -s "$1.out" "$2.out" && cmp -s "$1.err" "$2.err" &&
    cmp -s "$1.status" "$2.status"
}

while IFS='|' read -r name threads cmd; do
  name=$(echo $name)
  case "$name" in ''|'#'*) continue ;; esac
  threads=$(echo $threads)
  cmd=$(echo "$cmd" | sed "s#@#$dir/input#g")

  if [ "$threads" = batch ]; then
    runs="1 4"
  else
    runs="-"
  fi
  for t in $runs; do
    c="$cmd"
    label="$name"
    if [ "$t" != - ]; then
      c="$cmd --threads $t"
      label="$name (--threads $t)"
    fi
    run "$c" "$tmp/$name"
    if [ $update = 1 ]; then
      if [ "$t" = - ] || [ "$t" = 1 ]; then
        for e in out err status; do cp "$tmp/$name.$e" "$dir/expected/$name.$e"; done
      elif ! same "$tmp/$name" "$dir/expected/$name"; then
        echo "DIFFERS: $label from --threads 1"; fail=$((fail+1)); continue
      fi
      pass=$((pass+1))
    elif same "$tmp/$name" "$dir/expected/$name"; then
      pass=$((pass+1))
    else
      echo "FAIL: $label"
      for e in out err status; do
        diff "$dir/expected/$name.$e" "$tmp/$name.$e" | head -10
      done
      fail=$((fail+1))
    fi
  done
done < $dir/cases

echo "Command line fixtures: $pass passed, $fail failed"
[ $fail = 0 ]
