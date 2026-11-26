# benchmarks/bench_common.sh: what the benchmark harnesses share. Source it
# after setting ROOT (the repository root), NW and SW (the built tools):
#   . "$ROOT/benchmarks/bench_common.sh"
# It also puts this directory on PYTHONPATH, for the Python embedded in the
# harnesses (bench_common.py)

export PYTHONPATH="$ROOT/benchmarks${PYTHONPATH:+:$PYTHONPATH}"

# Stop unless the tools are built and GNU time (CPU time and peak RSS) is
# there
bench_require() {
  if [ ! -x "$NW" ] || [ ! -x "$SW" ]; then
    echo "Build Nucleics Align first with 'make' in $ROOT" >&2
    exit 1
  fi
  if [ ! -x /usr/bin/time ]; then
    echo "GNU time (/usr/bin/time) is required for CPU time and peak RSS" >&2
    exit 1
  fi
}

# The machine, compilers and builds the timings come from
bench_header() {
  echo "# host: $(uname -a)"
  echo "# cpu:  $(lscpu | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
  echo "# governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
  echo "# smt:  $(cat /sys/devices/system/cpu/smt/active 2>/dev/null || echo unknown)"
  echo "# gcc:  $(gcc --version | head -n1)"
  echo "# python: $(python3 --version)"
  echo "# nucleics-align: $(git -C "$ROOT" rev-parse --short HEAD)"
  echo "# nw_sha256: $(sha256sum "$NW" | awk '{print $1}')"
  echo "# sw_sha256: $(sha256sum "$SW" | awk '{print $1}')"
}

# fetch_ncbi <db> <accession> <file>: a FASTA record from NCBI (db nuccore
# or protein), unless the file is already there. Paced to NCBI's limit of
# three requests a second without a key
fetch_ncbi() {
  local db="$1" acc="$2" dst="$3"
  if [ ! -s "$dst" ]; then
    curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${db}&id=${acc}&rettype=fasta&retmode=text" -o "$dst"
    sleep 0.4
  fi
}
