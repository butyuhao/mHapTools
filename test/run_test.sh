#!/bin/bash
# Regression test of `mhaptools convert`; run from anywhere after building.
#   bash test/run_test.sh [path/to/mhaptools]
# test/data: Left_Ventricle_STL001.IGF2.bam is the wgbs_tools tutorial BAM
# (single-end, bwa-meth, hg19); sim_pe.bam is a simulated paired-end BAM with
# indels, soft clips, low MAPQ and missing mates.
set -uo pipefail
cd "$(dirname "$0")/.."
X=${1:-./mhaptools}
OUT=test/out
mkdir -p "$OUT"
fail=0

run() {   # name, expected file, convert options...
    local name=$1 expected=$2; shift 2
    if "$X" convert "$@" -o "$OUT/$name.mhap" > "$OUT/$name.log" 2>&1 &&
       diff "$expected" "$OUT/$name.mhap" >> "$OUT/$name.log" 2>&1; then
        echo "PASS  $name"
    else
        echo "FAIL  $name (see $OUT/$name.log)"
        fail=1
    fi
}

SE=(-i test/data/Left_Ventricle_STL001.IGF2.bam -c test/data/hg19_CpG.IGF2.gz)
PE=(-i test/data/sim_pe.bam -c test/data/sim_CpG.gz)
run igf2_se      test/expected/igf2_se.mhap      "${SE[@]}"
run sim_pe       test/expected/sim_pe.mhap       "${PE[@]}"
run sim_pe_split test/expected/sim_pe.split.mhap "${PE[@]}" --split
run sim_pe_nondir test/expected/sim_pe.nondir.mhap "${PE[@]}" -n
run sim_pe_region test/expected/sim_pe.region.mhap "${PE[@]}" -r chrS:3001-8000

# bgzipped output with its tabix index
if "$X" convert "${PE[@]}" -o "$OUT/sim_pe.mhap.gz" > "$OUT/sim_pe_gz.log" 2>&1 &&
   [ -s "$OUT/sim_pe.mhap.gz.tbi" ] && gzip -dc "$OUT/sim_pe.mhap.gz" | diff -q test/expected/sim_pe.mhap - > /dev/null; then
    echo "PASS  sim_pe_gz"
else
    echo "FAIL  sim_pe_gz (see $OUT/sim_pe_gz.log)"
    fail=1
fi
exit $fail
