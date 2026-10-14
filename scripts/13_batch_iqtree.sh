#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

IQTREE_BIN="${IQTREE_BIN:-iqtree}"
require_cmd "$IQTREE_BIN"
alignment="$PROJECT_ROOT/results/phylogeny/gubbins/campy.filtered_polymorphic_sites.fasta"
[[ -s "$alignment" ]] || {
  echo "Gubbins filtered polymorphic-sites alignment not found: $alignment (run 12_batch_gubbins.sh first)." >&2
  exit 1
}
out_dir="$PROJECT_ROOT/results/phylogeny/iqtree"
prefix="$out_dir/campy_final"
varsites="$out_dir/campy.filtered_polymorphic_sites.varsites.fasta"
tmp_varsites="${varsites}.tmp.$$"
mkdir -p "$out_dir" "$PROJECT_ROOT/results/final_trees"
trap 'rm -f "$tmp_varsites"' EXIT

awk '
  /^>/ {
    header = substr($0, 2)
    if (header == "") {
      print "Empty FASTA header in alignment" > "/dev/stderr"
      failed = 1
      exit 1
    }
    label = header
    sub(/[[:space:]].*$/, "", label)
    if (tolower(label) == "reference") {
      reference_count++
      current = 0
      seen_header = 1
      next
    }
    names[++n] = header
    current = n
    seen_header = 1
    next
  }
  {
    gsub(/[[:space:]]/, "")
    if (!seen_header && length($0) > 0) {
      print "Sequence data before first FASTA header" > "/dev/stderr"
      failed = 1
      exit 1
    }
    if (current > 0) {
      sequence[current] = sequence[current] toupper($0)
    }
  }
  END {
    if (failed) {
      exit 1
    }
    if (n < 2) {
      print "Alignment must contain at least two non-reference sequences" > "/dev/stderr"
      exit 1
    }
    length_first = length(sequence[1])
    if (length_first == 0) {
      print "Alignment contains no sequence sites" > "/dev/stderr"
      exit 1
    }
    for (i = 2; i <= n; i++) {
      if (length(sequence[i]) != length_first) {
        print "Alignment sequences have different lengths" > "/dev/stderr"
        exit 1
      }
    }
    variable_count = 0
    for (position = 1; position <= length_first; position++) {
      seen_a = seen_c = seen_g = seen_t = 0
      states = 0
      for (i = 1; i <= n; i++) {
        base = substr(sequence[i], position, 1)
        if (base == "A" && !seen_a) { seen_a = 1; states++ }
        else if (base == "C" && !seen_c) { seen_c = 1; states++ }
        else if (base == "G" && !seen_g) { seen_g = 1; states++ }
        else if (base == "T" && !seen_t) { seen_t = 1; states++ }
      }
      if (states > 1) {
        variable[++variable_count] = position
      }
    }
    if (variable_count == 0) {
      print "No variable nucleotide sites found after excluding Reference" > "/dev/stderr"
      exit 1
    }
    for (i = 1; i <= n; i++) {
      print ">" names[i]
      projected = ""
      for (j = 1; j <= variable_count; j++) {
        projected = projected substr(sequence[i], variable[j], 1)
      }
      print projected
    }
    printf "Excluded %d Reference record(s).\n", reference_count > "/dev/stderr"
    printf "Retained %d variable sites from %d alignment columns.\n", \
      variable_count, length_first > "/dev/stderr"
  }
' "$alignment" > "$tmp_varsites"

mv "$tmp_varsites" "$varsites"
"$IQTREE_BIN" \
  -s "$varsites" \
  -st DNA \
  -m MFP+ASC \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  --prefix "$prefix"

for tree in "$prefix".treefile "$prefix".contree; do
  if [[ -s "$tree" ]]; then
    cp "$tree" "$PROJECT_ROOT/results/final_trees/"
  fi
done
