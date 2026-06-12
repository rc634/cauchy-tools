#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BINARY="$SCRIPT_DIR/ahfinder"
BATCH_DIR="${1:-$SCRIPT_DIR/../polly/batch}"

if [[ ! -f "$BINARY" ]]; then
    echo "Error: binary not found at $BINARY — run make first"
    exit 1
fi

if [[ ! -d "$BATCH_DIR" ]]; then
    echo "Error: batch directory not found: $BATCH_DIR"
    exit 1
fi

BATCH_DIR="$(cd "$BATCH_DIR" && pwd)"
subdirs=( "$BATCH_DIR"/*/ )

echo "Batch directory : $BATCH_DIR"
echo "Runs found      : ${#subdirs[@]}"
echo ""

for subdir in "${subdirs[@]}"; do
    name="$(basename "$subdir")"
    echo "=== $name ==="
    "$BINARY" "$subdir" | tee "$subdir/output.txt"
    echo ""
done

echo "All runs complete."
