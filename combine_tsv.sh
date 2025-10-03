#!/bin/bash
# Combine all TSV files in a directory into one, keeping the header from the first file only

# Directory containing the TSV files
input_dir="./"  # change as needed
outfile="combined_all.tsv"

# Find all .tsv files (sorted)
files=($(find "$input_dir" -maxdepth 1 -type f -name "*.tsv" | sort))

if [[ ${#files[@]} -eq 0 ]]; then
  echo "❌ No TSV files found in $input_dir"
  exit 1
fi

echo "🔄 Combining ${#files[@]} TSV files from $input_dir..."

# Write first file completely (including header)
cat "${files[0]}" > "$outfile"

# Append remaining files skipping header
for f in "${files[@]:1}"; do
  tail -n +2 "$f" >> "$outfile"
done
