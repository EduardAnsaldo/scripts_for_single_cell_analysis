#!/bin/bash

BASE_DIR="."
DEST_BASE="processing"

mkdir -p "$DEST_BASE"

for folder in "$BASE_DIR"/*/; do
    folder_name=$(basename "$folder")
    echo "Processing folder: $folder_name"
    # Skip 'fastq_files' and 'processing' directories
    if [[ "$folder_name" == "fastq_files" || "$folder_name" == "$DEST_BASE" ]]; then
        continue
    fi

    dest_dir="$DEST_BASE"
    # No need to mkdir -p "$dest_dir" since it's already created above

    if [[ "$folder_name" == aggr* ]]; then
        src1="./$folder_name/outs/count/filtered_feature_bc_matrix"
        src2="./$folder_name/outs/vdj_t/filtered_contig_annotations.csv"
        src3="./$folder_name/outs/web_summary.html"

        [ -d "$src1" ] && cp -r "$src1" "$dest_dir/${folder_name}_filtered_feature_bc_matrix"
        [ -f "$src2" ] && cp "$src2" "$dest_dir/${folder_name}_filtered_contig_annotations.csv"
        [ -f "$src3" ] && cp "$src3" "$dest_dir/${folder_name}_web_summary.html"
    else
        src1="./$folder_name/outs/per_sample_outs/$folder_name/count/sample_filtered_feature_bc_matrix"
        src2="./$folder_name/outs/per_sample_outs/$folder_name/vdj_t/filtered_contig_annotations.csv"
        src3="./$folder_name/outs/per_sample_outs/$folder_name/web_summary.html"

        [ -d "$src1" ] && cp -r "$src1" "$dest_dir/${folder_name}_sample_filtered_feature_bc_matrix"
        [ -f "$src2" ] && cp "$src2" "$dest_dir/${folder_name}_filtered_contig_annotations.csv"
        [ -f "$src3" ] && cp "$src3" "$dest_dir/${folder_name}_web_summary.html"
    fi
done
