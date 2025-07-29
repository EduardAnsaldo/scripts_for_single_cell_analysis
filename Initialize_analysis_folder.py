#!/usr/bin/python3

import os
import shutil
import sys

# Paths to the source files
template_files = [
    'Analysis_template.ipynb',
    'function_template.r',
    'TCR_analysis_template.ipynb',
]
csv_file ='annotations.csv'

template_files_folder = '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/analyses_templates/scripts_for_single_cell_analysis/'

# New directory to copy files into
if len(sys.argv) < 2:
    print("Usage: python Initialize_analysis_folder.py <new_directory_path>")
    sys.exit(1)

new_dir = sys.argv[1]
scripts_dir = os.path.join(new_dir, "scripts")
os.makedirs(scripts_dir, exist_ok=True)

# Copy template files
for filename in template_files:
    shutil.copy(os.path.join(template_files_folder, filename), os.path.join(scripts_dir, filename))

# Copy CSV file
shutil.copy(os.path.join(template_files_folder, csv_file), os.path.join(scripts_dir, csv_file))

# Rename two of the template files in the new directory
os.rename(
    os.path.join(scripts_dir, 'Analysis_template.ipynb'),
    os.path.join(scripts_dir, 'Analysis.ipynb')
)
os.rename(
    os.path.join(scripts_dir, 'TCR_analysis_template.ipynb'),
    os.path.join(scripts_dir, 'TCR_analysis.ipynb')
)

print("Files copied and renamed successfully.")