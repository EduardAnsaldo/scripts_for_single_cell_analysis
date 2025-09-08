#!/usr/bin/python3

import os
import shutil
import sys

# Paths to the source files
template_files = [
    'function_template.r',
    'gProfiler2_functions.r',
    'ClusterProfiler_functions.r'
    
]
notebook_files = [
    'Analysis_template.qmd',
    'TCR_analysis_template.qmd',
]
csv_file = 'annotations.csv'

template_files_folder = '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/analyses_templates/scripts_for_single_cell_analysis/'

# New directory to copy files into
if len(sys.argv) < 2:
    print("Usage: python Initialize_analysis_folder.py <new_directory_path>")
    sys.exit(1)

new_dir = sys.argv[1]
scripts_dir = os.path.join(new_dir, "scripts")
data_dir = os.path.join(new_dir, "data")
results_dir = os.path.join(new_dir, "results")

os.makedirs(scripts_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# Copy template files to scripts_dir
for filename in template_files:
    shutil.copy(os.path.join(template_files_folder, filename), os.path.join(scripts_dir, filename))

# Copy CSV file to scripts_dir
shutil.copy(os.path.join(template_files_folder, csv_file), os.path.join(scripts_dir, csv_file))

# Copy notebook files to scripts_dir
for filename in notebook_files:
    shutil.copy(os.path.join(template_files_folder, filename), os.path.join(scripts_dir, filename))

# Rename the notebook files in scripts_dir
os.rename(
    os.path.join(scripts_dir, 'Analysis_template.qmd'),
    os.path.join(scripts_dir, 'Analysis.qmd')
)
os.rename(
    os.path.join(scripts_dir, 'TCR_analysis_template.qmd'),
    os.path.join(scripts_dir, 'TCR_analysis.qmd')
)

# Initialize renv in the new directory
os.chdir(new_dir)
os.system("R -e 'renv::init()'")

print("Files copied and renamed successfully.")