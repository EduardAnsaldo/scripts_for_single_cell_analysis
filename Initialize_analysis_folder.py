#!/usr/bin/python3

import os
import shutil
import sys
import argparse

# Paths to the source files
template_files = [
    'function_template.r',
    'gProfiler2_functions.r',
    'ClusterProfiler_functions.r'
]
notebook_files = [
    'Analysis_template.qmd',
]
csv_file = 'annotations.csv'

template_files_folder = '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/analyses_templates/scripts_for_single_cell_analysis/'

# Argument parsing
parser = argparse.ArgumentParser(description="Initialize analysis folder.")
parser.add_argument("new_directory_path", help="Path to the new analysis directory")
parser.add_argument("-T", "--tcr", action="store_true", help="Include TCR_analysis_template.qmd")
args = parser.parse_args()

new_dir = args.new_directory_path
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

# Optionally include TCR_analysis_template.qmd
if args.tcr:
    tcr_file = 'TCR_analysis_template.qmd'
    shutil.copy(os.path.join(template_files_folder, tcr_file), os.path.join(scripts_dir, tcr_file))
    # Rename TCR notebook
    os.rename(
        os.path.join(scripts_dir, 'TCR_analysis_template.qmd'),
        os.path.join(scripts_dir, 'TCR_analysis.qmd')
    )

# Rename the main notebook file in scripts_dir
os.rename(
    os.path.join(scripts_dir, 'Analysis_template.qmd'),
    os.path.join(scripts_dir, 'Analysis.qmd')
)

# Initialize renv in the new directory
os.chdir(new_dir)
os.system("R -e 'renv::init()'")
print("Files copied and renamed successfully.")

# Initialize git in the new directory
os.chdir(new_dir)
os.system('''R -e 'renv::install("usethis")' ''')
os.system("R -e 'usethis::use_git()'")
print("Git repository initialized successfully.")

# Edit .gitignore 
os.chdir(new_dir)
with open('.gitignore', 'a') as gitignore:
    gitignore.write('\nresults/**\ndata/**\n*.html\n*.png\n')
print(".gitignore updated successfully.")

# Perform first commit
os.system('''R -e 'usethis::use_readme_md(open=FALSE)' ''')
os.system("git add .")
os.system("git commit -m 'Initial commit: Added template scripts and notebooks'")
print("First commit made successfully.")

# Push changes to a new GitHub repository
os.system("R -e 'usethis::use_github(private=TRUE)'")