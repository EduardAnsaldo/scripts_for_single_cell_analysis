#!/usr/bin/python3

import os
import shutil
import sys
import argparse

# Paths to the source files
template_files = [
    'function_template.r',
    'gProfiler2_functions.r',
    'ClusterProfiler_functions.r',
    'install_packages.r',
    'package_list.r'
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
# os.makedirs(data_dir, exist_ok=True)
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
# os.system("Rscript scripts/install_packages.r")
print("Files copied and renamed successfully.")

# Initialize git in the new directory
os.chdir(new_dir)
os.system("git init")
print("Git repository initialized successfully.")

# Edit .gitignore 
os.chdir(new_dir)
with open('.gitignore', 'a') as gitignore:
    gitignore.write('\nresults/**\ndata/**\n*.html\n*.png\nscripts/*cache/**\nscripts/*files/**\n')
print(".gitignore updated successfully.")

# Perform first commit
os.system('''R -e 'usethis::use_readme_md(open=FALSE)' ''')
os.system("git add .")
os.system("git commit -m 'Initial commit: Added template scripts and notebooks'")
print("First commit made successfully.")

# Copy everything and create symbolic link to the data folder of the corresponding directory in dropbox
def copy_with_symlink(src, dst, exclude_dir='data'):
    if not os.path.isdir(src):
        print(f"Source directory '{src}' does not exist.")
        sys.exit(1)
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s_item = os.path.join(src, item)
        d_item = os.path.join(dst, item)
        if item == exclude_dir and os.path.isdir(s_item):
            # Create a symlink instead of copying
            os.symlink(os.path.abspath(s_item), d_item)
        elif os.path.isdir(s_item):
            shutil.copytree(s_item, d_item, symlinks=True)
        else:
            shutil.copy2(s_item, d_item)

path_map = {
    'DeCarvalho': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/De Carvalho lab/Lab-shared/Experiments/Eduard A. Gine (EG)',
    'Constantinides': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/Constantinides lab/Experiments/Eduard',
    'Ramanan': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/Ramanan'
}

parent_dir = os.path.dirname(os.getcwd())
parent_dir_name = os.path.basename(parent_dir)

# Check if parent_dir is one of the three paths
for key, path in path_map.items():
    if os.path.abspath(parent_dir) == os.path.abspath(path):
        # Create a new directory under /Users/eduardansaldo/Documents/<key>/<current_dir>
        current_dir_name = os.path.basename(os.getcwd())
        base_dst = os.path.join('/Users/eduardansaldo/Documents', key, current_dir_name)
        if not os.path.exists(base_dst):
            os.makedirs(base_dst)
        src = os.getcwd()
        dst = base_dst
        print(f"Source path: {src}")
        print(f"Destination path: {dst}")
        copy_with_symlink(src, dst, exclude_dir='data')
        break
else:
    print(f"Parent directory '{parent_dir}' does not match any known experiment paths.")


# Push changes to a new GitHub repository 
os.chdir(base_dst)
print(f"Current working directory for GitHub push: {os.getcwd()}")
os.system("git commit -a -m 'Initial commit: Added template scripts and notebooks'")
os.system("R -e 'usethis::use_github(private=TRUE)'") 