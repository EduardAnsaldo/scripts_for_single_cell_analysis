#!/usr/bin/python3

import os
import shutil
import sys

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

