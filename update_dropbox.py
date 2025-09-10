#!/usr/bin/python3

import os
import sys
import subprocess

# Map parent folder names to Dropbox paths
path_map = {
    'DeCarvalho': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/De Carvalho lab/Lab-shared/Experiments/Eduard A. Gine (EG)',
    'Constantinides': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/Constantinides lab/Experiments/Eduard',
    'Ramanan': '/Users/eduardansaldo/Scripps Research Dropbox/Eduard Ansaldo Gine/Ramanan'
}

# List of patterns to exclude
exclude_patterns = [
    'scripts/*cache/',
    'scripts/*files/',
    '*.DS_Store',
    '*.html.md'
]

def main():
    current_dir = os.path.abspath(os.getcwd())
    parent_dir = os.path.basename(os.path.dirname(current_dir))
    folder_name = os.path.basename(current_dir)

    if parent_dir not in path_map:
        print(f"Parent folder '{parent_dir}' not in path_map. Exiting.")
        sys.exit(1)

    dest_parent = path_map[parent_dir]
    dest_dir = os.path.join(dest_parent, folder_name)

    # Build rsync command
    cmd = [
        'rsync', '-avL', '--delete'
    ]
    for pattern in exclude_patterns:
        cmd.extend(['--exclude', pattern])
    cmd.extend([current_dir + '/', dest_dir + '/'])

    print("Running command:", ' '.join(cmd))
    subprocess.run(cmd)

if __name__ == "__main__":
    main()
