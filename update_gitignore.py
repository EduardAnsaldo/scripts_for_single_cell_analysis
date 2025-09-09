#!/usr/bin/python3

import os
import shutil
import sys


# Edit .gitignore 
with open('.gitignore', 'a') as gitignore:
    gitignore.write('\nresults/**\ndata/**\n*.html\n*.png\nscripts/*cache/**\nscripts/*files/**\n')
print(".gitignore updated successfully.")
