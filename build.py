# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   16:05:44, 17-Nov-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 2021-05-02 14:14:44

import tkinter as tk
from tkinter import filedialog
import os.path
from shutil import copy
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

root = tk.Tk()
root.withdraw()
rootDir = os.path.abspath(os.path.dirname(__file__)) # Root directory

files = filedialog.askopenfilenames(parent=root,title='Choose files', initialdir=rootDir)

name = re.match(r"^\d+_(.+?_.+?)_", os.path.basename(files[0])).group(1)
# name = input("Result Folder Name: ") # Get file to read

os.makedirs(os.path.join(rootDir, 'Results', name, 'Data'), exist_ok=True)
os.makedirs(os.path.join(rootDir, 'Results', name, 'Significant'), exist_ok=True)
os.makedirs(os.path.join(rootDir, 'Results', name, 'Raws'), exist_ok=True)


print(f'{bcolors.OKGREEN}Created folder {name}{bcolors.ENDC}')

for x in files:
	copy(x, os.path.join(rootDir, 'Results', name, 'Data', os.path.basename(x)))