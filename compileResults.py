# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   18:37:12, 28-Jan-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 14:53:20, 14-Oct-2020

import csv
import json
import base64
import os.path
import traceback

# Colored terminal text for python
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def title(s):
	return s[0].upper() + s[1:]

root = os.path.abspath(os.path.dirname(__file__)) # Root directory

#read in default html page

with open(os.path.join(root, 'webdefault', 'default.html')) as f:
	default = f.read()

with open(os.path.join(root, 'webdefault', 'default.js')) as f:
	default = default.replace('$$DEFAULTJS$$', f.read());

with open(os.path.join(root, 'webdefault', 'default.css')) as f:
	default = default.replace('$$DEFAULTCSS$$', f.read());

# Read in js/css tools

with open(os.path.join(root, 'webtools', 'jquery-3.4.1.min.js')) as f:
	default = default.replace("$$JQUERY$$", f.read())

with open(os.path.join(root, 'webtools', 'pickr.min.js')) as f:
	default = default.replace("$$PICKRJS$$", f.read())

with open(os.path.join(root, 'webtools', 'monolith.min.css')) as f:
	default = default.replace("$$PICKRCSS$$", f.read())

with open(os.path.join(root, 'webtools', 'font-awesome.min.css')) as f:
	default = default.replace("$$FONTAWESOME$$", f.read())

with open(os.path.join(root, 'webtools', 'bootstrap.min.js')) as f:
	default = default.replace("$$BOOTSTRAPJS$$", f.read())

with open(os.path.join(root, 'webtools', 'bootstrap-theme.min.css')) as f:
	default = default.replace("$$BOOTSTRAPTHEMECSS$$", f.read())

with open(os.path.join(root, 'webtools', 'bootstrap.min.css')) as f:
	default = default.replace("$$BOOTSTRAPCSS$$", f.read())

with open(os.path.join(root, 'webtools', 'clusterize.min.js')) as f:
	default = default.replace("$$CLUSTERIZEJS$$", f.read())

with open(os.path.join(root, 'webtools', 'clusterize.css')) as f:
	default = default.replace("$$CLUSTERIZECSS$$", f.read())

# with open(os.path.join(root, 'jquery-ui.min.js')) as f:
# 	jqueryui = f.read()

# with open(os.path.join(root, 'jquery-ui.min.css')) as f:
# 	jqueryui = f.read()

name = input("Result Folder Name: ") # Get file to read

def insertData(name):
	resultsFolder = os.path.join(root, 'Results', name)

	data = {};

	# fileIDs

	with open(os.path.join(resultsFolder, 'fileIDs.csv')) as f:
		reader = csv.DictReader(f)
		data['FileIds'] = []
		for row in reader:
			data['FileIds'].append(row)

	# get all raws
	rawsdir = os.fsencode(os.path.join(resultsFolder, 'Raws'))

	#loop through raws and put in data
	for file in os.listdir(rawsdir):
		filename = os.fsdecode(file)
		if filename.endswith(".json"):
			with open(os.path.join(resultsFolder, 'Raws', filename)) as f:
				print(f"Reading {filename}.");
				loaded = json.load(f);
				if isinstance(loaded, list):
					loaded = {title(filename.rsplit('.',1)[0]):loaded}
				data = {**data, **loaded};

	return default.replace('$$datainput$$', json.dumps(data)) # Place data into html file

def saveData(x):
	print(f"{bcolors.HEADER}Compiling files for {x}.{bcolors.ENDC}")
	with open(os.path.join(os.path.join(root, 'Results', x), f'{x}.html'), 'w') as f:
		try:
			f.write(insertData(x))
			print(f"{bcolors.OKGREEN}Compilation completed successfully.{bcolors.ENDC}")
		except Exception as e:
			print(f"{bcolors.FAIL}Compilation error occured. Results data may not be formatted properly.\nFull error:\n{traceback.format_exc()}{bcolors.ENDC}")
		finally:
			f.close()

if not name: # name is empty
	print(f"{bcolors.WARNING}No dataset supplied. Compiling all datasets by default.{bcolors.ENDC}")
	for x in os.listdir(os.path.join(root, 'Results')):
		if os.path.isdir(os.path.join(root, 'Results',x)):
			saveData(x);
else:
	with open(os.path.join(os.path.join(root, 'Results', name), f'{name}.html'), 'w') as f:
			saveData(name);