# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   18:37:12, 28-Jan-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 14:55:31, 09-Jun-2020

import csv
import json
import base64
import os.path

root = os.path.abspath(os.path.dirname(__file__)) # Root directory


name = input("Result Folder Name: ") # Get file to read

resultsFolder = os.path.join(root, 'Results', name)

with open(os.path.join(root, 'webdefault', 'default.html')) as f:
	default = f.read()

data = {};

# fileIDs

with open(os.path.join(resultsFolder, 'fileIDs.csv')) as f:
	reader = csv.DictReader(f)
	data['ids'] = []
	for row in reader:
		data['ids'].append(row)

# get all raws
rawsdir = os.fsencode(os.path.join(resultsFolder, 'Raws'))

#loop through raws and put in data
for file in os.listdir(rawsdir):
	filename = os.fsdecode(file)
	if filename.endswith(".json"):
		with open(os.path.join(resultsFolder, 'Raws', filename)) as f:
			data = {**data, **json.load(f)};

# Read in js tools

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

with open(os.path.join(root, 'webdefault', 'default.js')) as f:
	default = default.replace('$$DEFAULTJS$$', f.read());
	
with open(os.path.join(root, 'webdefault', 'default.css')) as f:
	default = default.replace('$$DEFAULTCSS$$', f.read());

# with open(os.path.join(root, 'jquery-ui.min.js')) as f:
# 	jqueryui = f.read()

# with open(os.path.join(root, 'jquery-ui.min.css')) as f:
# 	jqueryui = f.read()

default = default.replace('$$datainput$$', json.dumps(data)) # Place data into html file

with open(os.path.join(resultsFolder, f'{name}.html'), 'w') as f:
	f.write(default)
	f.close()