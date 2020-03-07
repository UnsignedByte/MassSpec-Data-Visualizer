# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   18:37:12, 28-Jan-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 16:41:11, 05-Mar-2020

import json
import os.path

root = os.path.abspath(os.path.dirname(__file__)) # Root directory


name = input("Result Folder Name: ") # Get file to read

with open(os.path.join(root, 'default.html')) as f:
	default = f.read()


with open(os.path.join(root, 'Results', name, 'output.json')) as f:
	# data = json.load(f)
	default = default.replace('$$datainput$$', f.read())

with open(os.path.join(root, 'jquery-3.4.1.min.js')) as f:
	default = default.replace("$$JQUERY$$", f.read())

with open(os.path.join(root, 'pickr.min.js')) as f:
	default = default.replace("$$PICKRJS$$", f.read())

with open(os.path.join(root, 'monolith.min.css')) as f:
	default = default.replace("$$PICKRCSS$$", f.read())

with open(os.path.join(root, 'font-awesome.min.css')) as f:
	default = default.replace("$$FONTAWESOME$$", f.read())

with open(os.path.join(root, 'bootstrap.min.js')) as f:
	default = default.replace("$$BOOTSTRAPJS$$", f.read())

with open(os.path.join(root, 'bootstrap-theme.min.css')) as f:
	default = default.replace("$$BOOTSTRAPTHEMECSS$$", f.read())

with open(os.path.join(root, 'bootstrap.min.css')) as f:
	default = default.replace("$$BOOTSTRAPCSS$$", f.read())

with open(os.path.join(root, 'clusterize.min.js')) as f:
	default = default.replace("$$CLUSTERIZEJS$$", f.read())

with open(os.path.join(root, 'clusterize.css')) as f:
	default = default.replace("$$CLUSTERIZECSS$$", f.read())

# with open(os.path.join(root, 'jquery-ui.min.js')) as f:
# 	jqueryui = f.read()

# with open(os.path.join(root, 'jquery-ui.min.css')) as f:
# 	jqueryui = f.read()

with open(os.path.join(root, 'Results', name, f'{name}.html'), 'w') as f:
	f.write(default)
	f.close()