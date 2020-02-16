# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   18:37:12, 28-Jan-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 15:51:49, 15-Feb-2020

import json
import os.path

root = os.path.abspath(os.path.dirname(__file__)) # Root directory


name = input("Result Folder Name: ") # Get file to read

with open(os.path.join(root, 'Results', name, 'output.json')) as f:
	# data = json.load(f)
	data = f.read()

with open(os.path.join(root, 'default.html')) as f:
	default = f.read()

with open(os.path.join(root, 'jquery-3.4.1.min.js')) as f:
	jquery = f.read()

# with open(os.path.join(root, 'jquery-ui.min.js')) as f:
# 	jqueryui = f.read()

# with open(os.path.join(root, 'jquery-ui.min.css')) as f:
# 	jqueryui = f.read()

with open(os.path.join(root, 'Results', name, f'{name}.html'), 'w') as f:
	f.write(default.replace('$$datainput$$', data).replace("$$JQUERY$$", jquery))
	f.close()