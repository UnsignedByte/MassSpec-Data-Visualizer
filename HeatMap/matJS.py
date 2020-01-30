# -*- coding: utf-8 -*-
# @Author: UnsignedByte
# @Date:   18:37:12, 28-Jan-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 16:31:46, 29-Jan-2020

import json
import os.path

root = os.path.abspath(os.path.dirname(__file__)) # Root directory


name = input("Result Folder Name: ") # Get file to read

with open(os.path.join(root, 'Results', name, 'output.json')) as f:
	# data = json.load(f)
	data = f.read()

with open(os.path.join(root, 'default.html')) as f:
	default = f.read()

with open(os.path.join(root, 'Results', name, f'{name}.html'), 'w') as f:
	f.write(default.replace('$$datainput$$', data))
	f.close()