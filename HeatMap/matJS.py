# @Author: Edmund Lam <edl>
# @Date:   17:22:19, 20-Nov-2019
# @Filename: jsonJS.py
# @Last modified by:   edl
# @Last modified time: 22:16:05, 20-Nov-2019

#pipe mat objects to js file
import tables

file = tables.open_file('output.mat');
print(file)
