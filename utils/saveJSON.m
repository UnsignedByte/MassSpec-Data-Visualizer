function json = saveJSON(fid, obj)
	% \n must become \\n in order for python to parse it properly
	ori = {'\n'};
	new = {'\\n'};
	json = jsonencode(obj);
	for i = 1:length(ori)
		json = strrep(json, ori{i}, new{i});
	end
	fprintf(fid, '%s', json);
end