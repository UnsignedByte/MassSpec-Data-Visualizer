function json = saveJSON(fid, obj)
	ori = {};
	new = {};
	json = jsonencode(obj);
	for i = 1:length(ori)
		json = strrep(json, ori{i}, new{i});
	end
	fprintf(fid, '%s', json);
end