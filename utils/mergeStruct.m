function original = mergeStruct(original, new)
	% Will add fields of new to original if they dont exist in original
	n = fieldnames(new);
	for i = 1:length(n)
		if ~isfield(original, n{i})
			original.(n{i}) = new.(n{i});
		end
	end
end