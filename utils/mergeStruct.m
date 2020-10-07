function new = mergeStruct(original, new)
	% Will add fields of new to original if they dont exist in original
	n = fieldnames(original);
	for i = 1:length(n)
		new.(n{i}) = original.(n{i});
	end
end