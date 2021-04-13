function new = mergeStruct(original, new)
	if isempty(new) %isempty(fieldnames(new)) does not work as if all fields are empty dot access would be banned
		new = original;
		return;
	end
	% Will add fields of new to original if they dont exist in original
	n = fieldnames(original);
	for i = 1:length(n)
		new.(n{i}) = original.(n{i});
	end
end