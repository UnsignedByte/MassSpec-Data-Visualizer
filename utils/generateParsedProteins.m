function ret = generateParsedProteins(list)
	parsed = parseProteins(list);
	names = cellfun(@(x) fieldnames(x), parsed, 'UniformOutput', false);
	names = unique(vertcat(names{:}));
	ret = table('Size', [length(list), length(names)], 'VariableTypes', repmat("string", 1, length(names)), 'VariableNames', names);
	for i = 1:length(list)
		for j = 1:length(names)
			if any(strcmp(fieldnames(parsed{i}), names{j}))
				ret{i,j} = string(parsed{i}.(names{j}));
			end
		end
	end
end