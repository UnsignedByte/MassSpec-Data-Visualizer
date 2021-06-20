function x = getfieldifexists(s, field, ranks)
	if nargin == 2
		ranks = 1:size(s,1);
	end
	if any(strcmp(fieldnames(s), field))
		x = s.(field)(ranks);
	else
		x = repmat("", length(ranks), 1);
	end
end