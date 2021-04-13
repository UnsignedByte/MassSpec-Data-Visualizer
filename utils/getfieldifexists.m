function x = getfieldifexists(s, field)
	if any(strcmp(fieldnames(s), field))
		x = s.(field);
	else
		x = [];
	end
end