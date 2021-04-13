function ret = getfieldFromStructCell(s, f)
	ret = cellfun(@(x) getfieldifexists(x, f), s, 'UniformOutput', false);
end