function inds = getProteinName(name, proteins, wantedparam) %get all proteins with same gene name
    proteins = parseProteins(proteins);
    inds = {};
    for i = 1:size(proteins,1)
        if strcmpi(name, proteins{i,wantedparam+1})
            inds{end+1} = proteins{i,1};
        end
    end
end