
function inds = getProteinName(name, proteins, wantedparams) %get all proteins with same gene name
    proteins = parseProteins(proteins);
    inds = {};

    for i = 1:size(proteins,1)
        if any(strcmpi(name, cellfun(@(x) getfieldifexists(proteins{i}, x) , wantedparams, 'UniformOutput', false)))
            inds{end+1} = proteins{i}.fullname;
        end
    end
end