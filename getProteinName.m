function inds = getProteinName(name, proteins, wantedparam) %get all proteins with same gene name
    function parsed = parseProteins() %return info about each protein in a list of protein names
        regStr = '^>?(.+?)\|(.+?)\|(.+?)_(.+?)\s((?:.+?)+?)\sOS=((?:.|\s|)+?)\sOX=(\d+?)(\sGN=.+?)?(\sPE=\d+)?(\sSV=\d+)?$';
        parsed = cellfun(@(x) [x, regexp(x, regStr, 'tokens'), {cell(1,10)}], proteins, 'UniformOutput', false);
        parsed = cellfun(@(x) [x{1}, x{2}], parsed, 'UniformOutput', false);
        parsed = vertcat(parsed{:});
        parsed(:, 9:end) = cellfun(@(x) x(5:end), parsed(:, 9:end), 'UniformOutput', false); %remove SV=, etc 
    end

    proteins = parseProteins();
    inds = {};
    for i = 1:size(proteins,1)
        if strcmpi(name, proteins{i,wantedparam+1})
            inds{end+1} = proteins{i,1};
        end
    end
end