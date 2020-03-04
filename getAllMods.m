function modlist = getAllMods(summarydat)
    modReg = '(.+)\s\/\s(\+|-)([0-9.]+)\s@\s([A-Za-z.,\s\-]+)\s\|\s(.+)'; %regex used to parse mod
    modlistStart = find(contains(summarydat.Var2, '%Fixed and variable modifications:'))+1; %index of start of list of mods
    if (numel(regexp(summarydat.Var2{modlistStart}, modReg, 'tokens')) == 0)
        modlistStart = modlistStart+1;
    end
    modlistEnd = modlistStart-1+find(contains(summarydat.Var2(modlistStart:end), '% Custom modification text below'))-1; %index of end of list of mods
    modlist = summarydat.Var2(modlistStart:modlistEnd); %get list of mods
    
    modlist = cellfun(@(x) regexp(x, modReg, 'tokens'), modlist); %split mod strings into chunks
    modlist = vertcat(modlist{:}); %convert to 2d
    modlist(:,3) = num2cell(cellfun(@(x) round(str2double(x),3), modlist(:,3))); %truncate masses to 3 decimals
    modlist(:,4) = cellfun(@parseStarts, modlist(:,4), 'UniformOutput', false);
    modlist(:,6) = modlist(:,1);
    modlist(:,1) = arrayfun(@(x) matlab.lang.makeValidName([modlist{x,1} '_' modlist{x,5}]),  1:size(modlist, 1),  'UniformOutput',  false);
    
    function matches = parseStarts(str) %parse list of starting positions into regex matching string
        matches = regexp(str, '(?<=\s|^)([a-zA-Z\s\-]+)(?=,|$)', 'match'); %match each starting pos
        %protein Nterm is special case as it has a space
        for ii = 1:length(matches)
            matches{ii} = regexprep(matches{ii}, ' ', '');% remove spaces
            matches{ii} = regexprep(matches{ii}, 'proteinn-?term', '^(M|-)\.', 'preservecase'); %protein nterm (can be M or -)
            matches{ii} = regexprep(matches{ii}, 'proteinc-?term', ' \.-$', 'preservecase');% protein cterm
            matches{ii} = regexprep(matches{ii}, 'n-?term', '\.', 'preservecase'); %replace NTerm with .
            matches{ii} = regexprep(matches{ii}, 'c-?term', ' \.', 'preservecase'); %Cterm with '<mod>.'
            if isempty(strfind(matches{ii}, ' '))
                matches{ii} = [matches{ii} ' '];
            end
        end
    end
end