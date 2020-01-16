function modlist = getAllMods(summarydat)
    modlistStart = find(contains(summarydat.Var2, '%Fixed and variable modifications:'))+1; %index of start of list of mods
    modlistEnd = find(contains(summarydat.Var2, '% Custom modification text below'))-1; %index of end of list of mods
    modlist = summarydat.Var2(modlistStart:modlistEnd); %get list of mods
    modReg = '(.+)\s\/\s(\+|-)([0-9.]+)\s@\s([A-Za-z.,\s]+)\s\|\s(.+)'; %regex used to parse mod
    modlist = cellfun(@(x) regexp(x, modReg, 'tokens'), modlist); %split mod strings into chunks
    modlist = vertcat(modlist{:}); %convert to 2d
    modlist(:,3) = num2cell(cellfun(@(x) round(str2double(x),3), modlist(:,3))); %truncate masses to 3 decimals
    modlist(:,4) = cellfun(@parseStarts, modlist(:,4), 'UniformOutput', false); 
    modlist(:,6) = arrayfun(@(x) matlab.lang.makeValidName([modlist{x,1} '_' modlist{x,5}]),  1:size(modlist, 1),  'UniformOutput',  false);
    
    function matches = parseStarts(str) %parse list of starting positions into regex matching string
        matches = regexp(str, '(?<=\s|^)([a-zA-Z\s]+)(?=,|$)', 'match'); %match each starting pos
        %protein Nterm is special case as it has a space
        for ii = 1:length(matches)
            matches{ii} = strrep(matches{ii}, ' ', '');% remove spaces
            matches{ii} = strrep(matches{ii}, 'ProteinNTerm', '^(M|-)\.'); %protein nterm (can be M or -)
            matches{ii} = strrep(matches{ii}, 'ProteinCTerm', ' \.-$');% protein cterm
            matches{ii} = strrep(matches{ii}, 'NTerm', '\.'); %replace NTerm with .
            matches{ii} = strrep(matches{ii}, 'CTerm', ' \.'); %Cterm with '<mod>.'
            if isempty(strfind(matches{ii}, ' '))
                matches{ii} = [matches{ii} ' '];
            end
        end
    end
end