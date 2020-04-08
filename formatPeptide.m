function [pep, mods] = formatPeptide(pep, modlist)
    function res = isMatch(loc, reg)
%         disp(loc);
%         disp(reg);
        [starts, matches] = regexp(pep, reg, 'start', 'match');
%         disp(starts);
%         disp(matches);
        for kk = 1:length(starts)
            if strfind(matches{kk}, ' ')+starts(kk)-1 == loc
                res = 1;
                return;
            end
        end
        res = 0;
    end
    regStr = '\[(\+|-)([0-9.]+)\]'; %regexp matching mod

    modspecs = regexp(pep, regStr, 'tokens');
    modspecs = vertcat(modspecs{:});
    if numel(modspecs) == 0 % if no mods on peptide
        mods = {};
        return;
    end
    modspecs(:,2) = num2cell(cellfun(@str2double, modspecs(:,2))); %truncate masses to 3 decimals

    pep = regexprep(pep, regStr, ' ');

    locs = strfind(pep, ' '); %get locations of all peptides
%     disp(pep);
%     locs = locs-(1:length(locs));

    mods = cell(size(modspecs,1),2);
    for ii = 1:size(modspecs,1)
        for jj = 1:size(modlist,1)
            if   modlist{jj,3} == modspecs{ii,2} ...%same daltons
              && strcmp(modlist{jj,2}, modspecs{ii,1}) ...%same +/-
              && any(cellfun(@(x) isMatch(locs(ii), x), modlist{jj,7}))
                mods{ii,1} = modlist{jj,1};
                mods{ii,2} = locs(ii)-ii;
            end
        end
    end
    pep = strrep(pep, ' ', '');
end