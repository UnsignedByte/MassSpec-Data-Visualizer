clear all

TempFile = uigetfile('.xlsx','Choose Datafile');
dat = readtable(TempFile, 'Sheet', 'Spectra'); %spectra tab
summarydat = readtable(TempFile, 'Sheet', 'Summary'); %summary data

modlistStart = find(contains(summarydat.Var2, '%Fixed and variable modifications:'))+1; %index of start of list of mods
modlistEnd = find(contains(summarydat.Var2, '% Custom modification text below'))-1; %index of end of list of mods
modlist = summarydat.Var2(modlistStart:modlistEnd); %get list of mods
modReg = '(.+)\s\/\s(\+|-)([0-9.]+)\s@\s([A-Za-z.,\s]+)\s\|\s(.+)'; %regex used to parse mod
modlist = cellfun(@(x) regexp(x, modReg, 'tokens'), modlist); %split mod strings into chunks
modlist = vertcat(modlist{:}); %convert to 2d
modlist(:,3) = num2cell(cellfun(@(x) round(str2double(x),3), modlist(:,3))); %truncate masses to 3 decimals
modlist(:,4) = cellfun(@parseStarts, modlist(:,4), 'UniformOutput', false); 
modlist(:,1) = arrayfun(@(x) matlab.lang.makeValidName([modlist{x,1} '_' modlist{x,5}]),  1:size(modlist, 1),  'UniformOutput',  false);

TempFastaFile = uigetfile('.fasta','Choose Datafile');
fastaFile = struct2table(fastaread(TempFastaFile));

proteinName = '>sp|P35802|GPM6A_MOUSE Neuronal membrane glycoprotein M6-a OS=Mus musculus OX=10090 GN=Gpm6a PE=1 SV=1'; %desired name for protein to analyze
Istart = 0;

for j = 1:size(dat, 1) %get range of rows with specified protein rank
    if strcmp(dat.ProteinName(j), proteinName) && Istart == 0
        Istart = j;
    elseif ~strcmp(dat.ProteinName(j), proteinName) && Istart ~= 0
        Iend = j-1;
        break;
    end
end

for i = 1:size(fastaFile,1) %find protein sequence from fasta
    if strcmp(fastaFile.Header{i},proteinName(2:end)) %2:end used to ignore the preceding >
        proteinSequence = fastaFile.Sequence{i};
        break;
    end
end

tablewidth = size(modlist, 1)+3;

alltypes    = cell(tablewidth, 1); %types of the thing
alltypes(:) = {'uint32'};
alltypes{2} = 'char';
allvarnames = cell(tablewidth, 1);
allvarnames(1:3) = {'Position', 'AA', 'Total_Peptides'};
allvarnames(4:end) = modlist(:,1);

resTable = table('Size', [length(proteinSequence), tablewidth],'VariableTypes',alltypes, 'VariableNames', allvarnames);

dqueue = zeros((Iend-Istart+1)*2,3);
peps = cell(2,Iend-Istart+1);

for i = Istart:Iend
    [pep, mods] = formatPeptide(dat.Peptide_ProteinMetricsConfidential_{i}, modlist);
    peps{1,i-Istart+1} = pep;
    peps{2,i-Istart+1} = mods;
    dqueue((i-Istart)*2+(1:2),:) = [dat.StartingPosition(i), i-Istart+1, 1; dat.StartingPosition(i)+length(pep)-4, i-Istart+1, -1];
end
dqueue = sortrows(dqueue);

for i = 1:length(proteinSequence)
    lastSlice = resTable(max(i-1,1),3);
    while dqueue(1,1) == i
        lastSlice.Total_Peptides(1) = lastSlice.Total_Peptides(1)+dqueue(1,3);
        for j = 1:size(peps{2,dqueue(1,2)},1)
%             disp(peps{2,dqueue(1,2)});
            if dqueue(1,3) == 1
                resTable.(peps{2,dqueue(1,2)}{j,1})(max(i+peps{2, dqueue(1,2)}{j,2}-3,1)) = resTable.(peps{2,dqueue(1,2)}{j,1})(max(i+peps{2, dqueue(1,2)}{j,2}-3,1)) + 1;
            end
%             lastSlice.(peps{2,dqueue(1,2)}{j,1})(1) = lastSlice.(peps{2,dqueue(1,2)}{j,1})(1)+dqueue(1,3);
        end
        dqueue(1,:) = [];
        if size(dqueue,1) == 0
            break;
        end
    end
%     disp(lastSlice);
    resTable.Position(i) = i;
    resTable.AA{i} = proteinSequence(i);
    resTable(i,3) = lastSlice;
    if size(dqueue,1) == 0
        break;
    end
end

filename = 'modformat';
writetable(resTable,[filename '.csv']);

function matches = parseStarts(str) %parse list of starting positions into regex matching string
    matches = regexp(str, '(?<=\s|^)([a-zA-Z\s]+)(?=,|$)', 'match'); %match each starting pos
    %protein Nterm is special case as it has a space
    for i = 1:length(matches)
        matches{i} = strrep(matches{i}, ' ', '');% remove spaces
        matches{i} = strrep(matches{i}, 'ProteinNTerm', '^(M|-)\.'); %protein nterm (can be M or -)
        matches{i} = strrep(matches{i}, 'ProteinCTerm', ' \.-$');% protein cterm
        matches{i} = strrep(matches{i}, 'NTerm', '\.'); %replace NTerm with .
        matches{i} = strrep(matches{i}, 'CTerm', ' \.'); %Cterm with '<mod>.'
        if isempty(strfind(matches{i}, ' '))
            matches{i} = [matches{i} ' '];
        end
    end
end

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
    for i = 1:size(modspecs,1)
        for j = 1:size(modlist,1)
            if   modlist{j,3} == modspecs{i,2} ...%same daltons
              && strcmp(modlist{j,2}, modspecs{i,1}) ...%same +/-
              && any(cellfun(@(x) isMatch(locs(i), x), modlist{j,4}))
                mods{i,1} = modlist{j,1};
                mods{i,2} = locs(i)-i;
            end
        end
    end
    pep = strrep(pep, ' ', '');
end