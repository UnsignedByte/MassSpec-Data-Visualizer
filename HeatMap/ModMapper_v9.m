clear all

TempFile = uigetfile('.xlsx','Choose Datafile');
dat = readtable(TempFile, 'Sheet', 'Spectra'); %spectra tab
summarydat = readtable(TempFile, 'Sheet', 'Summary'); %summary data

modlistStart = find(contains(summ.Var2, '%Fixed and variable modifications:'))+1; %index of start of list of mods
modlistEnd = find(contains(summ.Var2, '% Custom modification text below'))-1; %index of end of list of mods
modlist = summarydat.Var2(modlistStart:modlistEnd); %get list of mods
modReg = '(.+)\s\/\s(\+|-)([0-9.]+)\s@\s([A-Za-z.,\s]+)\s\|\s(.+)'; %regex used to parse mod
modlist = cellfun(@(x) regexp(x, modReg, 'tokens'), modlist); %split mod strings into chunks
modlist = vertcat(modlist{:}); %convert to 2d
modlist(:,3) = num2cell(cellfun(@(x) round(str2double(x),3), modlist(:,3))); %truncate masses to 3 decimals
modlist(:,4) = cellfun(@

TempFastaFile = uigetfile('.fasta','Choose Datafile');
fastaFile = struct2table(fastaread(TempFastaFile));

proteinRank = 1; %desired rank for protein to analyze
Istart = 0;

for j = 1:size(dat, 1) %get range of rows with specified protein rank
    if dat.ProteinRank(j) == proteinRank && Istart == 0
        Istart = j;
    elseif dat.ProteinRank(j) > proteinRank
        Iend = j-1;
        break;
    end
end

proteinName = dat.ProteinName{Istart(1)};

for i = 1:size(fastaFile,1) %find protein sequence from fasta
    if strcmp(fastaFile.Header{i},proteinName(2:end)) %2:end used to ignore the preceding >
        proteinSequence = fastaFile.Sequence{i};
        break;
    end
end

for i = Istart:Iend
    
end

function reg = parseStarts(str) %parse list of starting positions into regex matching string
    matches = regexp(str, '(?<=\s|^)(Protein NTerm|[a-zA-Z]+)(?=[\s,]|$)', 'match') %match each starting pos
    %protein Nterm is special case as it has a space
end

function [pep, mods] = formatPeptide(pep, possibleMods)
    pep = pep(3:end-2); %strip start and end
    regStr = '[A-Z]\[(?:\+|-)[0-9.]+\]'; %regexp matching mod
    
    [sI, eI] = regexp(pep, regStr);
    match = regexp(pep, regStr, 'match');
    match = cellfun(@(x) [x(2) {x(1), num2str(round(str2num(x(3:end-1)),3))}],match, 'UniformOutput', false); %get match nums
    lens = eI-sI+1; %get lengths of each matched str
    cSum = 0;
    for i = 1:length(sI)
        sI(i) = sI(i)-cSum; 
        sI(i) = sI(i)-cSum; %shift each mod by the value 
    end
end