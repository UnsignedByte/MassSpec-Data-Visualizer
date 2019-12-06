clear all;
%% Preference variables
UniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
UniqueCombineFunctions = {{@max, @nansum}}; %combined functions across datasets used
UniqueClassFunctions = {@nansum}; %functions to use for each dataset when combining into class
%Will save max and sum as well as values from the datasets
SingleColumns = [11]; %Columns equal across all datasets (will only take one column)
SingleClassFunctions = {@nansum}; %functions to use for each dataset when combining into class

%%

tic;
wantedMods = {'+80'};
wantedLocs = {'*'}; %locations wanted, * = all. If multiple, list together
%I.E. MS. = M, S, or .


wantedRegs = cell(length(wantedMods),1);
for i = 1:length(wantedMods)
    %get regex match pattern
    wantedRegs{i} = [strrep(regexptranslate('escape', wantedLocs{i}), '\*', '.') '\[' regexptranslate('escape', wantedMods{i}) '\]'];
end

%read datafiles
TempFile = uigetfile('.xlsx','Choose Data Files','MultiSelect','on');

NumFilesRead = length(TempFile);

datasets = cell(NumFilesRead,1);
datasetnames = cell(NumFilesRead,1);

%Building the main data structure here. All data is read in and retained.
%This should make it flexible to do alternative readouts on the fly. Should
%be scalable to any number of input file

for i = 1 : NumFilesRead
    [~, TempName, ~] = fileparts(TempFile{i});
    datasetnames{i} = TempName;
    sdat = readtable(TempFile{i}, 'Sheet', 'Spectra');
    pdat = readtable(TempFile{i}, 'Sheet', 'Proteins');
%    dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
    TempStruct(i).sdat = sdat;
    TempStruct(i).modMatches = zeros(size(sdat,1),length(wantedMods)+1);
    for j = 1:length(wantedMods)
        TempStruct(i).modMatches(:,j) = cellfun(@(x) length(regexp(x,wantedRegs{j})),TempStruct(i).sdat.ModificationType_s_);
    end
    TempStruct(i).modMatches(:,end) = sum(TempStruct(i).modMatches,2);
    
    curRank = TempStruct(i).sdat.ProteinRank(1);
    curI = 1;
    curJ = 1;
    for j = 1:size(sdat,1)
        if TempStruct(i).sdat.ProteinRank(j)~=curRank
            pdat{curI, 7} = nnz(TempStruct(i).modMatches(curJ:j-1,end));
            curRank =  TempStruct(i).sdat.ProteinRank(j);
            while pdat.ProteinRank(curI) ~= curRank
                curI = curI + 1;
            end
            curJ = j;
        end
    end
    datasets{i} = pdat;
end
toc;

getCombined(datasets, datasetnames, UniqueColumns, UniqueCombineFunctions, UniqueClassFunctions, SingleColumns, SingleClassFunctions);