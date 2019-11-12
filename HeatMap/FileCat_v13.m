%NEW VERSION
%
% Updated to Add in AA counts
%
%
% 20170517 UPDATING TO A CLEANER OUTPUT
%
%This version is based on FileCat_v6. The goal is to tighten up the script,
% make it scalable to any number of possible files (up to 100 files), and
%extract data directly from the .xlsx sheets. Also, to properly label the
%output columns. It should serve as a starting place for the transition
%script to R which Kratika will build for us in a few weeks.
%
clear all

%% Preference variables
UniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
UniqueCombineFunctions = {{@max, @nansum}}; %combined functions across datasets used
UniqueClassFunctions = {@nansum}; %functions to use for each dataset when combining into class
%Will save max and sum as well as values from the datasets
SingleColumns = [11]; %Columns equal across all datasets (will only take one column)
SingleClassFunctions = {@nansum}; %functions to use for each dataset when combining into class

%%

%This is not active right now
%[crapraw, craptxt, crapnum] = xlsread('Jacksoncrapome.xlsx');

%This should allow us to call
XcelSheet = 'Proteins';

tic;
%THIS LIST NO LONGER ACCURATE, NEED TO REBUILD 20170517
%num(:,1)  : Log Prob
%num(:,2)  : Best Prob
%num(:,3)  : Best Score
%num(:,4)  : Spectral Count NOW AT num(:,4)
%num(:,5)  : Unique Peptides
%num(:,6)  : Modified Peptides
%num(:,7) : Coverage
%num(:,8) : # of AAs in Protein
%num(:,9)  : Spectral Intensity

%Have a Dialog Box to ID how many files are being searched. Will later add
%the ability to define outputs here as well, for instance NSAF, or unique
%only
NumFilesRead = inputdlg('How Many .xlsx Files?'); 
NumFilesRead = str2num(NumFilesRead{1}); 

%This is where you set the number of files to be read. Currently set up to
%allow you to select all the files at once for time savings. Added a
%warning for cases where you give it a mismatched number of files and
%selected number of files.

TempFile = uigetfile('.xlsx','Choose First Datafile','MultiSelect','on');

NumFilesSelected = length(TempFile);

if NumFilesRead ~= length(TempFile)
    WarningValue = NumFilesSelected - NumFilesRead;
    TempWarning = ['FILE NUMBER MISMATCH: ', num2str(WarningValue)];
    disp(TempWarning)
end;

%Building the main data structure here. All data is read in and retained.
%This should make it flexible to do alternative readouts on the fly. Should
%be scalable to any number of input files.

ProteinNames = []; %all protein names

for i = 1 : NumFilesRead
    [TempPathStr, TempName, TempExt] = fileparts(TempFile{i});
    TempStruct(i).NameForFile = TempName;
    dat = readtable(TempFile{i}, 'Sheet', XcelSheet);
    dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
    TempStruct(i).dat = dat;
    ProteinNames = [ProteinNames; dat.Description];
end

ProteinNames = unique(ProteinNames); %get only unique

%Get all unique proteins over all datasets
ProteinNames(:,2) = num2cell(1:length(ProteinNames)); %Maps protein name to a unique ID

%Match ids if proteins are in groups
for i = 1:NumFilesRead
    for j = 1:TempStruct(i).dat.ProteinRank(end) %loop through all ranks in dataset
        %find ids of a "class" of proteins (proteins with the same rank)
        ids = cell2mat(cellfun(@(x) find(ismember(ProteinNames(:,1),x)), TempStruct(i).dat.Description(find(TempStruct(i).dat.ProteinRank==j)), 'UniformOutput', false));
        locs = [];
        for k = 1:numel(ids)
            locs = [locs find([ProteinNames{:,2}]==ids(k))];
        end
        ProteinNames(locs,2) = {min(cell2mat(ProteinNames(locs,2)))};
    end
end

ProteinNames = sortrows(ProteinNames,2); %sort names by id

toc;

totUniqueFuncs = zeros(length(UniqueColumns),1);
totUniqueFuncs(1) = length(UniqueCombineFunctions{1});
for i = 2:length(UniqueColumns)
    totUniqueFuncs(i) = totUniqueFuncs(i-1)+length(UniqueCombineFunctions{i});
end

Result = NaN(length(ProteinNames),NumFilesRead*length(UniqueColumns)+totUniqueFuncs(end)+length(SingleColumns)+2);
ClassResult = [];

realI = 1;
j = 1;
for i = 1:ProteinNames{end,2}
    startJ = j;
    while j <= size(ProteinNames, 1) && ProteinNames{j,2} == i
        for kk = 1:NumFilesRead
            ind = find(ismember(TempStruct(kk).dat.Description, ProteinNames{j,1})); %find index of peptide
            if ~isempty(ind)
                for k = 1:length(UniqueColumns)
                    Result(j,(k-1)*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+kk) = table2array(TempStruct(kk).dat(ind,UniqueColumns(k))); %paste data into result
                end
                %unique columns should be equal across all datasets, so only
                %run once
                if kk == 1
                    for k = 1:length(SingleColumns)
                        Result(j,length(SingleColumns)*NumFilesRead+totUniqueFuncs(end)+k) = table2array(TempStruct(1).dat(ind,SingleColumns(k))); %paste data into result
                    end
                end
            end
        end
        % add combined columns
        for k = 1:length(UniqueColumns)
            currSliceInd = k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k});
            for l = 1:length(UniqueCombineFunctions{k})
                Result(j, currSliceInd+l)=UniqueCombineFunctions{k}{l}(Result(j, currSliceInd-NumFilesRead+1:currSliceInd));
            end
        end
        
        % code to find the reverses and contaminants and flag them
        %for removal.
        Result(j, end-1) = 0;
        if       ~(isempty(strfind(ProteinNames{j,1},'>Reverse')) ...
                && isempty(strfind(ProteinNames{j,1},'Common contaminant')) ...
                && isempty(strfind(ProteinNames{j,1},'eratin, type')) ...
                && isempty(strfind(ProteinNames{j,1},' desmo')) ...
                && isempty(strfind(ProteinNames{j,1},'dermi')) ...
                && isempty(strfind(ProteinNames{j,1},'plak')))
            Result(j, end-1) = 1;
        end
        Result(j, end) = i;
        j = j + 1;
    end
    j = j - 1;
    if j >= startJ
        ClassResult(end+1,:) = NaN(1, size(Result,2));
        for k = 1:length(UniqueColumns)
            currSliceInd = (k-1)*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k});
            for kk = 1:NumFilesRead
                ClassResult(end, currSliceInd+kk) = UniqueClassFunctions{k}(Result(startJ:j,currSliceInd+kk));
            end
            Slice = ClassResult(end,currSliceInd+1:currSliceInd+NumFilesRead);
            for l = 1:length(UniqueCombineFunctions{k})
                ClassResult(end, k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+l)=UniqueCombineFunctions{k}{l}(Slice);
            end
        end
        for k = 1:length(SingleColumns)
            currInd = length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k;
            ClassResult(end,currInd) = SingleClassFunctions{k}(Result(startJ:j,currInd));
        end
        ClassResult(end,end-1) = any(Result(startJ:j, end-1)); %add flags
        ClassResult(end,end) = i;
        realI = realI + 1;
    end
    j = j + 1;
end

sortOrd = [size(ClassResult,2)-1, NumFilesRead+1, NumFilesRead+2];
sortDirs = {'ascend', 'descend', 'descend'};
ClassResult = sortrows(ClassResult, sortOrd, sortDirs);

Result = num2cell(Result);
Result = [Result(:,end) ProteinNames(:, 1) Result(:,1:end-1)];
ClassResult = num2cell(ClassResult);
ClassResult = [ClassResult(:,end) cell(size(ClassResult,1),1) ClassResult(:,1:end-1)];

CombinedRes = cell(size(ClassResult,1)+size(Result,1), size(Result,2)+1);
realI = 1;

for i = 1:size(ClassResult,1)
    CombinedRes(realI,:) = [i ClassResult(i,2:end) 1];
    found = Result(cell2mat(Result(:,1))==ClassResult{i,1},:);
    found = sortrows(found, sortOrd, sortDirs);
    found = found(:,2:end);
    CombinedRes(realI+1:realI+size(found,1),2:end-1) = found;
    CombinedRes(realI+1:realI+size(found,1),1) = {i};
    CombinedRes(realI+1:realI+size(found,1),end) = {0};
    realI = realI+size(found,1)+1;
end
%currently sorts by max then sum of spectra

toc;
% If this search is for the Jackson Lab, use the following segment:

%     for j = 1:length(craptxt(:,1))
%         if isempty(strcmp(FileOutCombined{i,1},craptxt{j,1})) == 0
%             reversefinder(i,1) = 1;
%         end;
%     end;
% end
toc;

% FileOutFinalClasses = cell(FileOutFinal{end,1},size(FileOutFinal,2)); %get header classes
% 
% j = 1;
% for i = 1:size(FileOutFinalClasses,1)
%     chunk = {};
%     while FileOutFinal{j,1} == i
%         chunk(end+1, :) = FileOutFinal(j,:);
%         j = j + 1;
%     end
%     class = cell(1,size(chunk,2));
%     for ll = 1:length(UniqueClassFunctions)
%         for k = 1:NumFilesRead
%             
%         end
%     end
%     for k = 2:size(chunk,2)
%         class{k} = UniqueClassFunctions{
%     end
% end

FileOutFinalDeComma = CombinedRes;

for i = 1:size(CombinedRes, 1)
    if isempty(strfind(CombinedRes{i,2},',')) == 0
        FileOutFinalDeComma{i,2} = strrep(CombinedRes{i,2}, ',', '-');
    end
end

%Here, I add in the file names in a readable format so they can be column
%headers. Note that these column headers are a tad different than before so
%they are readable, and should be based on the input file names as much as
%possible.

HeaderFileString = cell(size(CombinedRes,2), 1);
HeaderFileString(1:2) = {'Rank_Number'; 'Protein_Name'};
HeaderFileString{end-1} = 'Containment';
HeaderFileString{end} = 'Row_Type';

for i = 1 : NumFilesRead
    
    CleanedFileName = strrep(TempStruct(i).NameForFile, '.raw_20', '_');
    CleanedFileName = strrep(CleanedFileName, '_Byonic','');
    CleanedFileName = strrep(CleanedFileName, '_Elite','');
    CleanedFileName = strrep(CleanedFileName, '_control',''); %Remove after I'm done with Janos' dataset
    CleanedFileName = strrep(CleanedFileName, '25fmol_ul_6x5spike_',''); %Remove after I'm done with Janos' dataset
    %label by dataset
    for j = 1:length(UniqueColumns)
        HeaderFileString{2+(j-1)*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+i} = [TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)} '_' CleanedFileName];
    end
end

%label single columns
for j = 1:length(SingleColumns)
    HeaderFileString{2+length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+j} = TempStruct(1).dat.Properties.VariableNames{SingleColumns(j)};
end

%label max and sum
for j = 1:length(UniqueColumns)
    for jj = 1:length(UniqueCombineFunctions{j})
        HeaderFileString{2+j*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+jj} = [func2str(UniqueCombineFunctions{j}{jj}) '_' TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)}];
    end
end

%Now we make the final output file for printing to a csv. we first make a
%table rearrange it appropriately, and then have it fill in the rank number
%values before printing_ it out.

FinalFileOut = cell2table(FileOutFinalDeComma);
FinalFileOut.Properties.VariableNames = HeaderFileString;

filename = 'tableformat';
writetable(FinalFileOut,[filename '.csv']);
fid = fopen([filename '.json'], 'w');
fprintf(fid, jsonencode(FinalFileOut,'ConvertInfAndNaN', false));
fclose(fid);
toc;