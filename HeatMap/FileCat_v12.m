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

UniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
%Will save max and sum as well as values from the datasets
SingleColumns = [11]; %Columns equal across all datasets (will only take one column)

Result = NaN(length(ProteinNames),(NumFilesRead+2)*length(UniqueColumns)+length(SingleColumns));

i = 1;
j = 1;
while i <= ProteinNames{end,2}
    while j <= size(ProteinNames, 1) && ProteinNames{j,2} == i
        for kk = 1:NumFilesRead
            ind = find(ismember(TempStruct(kk).dat.Description, ProteinNames{j,1})); %find index of peptide
            if ~isempty(ind)
                for k = 1:length(UniqueColumns)
                    Result(j,(k-1)*(NumFilesRead+2)+kk) = table2array(TempStruct(kk).dat(ind,UniqueColumns(k))); %paste data into result
                end
            
                %unique columns should be equal across all datasets, so only
                %run once
                if kk == 1
                    for k = 1:length(SingleColumns)
                        Result(j,length(SingleColumns)*(NumFilesRead+2)+k) = table2array(TempStruct(1).dat(ind,SingleColumns(k))); %paste data into result
                    end
                end
            end
        end
        % add sum and max columns
        for k = 1:length(UniqueColumns)
            slice = Result(j,(k-1)*(NumFilesRead+2)+1:k*(NumFilesRead+2)-2);
            if ~all(isnan(slice)) %check if all elems are NaN, if so, skip
                Result(j, k*(NumFilesRead+2)-1)=max(slice);
                Result(j, k*(NumFilesRead+2))=nansum(slice);
            end
        end
        j = j + 1;
    end
    i = i + 1;
end

%New little bit of code to find the reverses and contaminants and flag them
%for removal.
toc;
ReverseFinder = zeros(size(Result,1),1);

%Disabled crapome file, may remove entirely in the next version.
% [num,craptxt,raw] = xlsread(crapome.file);
%this file txt should be looped for each target

for i = 1:size(Result,1)
    if ~isempty(strfind(ProteinNames{i,1},'>Reverse'))
        ReverseFinder(i,1) = 1;
    elseif ~isempty(strfind(ProteinNames{i,1},'Common contaminant'))
        ReverseFinder(i,1) = 1;
    elseif ~isempty(strfind(ProteinNames{i,1},'eratin, type'))
        ReverseFinder(i,1) = 1;
    elseif ~isempty(strfind(ProteinNames{i,1},' desmo'))
        ReverseFinder(i,1) = 1;
    elseif ~isempty(strfind(ProteinNames{i,1},'dermi'))
        ReverseFinder(i,1) = 1;
    elseif ~isempty(strfind(ProteinNames{i,1},'plak'))
        ReverseFinder(i,1) = 1;
    end
% If this search is for the Jackson Lab, use the following segment:

%     for j = 1:length(craptxt(:,1))
%         if isempty(strcmp(FileOutCombined{i,1},craptxt{j,1})) == 0
%             reversefinder(i,1) = 1;
%         end;
%     end;
end

ReverseFinder = num2cell(ReverseFinder);

Result(isnan(Result))=-Inf;

Result = num2cell(Result);

FileOutFinal = [ProteinNames(:, [2 1]), Result, ReverseFinder];
toc;

FileOutFinalDeComma = FileOutFinal;

for i = 1:size(FileOutFinal, 1)
    if isempty(strfind(FileOutFinal{i,2},',')) == 0
        FileOutFinalDeComma{i,2} = strrep(FileOutFinal{i,2}, ',', '-');
    end
end

%Here, I add in the file names in a readable format so they can be column
%headers. Note that these column headers are a tad different than before so
%they are readable, and should be based on the input file names as much as
%possible.

HeaderFileString = cell(3+size(Result,2), 1);
HeaderFileString(1:2) = {'Rank_Number'; 'Protein_Name'};
HeaderFileString{end} = 'Contaminant';

for i = 1 : NumFilesRead
    
    CleanedFileName = strrep(TempStruct(i).NameForFile, '.raw_20', '_');
    CleanedFileName = strrep(CleanedFileName, '_Byonic','');
    CleanedFileName = strrep(CleanedFileName, '_Elite','');
    CleanedFileName = strrep(CleanedFileName, '_control',''); %Remove after I'm done with Janos' dataset
    CleanedFileName = strrep(CleanedFileName, '25fmol_ul_6x5spike_',''); %Remove after I'm done with Janos' dataset
    %label by dataset
    for j = 1:length(UniqueColumns)
        HeaderFileString{2+(j-1)*(NumFilesRead+2)+i} = [TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)} '_' CleanedFileName];
    end
end

%label single columns
for j = 1:length(SingleColumns)
    HeaderFileString{2+length(UniqueColumns)*(NumFilesRead+2)+j} = TempStruct(1).dat.Properties.VariableNames{SingleColumns(j)};
end

%label max and sum
for j = 1:length(UniqueColumns)
    HeaderFileString{2+j*(NumFilesRead+2)-1} = ['Max_' TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)}];
    HeaderFileString{2+j*(NumFilesRead+2)} = ['Sum_' TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)}];
end

%Now we make the final output file for printing to a csv. we first make a
%table rearrange it appropriately, and then have it fill in the rank number
%values before printing_ it out.

FinalFileOut = cell2table(FileOutFinalDeComma);
FinalFileOut.Properties.VariableNames = HeaderFileString;

%Sort data by max, sum of first column by grouping
UniqueIDS = unique([ProteinNames{:,2}]);
maxes = zeros(length(UniqueIDS), 3);
maxes(:,1) = UniqueIDS';

i = 1;
j = 1;
while i <= size(maxes,1)
    cMax = [NaN NaN];
    while j <= size(ProteinNames, 1) && ProteinNames{j,2} == maxes(i,1)
        cMax = max(cMax, table2array(FinalFileOut(j,3+NumFilesRead+(0:1))));
        j = j + 1;
    end
    maxes(i,2:end) = cMax;
    i = i + 1;
end
maxes = sortrows(maxes, [2 3], {'descend', 'descend'});

for i = 1:size(maxes,1)
    FinalFileOut.Rank_Number(FinalFileOut.Rank_Number==maxes(i,1), 1) = i;
end

FirstName = TempStruct(1).dat.Properties.VariableNames{UniqueColumns(1)}; %name of first column
tempArr = table2array(FinalFileOut(:,3+NumFilesRead+(0:1)));
FinalFileOut = sortrows(FinalFileOut,{'Contaminant', 'Rank_Number', ['Max_' FirstName],['Sum_' FirstName]},{'ascend','ascend','descend','descend'});

currI = 0;
currR = 0;

for i = 1:size(FileOutFinalDeComma, 1)
    if currI ~= FinalFileOut.Rank_Number(i)
        currI = FinalFileOut.Rank_Number(i);
        currR = currR + 1;
    end
    FinalFileOut.Rank_Number(i) = currR; %update rank
end
filename = 'tableformat.csv';
writetable(FinalFileOut,filename)
toc;