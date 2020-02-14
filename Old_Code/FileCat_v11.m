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

for i = 1 : NumFilesRead
    [TempPathStr, TempName, TempExt] = fileparts(TempFile{i});
    TempStruct(i).NameForFile = TempName;
    [tempnum,temptxt,tempraw] = xlsread(TempFile{i},XcelSheet);
    TempStruct(i).Proteins = temptxt(2:length(temptxt(:,2)),2);
    TempStruct(i).SpecCounts = tempnum(:,7); % currently set to spectral counts
    TempStruct(i).AALength = tempnum(:,11); % currently set to AA in Protein
    TempStruct(i).AALength = num2cell(TempStruct(i).AALength);
    TempStruct(i).Max = length(TempStruct(i).Proteins);
    TempStruct(i).SpecNum = num2cell(TempStruct(i).SpecCounts);
    TempStruct(i).CountsOut = cell(TempStruct(i).Max,2);
    TempStruct(i).CountsOut(:,1) = TempStruct(i).Proteins;
    TempStruct(i).CountsOut(:,2) = TempStruct(i).SpecNum;
end

%This is where we remove the extra lines corresponding to isoforms that are
%not uniquely identified. It could, in principle, be disabled if you want
%to give the client all possible assignable names, but they won't have
%values associated with them as this is set up. I think this is a bad idea
%- would rather limit another way.

for i = 1:NumFilesRead
    for j = 1:length(TempStruct(i).CountsOut(:,1))
        if ~isnan(TempStruct(i).CountsOut{j,2})
            TempStruct(i).CountsUsed(j,:) = TempStruct(i).CountsOut(j,:);
            TempStruct(i).AAsUsed(j,1) = TempStruct(i).AALength(j,1);
        end;
    end;
     TempStruct(i).CountsUsed(cellfun('isempty',TempStruct(i).CountsUsed)) = [];
     TempStruct(i).AAsUsed(cellfun('isempty',TempStruct(i).AAsUsed)) = [];
     TempStruct(i).CountsUsed = reshape(TempStruct(i).CountsUsed,[],2);
     TempStruct(i).ActiveProteinNames = TempStruct(i).CountsUsed(:,1);
     
end;
            
%Now I'll generate a string with all of the concatenated protein names and
%limit to only unique names

FileOutNames =[];

for i = 1 : NumFilesRead
    FileOutNames = [FileOutNames; TempStruct(i).ActiveProteinNames];
end;

FileOutUni = unique(FileOutNames);


%There is a string comparison loop that starts at the top of the unique
%list, and goes through each protein name list. If it finds it, it
%transfers the number of 'counts' into the appropriate column. otherwise,
%it leaves in a NaN

FileOutTemp = cell(length(FileOutUni),NumFilesRead);
FileOutTemp(:) = {NaN}; 
FileOutAAs = cell(length(FileOutUni),1);
FileOutCombined = [FileOutUni, FileOutTemp];
toc;

for k = 1 : NumFilesRead
    for i = 1:length(FileOutCombined(:,1))
        for j = 1:length(TempStruct(k).CountsUsed(:,1))
            if strcmp(FileOutCombined(i,1),TempStruct(k).CountsUsed(j,1)) == 1
                FileOutCombined(i,k+1) = TempStruct(k).CountsUsed(j,2);
                FileOutAAs(i,1) = TempStruct(k).AAsUsed(j);
            end;
        end;
    end;
    toc;
end;

%New little bit of code to find the reverses and contaminants and flag them
%for removal.
toc;
ReverseFinder = zeros(length(FileOutCombined(:,1)),1);

%Disabled crapome file, may remove entirely in the next version.
% [num,craptxt,raw] = xlsread(crapome.file);
%this file txt should be looped for each target

for i = 1:length(FileOutCombined(:,1))
    if isempty(strfind(FileOutCombined{i,1},'>Reverse')) == 0
        ReverseFinder(i,1) = 1;
    elseif isempty(strfind(FileOutCombined{i,1},'Common contaminant')) == 0
        ReverseFinder(i,1) = 1;
    elseif isempty(strfind(FileOutCombined{i,1},'eratin, type')) == 0
        ReverseFinder(i,1) = 1;
    elseif isempty(strfind(FileOutCombined{i,1},' desmo')) == 0
        ReverseFinder(i,1) = 1;
    elseif isempty(strfind(FileOutCombined{i,1},'dermi')) == 0
        ReverseFinder(i,1) = 1;
    elseif isempty(strfind(FileOutCombined{i,1},'plak')) == 0
        ReverseFinder(i,1) = 1;
    end;
% If this search is for the Jackson Lab, use the following segment:

%     for j = 1:length(craptxt(:,1))
%         if isempty(strcmp(FileOutCombined{i,1},craptxt{j,1})) == 0
%             reversefinder(i,1) = 1;
%         end;
%     end;
end;

ReverseFinder = num2cell(ReverseFinder);

%Add in the Max, Sum, and Rank columns here

MaxFinder = zeros(length(FileOutCombined(:,1)),1);
SumFinder = zeros(length(FileOutCombined(:,1)),1);

for i = 1 : length(FileOutCombined(:,1))
    MaxFinder(i,1) = max(cell2mat(FileOutCombined(i,2:NumFilesRead+1)),[],'omitnan');
    SumFinder(i,1) = sum(cell2mat(FileOutCombined(i,2:NumFilesRead+1)),[],'omitnan');
end;

MaxFinder = num2cell(MaxFinder);
SumFinder = num2cell(SumFinder);
RankFinder = num2cell(zeros(length(FileOutCombined(:,1)),1));

FileOutFinal = [RankFinder, FileOutCombined, ReverseFinder, MaxFinder, SumFinder, FileOutAAs];
toc;

FileOutFinalDeComma = FileOutFinal;

for i = 1:length(FileOutFinal(:,1))
    if isempty(strfind(FileOutFinal{i,1},',')) == 0
        FileOutFinalDeComma{i,1} = strrep(FileOutFinal{i,1}, ',', '-');
    end;
end;

%Here, I add in the file names in a readable format so they can be column
%headers. Note that these column headers are a tad different than before so
%they are readable, and should be based on the input file names as much as
%possible.

HeaderFileString = {'Rank_Number'; 'Protein_Name'};

for i = 1 : NumFilesRead
    CleanedTempString = strrep(TempStruct(i).NameForFile, '.raw_20', '_');
    CleanedTempString = strrep(CleanedTempString, '_Byonic','');
    CleanedTempString = strrep(CleanedTempString, '_Elite','');
    CleanedTempString = strrep(CleanedTempString, '_control',''); %Remove after I'm done with Janos' dataset
    CleanedTempString = strrep(CleanedTempString, '25fmol_ul_6x5spike_',''); %Remove after I'm done with Janos' dataset
    CleanedTempString = strcat('Out', CleanedTempString);
    HeaderFileString = [HeaderFileString; CleanedTempString];
end;

HeaderFileString = [HeaderFileString; 'Contaminant'; 'Max'; 'Sum'; 'AAs'];

%Now we make the final output file for printing to a csv. we first make a
%table rearrange it appropriately, and then have it fill in the rank number
%values before printing_ it out.

FinalFileOut = cell2table(FileOutFinalDeComma);
FinalFileOut.Properties.VariableNames = HeaderFileString;
FinalFileOut = sortrows(FinalFileOut,{'Contaminant','Max','Sum'},{'ascend','descend','descend'});

for i = 1:length(FileOutFinalDeComma)
    FinalFileOut.Rank_Number(i) = i;
end;
writetable(FinalFileOut,'tableformat.csv');