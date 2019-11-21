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
UniqueColumns = [7, 10]; %Columns to take from each dataset (sorted by max and sum of first elem)
UniqueCombineFunctions = {{@max, @nansum}, {@max, @nanmean}}; %combined functions across datasets used
UniqueClassFunctions = {@nansum, @nanmean}; %functions to use for each dataset when combining into class
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

%This is where you set the number of files to be read. Currently set up to
%allow you to select all the files at once for time savings. Added a
%warning for cases where you give it a mismatched number of files and
%selected number of files.

TempFile = uigetfile('.xlsx','Choose First Datafile','MultiSelect','on');

NumFilesRead = length(TempFile);

datasets = cell(NumFilesRead,1);
datasetnames = cell(NumFilesRead,1);

for i = 1 : NumFilesRead
    [~, TempName, ~] = fileparts(TempFile{i});
    datasetnames{i} = TempName;
    datasets{i} = readtable(TempFile{i}, 'Sheet', XcelSheet);
end

getCombined(datasets, datasetnames, UniqueColumns, UniqueCombineFunctions, UniqueClassFunctions, SingleColumns, SingleClassFunctions);
