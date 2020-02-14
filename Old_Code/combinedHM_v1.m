clear all;
%% Preference variables
UniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
UniqueCombineFunctions = {{@max, @nansum}}; %combined functions across datasets used
UniqueClassFunctions = {@nansum}; %functions to use for each dataset when combining into class
%Will save max and sum as well as values from the datasets
SingleColumns = [11]; %Columns equal across all datasets (will only take one column)
SingleClassFunctions = {@max}; %functions to use for each dataset when combining into class

%%

tic;

TempFiles = uigetfile('.xlsx','Choose Data Files', 'Multiselect', 'on');

wantedGenes = splitlines(fileread('wantedProteins.txt'));
if numel(wantedGenes) > 0
    TempFastaFile = uigetfile('.fasta','Choose Fasta File');
    fastaFile = struct2table(fastaread(TempFastaFile));
end


NumFilesRead = length(TempFiles);
if NumFilesRead == 1
    TempFiles = {TempFiles};
end

datasetnames = cell(NumFilesRead,1);

origsumdats = cell(NumFilesRead, 1);
origsdats = cell(NumFilesRead, 1);
origpdats = cell(NumFilesRead, 1);
modlists = cell(NumFilesRead, 1);

for i = 1 : NumFilesRead
    [~, TempName, ~] = fileparts(TempFiles{i});
    datasetnames{i} = TempName;
    origsdats{i} = readtable(TempFiles{i}, 'Sheet', 'Spectra');
    origpdats{i} = readtable(TempFiles{i}, 'Sheet', 'Proteins');
    origsumdats{i} = readtable(TempFiles{i}, 'Sheet', 'Summary');
    modlists{i} = getAllMods(origsumdats{i});
end

for kk = 1:NumFilesRead
    TempFile = TempFiles{kk};
    
    resfolder = fullfile('Results', getResultFolder(TempFile), 'ModMapper');

    if ~isfolder(resfolder)
        mkdir(resfolder);
    end

    dat = origsdats{kk}; %spectra tab
    summarydat = origsumdats{kk}; %summary data
    proteindat = origpdats{kk}; %protein data
    for i = 1:length(wantedGenes)
        proteinName = getProteinName(wantedGenes{i}, proteindat.Description, 3);
        if isempty(proteinName)
            continue;
        end
        proteinName = proteinName{1};
        disp(['Processing Gene ' wantedGenes{i}]);
        toc;
        getPeptideMap(proteinName, dat, summarydat, fastaFile, fullfile(resfolder, wantedGenes{i}), makeValidSheetName(getResultFile(TempFile)));
    end
    % proteinName = '>sp|P14873|MAP1B_MOUSE Microtubule-associated protein 1B OS=Mus musculus OX=10090 GN=Map1b PE=1 SV=2'; %desired rank for protein to analyze
end

resfolder = fullfile('Results', getResultFolder(TempFiles{1}), 'HeatMap');

if ~isfolder(resfolder)
    mkdir(resfolder);
end

wantedMods = splitlines(fileread('wantedMods.txt'));
for j = 1:numel(wantedMods)
    
    wantedMod = wantedMods{j};

    %Building the main data structure here. All data is read in and retained.
    %This should make it flexible to do alternative readouts on the fly. Should
    %be scalable to any number of input file

    datasets = cell(NumFilesRead, 1);
    for i = 1 : NumFilesRead
        pdat = origpdats{i};
        if ~strcmp(wantedMod, 'All')
            sdat = origsdats{i};
            modlist = modlists{i};
            modlist = modlist(strcmp(modlist(:,6),wantedMod),:);

        %    dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
            TempStruct(i).sdat = sdat;
            TempStruct(i).modMatches = cellfun(@(x) contains(x, ['[' modlist{1,2} num2str(modlist{1,3}) ']']),TempStruct(i).sdat.Peptide_ProteinMetricsConfidential_);

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
        end
        datasets{i} = pdat;
    end
    toc;

    getCombined(datasets, datasetnames, UniqueColumns, UniqueCombineFunctions, UniqueClassFunctions, SingleColumns, SingleClassFunctions, fullfile(resfolder, wantedMod));
end