clear all;
% allow access to utils functions
addpath('utils')

% Add parse parameters thing
mex -setup c++
mex -v CXXFLAGS="\$CXXFLAGS -std=c++17" utils/parseParams.cpp

params = struct;

% turn off table reading warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% used for debug purposes (saves console output to a file)
diary on

%% Preference Variables
params.UniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
params.UniqueCombineFunctions = {{@max, @nansum}}; %combined functions across datasets used
params.UniqueClassFunctions = {@nansum}; %functions to use for each dataset when combining into class
%Will save max and sum as well as values from the datasets
params.SingleColumns = [11]; %Columns equal across all datasets (will only take one column)
params.SingleClassFunctions = {@max}; %functions to use for each dataset when combining into class
params.proteins = {"All"};
params.mods = {};

params = mergeStruct(parseParams([mfilename '.m']), params);

%select data files
[TempFiles, folder] = uigetfile('.xlsx','Choose Data Files', 'Multiselect', 'on');

% combine folder with basename to get full path (allows selection of files anywhere)
TFPath = fullfile(folder, TempFiles);

%read protein param file
% params.proteins = splitlines(strtrim(fileread(fullfile('Params', 'proteins.txt'))));
% params.proteins = params.proteins(~cellfun('isempty',params.proteins));
if numel(params.proteins) > 0
    % Read in fasta data if genes are wanted (modmapper)
    [baseName, folder] = uigetfile('.fasta','Choose Fasta File');
    fastaFile = struct2table(fastaread(fullfile(folder, baseName)));
end

% Start timing
tic;

if isa(TempFiles, 'char')
    TempFiles = {TempFiles};
end
NumFilesRead = length(TempFiles);

if isfield(params, 'test_groups') && length(params.test_groups) == NumFilesRead
    params.test_groups = cell2mat(params.test_groups);
else
    params.test_groups = 1:NumFilesRead;
end

datasetnames = cell(NumFilesRead,1);

origsumdats = cell(NumFilesRead, 1);
origsdats = cell(NumFilesRead, 1);
origpdats = cell(NumFilesRead, 1);
modlists = cell(NumFilesRead, 1);

for i = 1 : NumFilesRead
    [~, TempName, ~] = fileparts(TempFiles{i});
    datasetnames{i} = cleanFileName(TempName);
    origsdats{i} = readtable(TFPath{i}, 'Sheet', 'Spectra');
    origpdats{i} = readtable(TFPath{i}, 'Sheet', 'Proteins');
    origsumdats{i} = readtable(TFPath{i}, 'Sheet', 'Summary');
    modlists{i} = getAllMods(origsumdats{i});
    disp(['Read File ' num2str(i)]);
    toc;
end

fileidtable = table([1:NumFilesRead]', datasetnames, params.test_groups', 'VariableNames', {'ID', 'Filename', 'Test_Group'});

writetable(fileidtable,fullfile('Results', getResultFolder(TempFiles{1}), 'fileIDs.csv'));

resTables = cell(1, length(params.proteins));

for kk = 1:NumFilesRead
    TempFile = TempFiles{kk};
    
    resfolder = fullfile('Results', getResultFolder(TempFile), 'ModMapper');

    if ~isfolder(resfolder)
        mkdir(resfolder);
    end

    dat = origsdats{kk}; %spectra tab
    summarydat = origsumdats{kk}; %summary data
    for i = 1:length(params.proteins)
        if kk == 1
            resTables{i} = struct;
            resTables{i}.Summary = table('Size', [1,3], 'VariableTypes', {'uint32', 'string', 'string'}, 'VariableNames', {'SheetNumber', 'Filenames', 'Protein_Name'});
            resTables{i}.Sheets = {};
            resTables{i}.Name = params.proteins{i};
        end
        proteinName = getProteinName(params.proteins{i}, dat.ProteinName, 3);
        if isempty(proteinName)
            continue;
        end
        proteinName = proteinName{1};
        filename = fullfile(resfolder, params.proteins{i});
        sheetname = makeValidSheetName(getResultFile(TempFile));
        disp(['Processing Gene ' params.proteins{i} ' for file ' sheetname]);
        toc;
        resTable = getPeptideMap(proteinName, dat, summarydat, fastaFile, sheetname);
        if numel(resTable)==0
            continue;
        end
        if ~isempty(resTables{i}.Sheets)
            if ismember(proteinName, resTables{i}.Summary.Protein_Name)
                sheetID = resTables{i}.Summary.SheetNumber(contains(resTables{i}.Summary.Protein_Name,proteinName));
                resTables{i}.Summary.Filenames{sheetID} = [resTables{i}.Summary.Filenames{sheetID} ', ' sheetname];
                resTable = [resTables{i}.Sheets{sheetID} resTable(:,3:end)];
            else
                sheetID = size(resTables{i}.Summary, 1)+1;
                resTables{i}.Summary = [resTables{i}.Summary;{sheetID, sheetname, proteinName}];
%                 disp(resTables{i}.Summary)
            end
%             disp(sheetID);
%             disp(sheetname);
        else
            resTables{i}.Summary.SheetNumber(1) = 1;
            resTables{i}.Summary.Filenames{1} = sheetname;
            resTables{i}.Summary.Protein_Name{1} = proteinName;
            sheetID = 1;
        end
        resTables{i}.Sheets{sheetID} = resTable;
    end
    % proteinName = '>sp|P14873|MAP1B_MOUSE Microtubule-associated protein 1B OS=Mus musculus OX=10090 GN=Map1b PE=1 SV=2'; %desired rank for protein to analyze
end

for i = 1:length(params.proteins)
    filename = fullfile(resfolder, params.proteins{i});
    writetable(resTables{i}.Summary,[filename '.xlsx'], 'Sheet', 'Summary');
    for j = 1:size(resTables{i}.Sheets, 2)
        writetable(resTables{i}.Sheets{j},[filename '.xlsx'], 'Sheet', num2str(j));
    end
end
% fid = fopen([fullfile(resfolder, 'ModMapper') '.json'], 'w');
% fprintf(fid, jsonencode(resTables,'ConvertInfAndNaN', false));
% fclose(fid);

resfolder = fullfile('Results', getResultFolder(TempFiles{1}), 'HeatMap');

if ~isfolder(resfolder)
    mkdir(resfolder);
end

% params.mods = splitlines(strtrim(fileread(fullfile('Params', 'mods.txt'))));

rt2s = cell(numel(params.mods), 1);

for kk = 1:numel(params.mods)
    
    wantedMod = params.mods{kk};

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
            
            TempStruct(i).sdat = sdat;
            if numel(modlist) == 0
                TempStruct(i).modMatches = zeros(size(TempStruct(i).sdat.Peptide_ProteinMetricsConfidential_));
            else

            %    dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
                TempStruct(i).modMatches = cellfun(@(x) contains(x, ['[' modlist{1,2} num2str(modlist{1,3}) ']']),TempStruct(i).sdat.Peptide_ProteinMetricsConfidential_);
            end

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
    disp(['Mod ' wantedMod ' finished.'])
    toc;

    rt2s{kk} = struct;
    rt2s{kk}.Data = getCombined(datasets, cellstr(num2str([1:NumFilesRead]')), params.UniqueColumns, params.UniqueCombineFunctions, params.UniqueClassFunctions, params.SingleColumns, params.SingleClassFunctions);
    rt2s{kk}.Name = wantedMod;
    writetable(rt2s{kk}.Data,fullfile(resfolder, [wantedMod '.csv']));
end

disp('Saving final files...')

% Time of completion
completeTime = datestr(now,'dd-mm-yyyy_HH:MM:SS');

Output = struct;
Output.ModMapper = resTables;
Output.HeatMap = rt2s;
% Output.fileIDs = fileidtable;

if ~isfolder(fullfile('Results', getResultFolder(TempFile), 'Raws'))
    mkdir(fullfile('Results', getResultFolder(TempFile), 'Raws'));
end

fid = fopen(fullfile('Results', getResultFolder(TempFile), 'Raws', 'combinedHM.json'), 'w');
saveJSON(fid, Output);
fclose(fid);

if ~isfolder(fullfile('Results', getResultFolder(TempFile), 'Params'))
    mkdir(fullfile('Results', getResultFolder(TempFile), 'Params'));
end

fid = fopen(fullfile('Results', getResultFolder(TempFile), 'Params', ['proteins_' completeTime '.txt']),'w');
for i = 1 : numel(params.proteins)
    fprintf(fid, '%s\n', params.proteins{i});
end
fclose(fid);

fid = fopen(fullfile('Results', getResultFolder(TempFile), 'Params', ['mods_' completeTime '.txt']),'w');
for i = 1 : numel(params.mods)
    fprintf(fid, '%s\n', params.mods{i});
end
fclose(fid);

disp('Done.');
toc;

diary off