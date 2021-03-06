clear all;
% allow access to utils functions
addpath('utils')

% turn off table reading warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% used for debug purposes (saves console output to a file)
diary on

% Initialize parameter parsing script
mex -setup c++
mex CXXFLAGS="\$CXXFLAGS -std=c++17" utils/parseParams.cpp

params = struct;

%% Preference Variables
params.uniqueColumns = [7]; %Columns to take from each dataset (sorted by max and sum of first elem)
params.uniqueCombineFunctions = {{'max', 'nansum'}}; %combined functions across datasets used
params.uniqueClassFunctions = {'nansum'}; %functions to use for each dataset when combining into class
%Will save max and sum as well as values from the datasets
params.singleColumns = [11]; %Columns equal across all datasets (will only take one column)
params.singleClassFunctions = {'max'}; %functions to use for each dataset when combining into class
params.proteins = {};
params.significantProteins = {};
params.mods = {'All'};
params.wantedCol="x_OfSpectra";

params = mergeStruct(parseParams([mfilename '.m']), params);

%select data files
folder = fullfile('Results', params.name, 'Data');
% disp(dir(fullfile(folder, '*.xlsx')));
TempFiles = extractfield(dir(fullfile(folder, '*.xlsx')), 'name');
TempFiles = TempFiles(~startsWith(TempFiles(:),'~$')); %Ignore tempsave files

% [TempFiles, folder] = uigetfile('.xlsx','Choose Data Files', 'Multiselect', 'on');

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

if isfield(params, 'testGroups') && length(params.testGroups) == NumFilesRead
    disp('Using existing test groups...')
    params.testGroups = cell2mat(params.testGroups);
else
    params.testGroups = 1:NumFilesRead;
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

fileidtable = table([1:NumFilesRead]', datasetnames, params.testGroups', 'VariableNames', {'ID', 'Filename', 'Test_Group'});

if ~isfile(fullfile('Results', getResultFolder(TempFiles{1}), 'fileIDs.csv'))
    writetable(fileidtable,fullfile('Results', getResultFolder(TempFiles{1}), 'fileIDs.csv'));
end

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
        % disp(i);
        if kk == 1
            resTables{i} = struct;
            resTables{i}.Summary = table('Size', [1,3], 'VariableTypes', {'uint32', 'string', 'string'}, 'VariableNames', {'SheetNumber', 'Filenames', 'Protein_Name'});
            resTables{i}.Sheets = {};
            resTables{i}.Name = params.proteins{i};
        end
        proteinName = getProteinName(params.proteins{i}, dat.ProteinName, {'dbname', 'genename', 'proteinname'});
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
%             disp(resTables{i}.Summary)
            end
%          disp(sheetID);
%          disp(sheetname);
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
if ~isfolder(fullfile(resfolder, 'Files'))
    mkdir(fullfile(resfolder, 'Files'));
end
if ~isfolder(fullfile(resfolder, 'TestGroups'))
    mkdir(fullfile(resfolder, 'TestGroups'));
end

% params.mods = splitlines(strtrim(fileread(fullfile('Params', 'mods.txt'))));

rt2s = cell(numel(params.mods), 1);

ProteinNames = [];
for i = 1:NumFilesRead
    ProteinNames = [ProteinNames; origpdats{i}.Description];
end
ProteinNames = unique(ProteinNames); % get only unique

parsedProteins = generateParsedProteins(ProteinNames);
% disp(parsedProteins);

disp('Protein names parsed');
toc;

allgroups = unique(params.testGroups);
for kk = 1:numel(params.mods)
    
    wantedMod = params.mods{kk}; %Get list of modifications for heatmap

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

                % dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
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
        % dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero

        for ii = 3:size(pdat,2)
            for jj = 1:size(pdat,1)
                if ~isnumeric(pdat{jj,ii})
                    newvalue = str2double(pdat{jj,ii});
                    if isnan(newvalue) & ~isempty(char(pdat{jj,ii}))
                        warning(['Position (' num2str(jj) ', ' num2str(ii) ') expected number, found "' char(pdat{jj,ii}) '"'])
                    end
                    pdat{jj,ii} = {newvalue};
                elseif isnan(pdat{jj, ii})
                    pdat{jj, ii} = 0; % replace NaN values with 0
                end
            end
        end
        datasets{i} = pdat;
    end
    disp(['Mod ' wantedMod ' loaded.'])
    toc;

    rt2s{kk} = struct;
    rt2s{kk}.Data = getCombined(datasets, cellstr(num2str([1:NumFilesRead]')), parsedProteins, params.uniqueColumns, params.uniqueCombineFunctions, params.uniqueClassFunctions, params.singleColumns, params.singleClassFunctions);
    rt2s{kk}.Name = wantedMod;
    disp(['File data for' wantedMod ' generated.'])
    toc;
    writetable(rt2s{kk}.Data,fullfile(resfolder, 'Files', [wantedMod '.csv']));

    rt2s{kk}.GroupData = rt2s{kk}.Data(:,1:3);
    for uniqueCol = 1:length(params.uniqueColumns)
        tmpcols = table;
        for group = 1:length(allgroups)
            avginds = find(params.testGroups==allgroups(group));
            tmpcols = [ ...
                tmpcols ...
                table( ...
                    mean(rt2s{kk}.Data{:,avginds+3+(uniqueCol-1)*(NumFilesRead+length(params.uniqueCombineFunctions{uniqueCol}))}, 2), ...
                    'VariableNames', {[datasets{1}.Properties.VariableNames{params.uniqueColumns(uniqueCol)} '_' num2str(allgroups(group))]} ...
                    ) ...
                ];
        end
        tmpcols = [
            tmpcols ...
            table( ...
                'Size', [size(rt2s{kk}.Data, 1), 2], ...
                'VariableTypes', {'double', 'double'}, ...
                'VariableNames', cellfun(@(x) [x '_' datasets{1}.Properties.VariableNames{params.uniqueColumns(uniqueCol)}], params.uniqueCombineFunctions{uniqueCol}, 'UniformOutput', false)...
            ) ...
        ];
        startI = 1;
        for endI = 1:size(rt2s{kk}.Data,1)
            for combineFunc = 1:length(params.uniqueCombineFunctions{uniqueCol})
                tmpcols{endI, length(allgroups)+combineFunc} = feval(params.uniqueCombineFunctions{uniqueCol}{combineFunc}, tmpcols{endI,1:length(allgroups)});
            end
        end
        rt2s{kk}.GroupData = [rt2s{kk}.GroupData tmpcols];
    end

    rt2s{kk}.GroupData = [rt2s{kk}.GroupData rt2s{kk}.Data(:,(end-length(params.singleColumns)-1):end)];

    disp(['Testgroup data for ' wantedMod ' done.'])
    toc;

    writetable(rt2s{kk}.GroupData,fullfile(resfolder, 'TestGroups', [wantedMod '.csv']));
end

% Time of completion
completeTime = datestr(now,'dd_mm_yyyy_HH_MM_SS');

Output = struct;
Output.ModMapper = resTables;
Output.HeatMap = rt2s;
% Output.fileIDs = fileidtable;

if ~isfolder(fullfile('Results', getResultFolder(TempFile), 'Raws'))
    mkdir(fullfile('Results', getResultFolder(TempFile), 'Raws'));
end

if ~isfolder(fullfile('Results', getResultFolder(TempFile), 'Significance'))
    mkdir(fullfile('Results', getResultFolder(TempFile), 'Significance'));
end

disp('Generating Significance...')
toc;
for i = 1:length(params.mods)
    sigFile = fullfile('Results', getResultFolder(TempFile), 'Significance', params.mods{i});
    mkdir(sigFile);
    sigData = rt2s{i}.Data(find(rt2s{i}.Data.Row_Type),1:3);
    sigData.provided = cell(size(sigData, 1), 1);
    for j = 1:size(sigData,1)
        strs = regexp([sigData.Protein_Name{j} '/' sigData.Gene_Name{j}], '(?<=^|\/)(?<name>.+?)(?=$|\/)', 'match');
        for k = 1:length(params.significantProteins)
            if any(strcmp(strs,params.significantProteins(k)))
                sigData.provided{j} = 'listed';
            end
        end
    end
    writetable(sigData, fullfile(sigFile, 'raw.csv'));
end

% writetable(rt2s{kk}.(:,1:3))

disp('Saving final files...')
toc;

fid = fopen(fullfile('Results', getResultFolder(TempFile), 'Raws', 'combinedHM.json'), 'w');
saveJSON(fid, Output);
fclose(fid);

if ~isfolder(fullfile('Results', getResultFolder(TempFile), 'Params'))
    mkdir(fullfile('Results', getResultFolder(TempFile), 'Params'));
end

copyfile('params.p', fullfile('Results', getResultFolder(TempFile), 'Params', ['params_' completeTime '.p']))

disp('Done.');
toc;

diary off