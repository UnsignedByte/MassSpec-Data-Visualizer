clear all;
addpath('utils');

name = input('Dataset Name:', 's');
% get id of base file

% read parameters
pfolder = fullfile('Results', name, 'Params');
allParams = getLatestParams(pfolder);
wantedMods = splitlines(fileread(fullfile(pfolder, allParams('mods'))));
wantedMods = wantedMods(~cellfun('isempty', wantedMods)); %remove empty

if ~isfolder(fullfile('Results', name))
    disp(['Missing HeatMap for ' name '. combinedHM_v2 must be run on the dataset before linear regression analysis can be done.']);
    return;
end

if ~isfolder(fullfile('Results', name, 'Regression'))
    mkdir(fullfile('Results', name, 'Regression'));
end

parentF = fullfile('Results', name, 'Regression');

for i = 1:numel(wantedMods)
    % read in heatmap data
    sumdat = readtable(fullfile('Results', name, 'HeatMap', [wantedMods{i} '.csv']));
    % parse out only classes with no contaminant
    sumdat = sumdat(((sumdat.Row_Type-sumdat.Contaminant) == 1),1:end-5);
    % save only data columns
    sumdat = sumdat(:,cellfun(@(x) startsWith(x, 'x_OfSpectra_'), sumdat.Properties.VariableNames));
    mkdir(fullfile(parentF, wantedMods{i}));
    for file = 1:size(sumdat, 2) % loop through base files
        mkdir(fullfile(parentF, wantedMods{i}), ['File_' num2str(file)]);
        for j = 1:size(sumdat,2)
            if file == j
                continue;
            end
            f = figure('visible', 'off');
            hold on
                scatter(table2array(sumdat(:,file)),table2array(sumdat(:,j)), 25, linspace(0, 1, size(sumdat,1)), 'filled');
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['norm_File_' num2str(j) '.svg']));
                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['log_File_' num2str(j) '.svg']));
            hold off
        end
    end
end