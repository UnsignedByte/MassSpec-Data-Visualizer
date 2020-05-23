clear all;
addpath('utils')

name = input('Dataset Name:', 's');
file = str2double(input('Number for Base File:', 's'));

wantedMods = splitlines(fileread(fullfile('Results', name, 'Params', 'mods.txt')));
wantedMods = wantedMods(~cellfun('isempty', wantedMods)); %remove empty

if ~isfolder(fullfile('Results', name))
    disp(['Missing HeatMap for ' name '. combinedHM_v2 must be run on the dataset before linear regression analysis can be done.']);
    return;
end

if ~isfolder(fullfile('Results', name, 'Regression'))
    mkdir(fullfile('Results', name, 'Regression'));
end

for i = 1:numel(wantedMods)
    sumdat = readtable(fullfile('Results', name, 'HeatMap', [wantedMods{i} '.csv']));
    sumdat = sumdat(((sumdat.Row_Type-sumdat.Contaminant) == 1),1:end-5); % parse out only classes with no contaminant
    sumdat = sumdat(:,cellfun(@(x) startsWith(x, 'x_OfSpectra_'), sumdat.Properties.VariableNames));
    for j = 1:size(sumdat,2)
        if file == j
            continue;
        end
        f = figure('visible', 'off');
        hold on
            scatter(table2array(sumdat(:,file)),table2array(sumdat(:,j)), 25, linspace(0, 1, size(sumdat,1)), 'filled');
            saveas(f, fullfile('Results', name, 'Regression', [wantedMods{i} '_File_' num2str(file) '_vs_' num2str(j) '_normal.svg']));
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            saveas(f, fullfile('Results', name, 'Regression', [wantedMods{i} '_File_' num2str(file) '_vs_' num2str(j) '_log.svg']));
        hold off
    end
end