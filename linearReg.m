clear all;
addpath('utils');

% Turn off file exists warning
warning('OFF', 'MATLAB:mkdir:DirectoryExists')

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
    disp(['Parsing mod ' wantedMods{i}]);
    cfig = figure('visible', 'off');
    nm = size(sumdat, 2);
    for file = 1:nm % loop through base files
        disp(['Parsing file ' num2str(file)])
        mkdir(fullfile(parentF, wantedMods{i}), ['File_' num2str(file)]);
        for j = 1:nm
            if file == j
                continue;
            end

            X = table2array(sumdat(:,file));
            Y = table2array(sumdat(:,j));

            hold on
                set(0, 'CurrentFigure', cfig);
                subplot(nm, nm, (j-1)*nm+file);
                scatter(X,Y, 10, linspace(0, 1, size(sumdat,1)), 'filled');
                
                set(gca, 'units', 'normalized'); %Just making sure it's normalized
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);

                Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                                 %[Left Bottom Right Top] spacing
                NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
                set(gca, 'Position', NewPos);
            hold off

            f = figure('visible', 'off');
            hold on
                scatter(X,Y, 25, linspace(0, 1, size(sumdat,1)), 'filled');
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['norm_File_' num2str(j) '.svg']));
                set(gca, 'YScale', 'log')
                set(gca, 'XScale', 'log')
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['log_File_' num2str(j) '.svg']));
            hold off
        end
    end
    saveas(cfig, fullfile(parentF, wantedMods{i}, 'combined.svg'));
end