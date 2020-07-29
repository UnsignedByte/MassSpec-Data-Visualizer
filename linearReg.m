clear all;
close all;
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

            X = table2array(sumdat(:,file));
            Y = table2array(sumdat(:,j));

            hold on
                set(0, 'CurrentFigure', cfig);
                ax = axes('Parent',cfig);

                if file < j
                    scatter(X,Y, 5, 'k', 'filled');
                elseif file > j
                    % set(ax, 'Color', [1 1 0.0667]);
                    xlim([0 1])
                    ylim([0 1])
                    text(0.5, 0.5, num2str(corr(X, Y)), 'HorizontalAlignment','center','VerticalAlignment','middle')
                end
                % plot(rand(1,(j-1)*nm+file));

%                 set(gca, 'units', 'normalized'); %Just making sure it's normalized
                set(ax,'xtick',[], 'ytick', []);
                
                % disp([file-1 nm-j 1 1]/nm)
                % disp((j-1)*nm+file)
                set(ax, 'position', [file-1 nm-j 1 1]/nm);
            hold off

            if file == j
                continue;
            end
            f = figure('visible', 'off');
            hold on
                % set(0, 'CurrentFigure', f);
                scatter(X,Y, 25, linspace(0, 1, size(sumdat,1)), 'filled');
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['norm_File_' num2str(j) '.svg']));
                set(gca, 'YScale', 'log', 'XScale', 'log')
                saveas(f, fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['log_File_' num2str(j) '.svg']));
            hold off
        end
    end
    saveas(cfig, fullfile(parentF, wantedMods{i}, 'combined.svg'));
end

% cfig = figure()
% nm = 4
% for j = 1:nm
%     for k = 1:nm
%         set(0, 'CurrentFigure', cfig);
%         ax = axes('Parent',cfig);
%         if k < j
%             plot(rand(1,(j-1)*nm+k));
%         elseif k > j
%             disp('hi');
%             set(ax, 'Color', [1 1 0.0667]);
%         end
%         set(ax,'xtick',[], 'ytick', []);
%         disp([k-1 nm-j 1 1]/nm)
%         disp((j-1)*nm+k)
%         set(ax, 'position', [k-1 nm-j 1 1]/nm);
%     end
% end

% f = figure
% set(gca, 'Color', [1 1 0.0667]);
% saveas(f, 'test.svg')
