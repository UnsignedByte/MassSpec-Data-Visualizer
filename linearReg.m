clear all;
close all;
addpath('utils');

% Turn off file exists warning
warning('OFF', 'MATLAB:mkdir:DirectoryExists')

mex -setup c++
mex CXXFLAGS="\$CXXFLAGS -std=c++17" utils/parseParams.cpp

% params = struct;
% params = mergeStruct(parseParams([mfilename '.m']), params);

params = parseParams([mfilename '.m']);

if ~isfield(params, "name")
    params.name = input('Dataset Name:', 's');
end
% get id of base file

% read parameters
pfolder = fullfile('Results', params.name, 'Params');
allParams = getLatestParams(pfolder);
wantedMods = splitlines(fileread(fullfile(pfolder, allParams('mods'))));
wantedMods = wantedMods(~cellfun('isempty', wantedMods)); %remove empty

if ~isfolder(fullfile('Results', params.name))
    disp(['Missing HeatMap for ' params.name '. combinedHM_v2 must be run on the dataset before linear regression analysis can be done.']);
    return;
end

if ~isfolder(fullfile('Results', params.name, 'Regression'))
    mkdir(fullfile('Results', params.name, 'Regression'));
end

parentF = fullfile('Results', params.name, 'Regression');

jsonOut = cell(1, numel(wantedMods));

for i = 1:numel(wantedMods)
    % read in heatmap data
    sumdat = readtable(fullfile('Results', params.name, 'HeatMap', [wantedMods{i} '.csv']));
    % parse out only classes with no contaminant
    sumdat = sumdat(((sumdat.Row_Type-sumdat.Contaminant) == 1),1:end-5);
    % save only data columns
    sumdat = sumdat(:,cellfun(@(x) startsWith(x, 'x_OfSpectra_'), sumdat.Properties.VariableNames)); % # spectra from each file
    mkdir(fullfile(parentF, wantedMods{i}));
    disp(['Parsing mod ' wantedMods{i}]);
    cfig = figure('visible', 'off');
    nm = size(sumdat, 2);
    jsonOut{i} = struct;
    jsonOut{i}.raw = cell(nm, nm);
    tiles = cell(nm, nm);
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
                    cor = corr(X, Y);
                    col = abs(cor)*100;
                    % Desmos graphs for color formula
                    % g\left(x\right)=\min\left(\max\left(x,\ \frac{90}{30}\left(x-30\right)+30,\frac{105}{25}\left(x-70\right)+150\right),255\right)
                    % b\left(x\right)=\max\left(0,\frac{255}{5}\left(x-95\right)\right)
                    col = [255, min(max([col, 3*col-60, 4.2*col-144]), 255), max(0, 51*col-4845)]/255;
                    % set(ax, 'Color', [1 1 0.0667]);
                    d = (1/8:1/4:1)*2*pi;
                    fill(cos(d), sin(d), col)
                    xlim([-sqrt(2) sqrt(2)]/2)
                    ylim([-sqrt(2) sqrt(2)]/2)
                    text(0.5, 0.5, num2str(cor), 'HorizontalAlignment','center','VerticalAlignment','middle')
                end
                % plot(rand(1,(j-1)*nm+file));

                % set(gca, 'units', 'normalized'); %Just making sure it's normalized
                set(ax,'xtick',[], 'ytick', []);
                
                % disp([file-1 nm-j 1 1]/nm)
                % disp((j-1)*nm+file)
                set(ax, 'position', [file-1 nm-j 1 1]/nm);
                tiles{file, j} = ax;
            hold off

            if file == j
                continue;
            end
            f = figure('visible', 'off');
            jsonOut{i}.raw{file, j} = struct;
            hold on
                % set(0, 'CurrentFigure', f);
                scatter(X,Y, 25, linspace(0, 1, size(sumdat,1)), 'filled');
                fname = fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['norm_File_' num2str(j) '.svg']);
                saveas(f, fname);
                jsonOut{i}.raw{file, j}.norm = fileread(fname);
                set(gca, 'YScale', 'log', 'XScale', 'log')
                fname = fullfile(parentF, wantedMods{i}, ['File_' num2str(file)], ['log_File_' num2str(j) '.svg']);
                saveas(f, fname);
                jsonOut{i}.raw{file, j}.log = fileread(fname);
                jsonOut{i}.raw{file, j}.name = [num2str(file) '_' num2str(j)];
            hold off
        end
    end
    fname = fullfile(parentF, wantedMods{i}, 'combined');
    saveas(cfig, [fname '_norm.svg']);
    for ii = 1:nm
        for j = 1:(ii-1)
            set(tiles{j, ii}, 'YScale', 'log', 'XScale', 'log')
        end
    end
    saveas(cfig, [fname '_log.svg']);
    jsonOut{i}.name = wantedMods{i};
    jsonOut{i}.combined = struct;
    jsonOut{i}.combined.norm = fileread([fname '_norm.svg']);
    jsonOut{i}.combined.log = fileread([fname '_log.svg']);
    jsonOut{i}.raw = jsonOut{i}.raw(~cellfun('isempty',jsonOut{i}.raw));
end

fid = fopen(fullfile('Results', params.name, 'Raws', 'linearReg.json'), 'w');
saveJSON(fid, jsonOut);
fclose(fid);