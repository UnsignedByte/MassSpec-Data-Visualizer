clear all
tic;

TempFiles = uigetfile('.xlsx','Choose Datafile', 'Multiselect', 'on');
TempFastaFile = uigetfile('.fasta','Choose Datafile');
fastaFile = struct2table(fastaread(TempFastaFile));


wantedGenes = splitlines(fileread('wantedProteins.txt'));

for kk = 1:numel(TempFiles)
    TempFile = TempFiles{kk};
    
    resfolder = fullfile('Results', getResultFolder(TempFile), 'ModMapper');

    if ~isfolder(resfolder)
        mkdir(resfolder);
    end

    dat = readtable(TempFile, 'Sheet', 'Spectra'); %spectra tab
    summarydat = readtable(TempFile, 'Sheet', 'Summary'); %summary data
    proteindat = readtable(TempFile, 'Sheet', 'Proteins'); %summary data
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