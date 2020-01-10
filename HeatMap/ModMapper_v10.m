clear all

TempFile = uigetfile('.xlsx','Choose Datafile');
dat = readtable(TempFile, 'Sheet', 'Spectra'); %spectra tab
summarydat = readtable(TempFile, 'Sheet', 'Summary'); %summary data
proteindat = readtable(TempFile, 'Sheet', 'Proteins'); %summary data

TempFastaFile = uigetfile('.fasta','Choose Datafile');
fastaFile = struct2table(fastaread(TempFastaFile));

wantedGenes = splitlines(fileread('wantedProteins.txt'));

for i = 1:length(wantedGenes)
    proteinName = getProteinName(wantedGenes{i}, proteindat.Description, 3);
    if isempty(proteinName)
        continue;
    end
    proteinName = proteinName{1};
    getPeptideMap(proteinName, dat, summarydat, fastaFile, ['modformat_' wantedGenes{i}]);
end
% proteinName = '>sp|P14873|MAP1B_MOUSE Microtubule-associated protein 1B OS=Mus musculus OX=10090 GN=Map1b PE=1 SV=2'; %desired rank for protein to analyze