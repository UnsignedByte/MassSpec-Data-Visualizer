# Data Visualizer

Detailed instructions on how to interpret output data.

## FileIds

Table containing File ID, name, and test group information. `Test_Group` column can be manually altered to change test groups.

| ID | Filename  | Test\_Group |
|----|-----------|-------------|
| 1  | FileName1 | 1           |
| 2  | FileName2 | 1           |
| 3  | FileName3 | 2           |
| 4  | FileName4 | 2           |

## StatTests

Results from running statistics tests comparing pairs of test groups.

### Files

Column names are an underscore-separated list of Test Group IDs (taken from FileIds.csv). Rows correspond to the same row from HeatMap. Each cell represents the resulting p-value from the stat test.

Data separated by Modification Name (I.E. All) and  Test Type (I.E. wilcoxon).


### Details

Filename formatted using `<modification-filter>_<test-name>.csv`.

## LinearReg

Linear regression data on # of spectral counts.

### Files

Data separated by Modification Name (I.E. All), compared files, and axes type (I.E. log). Also includes a combined svg image with all regression plots and pearson's r correlation coefficients.

### Details

Filepath formatted using `<modification-filter>/File_<File1_ID>/<axes_scale>_File_<File2_ID>.svg`.

## ModMapper

Peptide-specific modification information.

### Files

Data separated by Gene Name.

Each file contains a **Summary** Sheet with summary information. This contains information on which protein sequence is being used in each subsequent sheet, as well as the files included in that sheet.

| SheetNumber | Filenames | Protein\_Name |
|-------------|-----------|--------------|
| 1           | FileName1 | ProteinName1 |
| 2           | FileName2 | ProteinName2 |

Subsequent sheets will contain information on the existence of different modifications at each amino acid and on each file.

| Amino Acid | ModA\_File1 | ModB\_File1 | Total\_File1 | ModA\_File2 | ModB\_File2 | Total\_File2 | ... |
|------------|------------|------------|-------------|------------|------------|-------------|-----|
| A          | 0          | 1          | 1           | 0          | 0          | 0           | ... |
| D          | 1          | 2          | 3           | 2          | 0          | 2           | ... |

Each cell represents the number of the observed modification found at the amino acid.

### Details

Filename formatted using `<Gene\_Name>.xlsx`.

Files must occasionally be separated into different sheets as isotopes of proteins may appear in each file, resulting in a slightly different protein sequence that cannot be directly compared.

## HeatMap

HeatMap containing spectral data and basic analysis.

### Files

Data separated by Modification Name.

Each file contains a sheet with HeatMap data. Data is separated into groups of proteins in each rank, with `class` rows containing grouped information.

| Rank\_Number | Protein\_Name | Gene\_Name | x\_OfSpectra\_1 | x\_OfSpectra\_2 | max\_x\_OfSpectra | x\_AA\_sInProtein | Contaminant | Row\_Type |
|--------------|---------------|------------|-----------------|-----------------|-------------------|-------------------|-------------|-----------|
| 1            | A/B           | C/D        | 100             | 51              | 100               | 987               | 0           | 1         |
| 1            | A             | C          | 100             | 1               | 100               | 975               | 0           | 0         |
| 1            | B             | D          | NaN             | 50              | 50                | 987               | 0           | 0         |

Each cell represents the number of the observed modification found at the amino acid.

### Details

Filename formatted using `<Modification_Name>.csv`.

Protein Names for classes are `/` separated lists of the parsed protein name taken from `>sp|<Protein_Name>|rest`.
Gene Names are similarly `/` separated in classes, and each gene name is taken from the `GN=` tab.

## VennDiagram

Data separated by Modification Name.

Each file contains an svg with the venn diagram, as well as two raw files. Column names of raw files represent the test groups included in the column (X1011 means test groups 1, 3, and 4 are included). Each column contains a list of IDs or Gene Names of all the proteins contained within the test group. Total counts are displayed within the Venn Diagram.

### Details

A protein is considered *included* in a test group if its spectral count is greater than a threshold proportion of the max spectral count across all groups (default `vennCutoff=0.8`).

## ClusterHeatMap

Data separated by Modification Name and sorting function.

Each file contains a svg with two heatmaps, one with files and another with averaged testgroup data. Heatmaps rows will be labeled with Protein Ranks or Gene names.

### Details

Custom functions can be defined in `vennclustermap.r`, and heatmap row count can be determined in `params.p`; default `hmcount=64`.