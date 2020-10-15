# Data Visualizer

A data visualizer for Mass Spectrometry data.

## Running

### params.p

Parameter file created within the home directory. 

#### Format
```py
[GLOBAL]
# Global parameters here

name = Dataset

[filename.fileextension]
# File specific parameters here

numericobject = 1.0

stringobject = Hello World!

listobject = [
stringvalue1
stringvalue2
]
```

### [combinedHM.m](combinedHM.m)

Run using `combinedHM` within matlab

#### Function

Generates Mod Mapper and Heatmap data. Heatmap data necessary to run most later programs.

#### Parameters

 * `mods`: list of modifications for the heatmap to generate filtered versions for. `All` refers to no filter.
 * `proteins`: list of protein names to use for ModMapper.

### [linearReg.m](linearReg.m)

Run using `linearReg` within matlab

#### Parameters

 * `name`: string referring to name of the wanted dataset.

#### Function

Linear regression analysis between data files comparing protein spectral counts.

### [vennclustermap.r](vennclustermap.r)

Run using `rscript vennclustermap.r`

#### Parameters

 * `name`: string referring to name of the wanted dataset.
 * `vennImgSize`: resolution of the venn diagram image in pixels. [default 2800]
 * `heatmapcount`: total number of rows to select from heatmap [default 64]
 * `heatmapcolors`: color palette used to generate gradient for heatmap (low to high).

#### Function

Venn diagrams and set overlap calculations between data sets, as well as cluster heatmaps filtered and sorted by various functions including rank, standard deviation, and the coefficient of variation.

## Developing
 * `source build.sh` to set up virtual enviroment and matlab CLI alias (default uses `MATLAB R2019b`)

## Changelog

* 2020/10/14: Add `clear colors` button for column coloring in web viewer and implement remaining statistics tests.
* 2020/10/13: Add stat tests applied to test groups (`statTests.r`).
* 2020/10/12: Fix JSON saving and loading bugs, allow nested lists in `params.p`.
* 2020/10/06: Finished custom `params.p` parsing for R and MATLAB
* 2020/08/26: Added custom `params.p` for customization
* 2020/07/28: Added combined plot without colour
* 2020/07/27: Fixed missing sheets in ModMapper xlsx raws