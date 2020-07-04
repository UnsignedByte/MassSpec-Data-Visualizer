# Data Visualizer

A data visualizer for Mass Spectrometry data.

---

## Running

### `Params/`

Contains parameter files for the various files.

### `combinedHM.m`

#### Function

Generates Mod Mapper and Heatmap data. Heatmap data necessary to run most later programs.

#### Parameters

 * `mods.txt` contains a list of newline-separated modifications for the heatmap to generate filtered versions for. `All` refers to no filter.
 * `proteins.txt` contains a list of protein names to use for ModMapper.

### `linearReg.m`

#### Function

Linear regression analysis between data files comparing protein spectral counts.

### `vennclustermap.r`

#### Function

Venn diagrams and set overlap calculations between data sets, as well as cluster heatmaps filtered and sorted by various functions including rank, standard deviation, and the coefficient of variation.

---

## Developing
 * `source build.sh` to set up virtual enviroment and matlab CLI alias (default uses `MATLAB R2019b`)
