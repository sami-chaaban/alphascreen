# alphascreen

Use this package to generate fastas for a set of interaction partners to run Alphafold predictions.

## Installation<a name="installation"></a>

* Set up a fresh conda environment with Python >= 3.8: `conda create -n alphascreen python=3.8`

* Activate the environment: `conda activate alphascreen`.

* Install alphascreen: **`pip install alphascreen`**

* Install pymol dependancies: **`conda install -c schrodinger pymol-bundle`**

## Usage<a name="installation"></a>

### Job setup

```
alphascreen --parse filename [options]
```

**Options**

**```--fragment```** *```length```*

Approximate fragment length. Default is 500.

**```--overlap```** *```length```*

Sequence is extended by this amount on either side of slices. Default is 50.

**```--dimerize```** *```uniprot-id```*

Dimerize this uniprot ID whenever it is found.

**```--consider```** *```uniprot/start/end```*

Uniprot ID and sequence range to consider. Example: "Q86VS8/1/200" only considers amino acids 1-200 for uniprot ID Q86VS8.

**```--alphafold_exec```** *```path-to-colabfold-executable```*

Path to script that runs Colabfold for writing the commands. Default is "colabfold2".

**```--columnA```** *```columnA-name```*

Name of column heading for uniprot IDs for the first set of interactors. Default is "SWISS-PROT Accessions Interactor A".

**```--columnB```** *```columnB-name```*

Name of column heading for uniprot IDs for the second set of interactors. Default is "SWISS-PROT Accessions Interactor B".

### Check runs

```
alphascreen --check
```

Checks how many runs are finished so far.

```
alphascreen --write_unfinished
```

Checks how many runs are finished so far and writes out a new bash script with the remaining Colabfold commands.

### Analyze results

```
alphascreen --show_top threshold (--overwrite)
```

Generate summary files for the runs so far. Only iptms above the threshold value provided will be considered (e.g. ```alphascreen --show_top 0.3``` for iptms above 0.3). Pass the "overwrite" option to overwrite snapshots that have already been generated, otherwise it will skip those to save time.

```
alphascreen --write_table
```

Output all the results into a table ranked by iptm score.