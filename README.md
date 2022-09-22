# alphascreen

Use this package to generate fastas for a set of interaction partners to run Alphafold predictions. The input is a table which includes two columns containing uniprot IDs for the interaction partners. The sequences are fetched from Uniprot and fragmented before generating fasta files. Fragmenting helps keep the total sequence length short enough so the jobs don't run out of memory.

The output is a bash script that allows you to run Alphafold on all the generated fasta files. The syntax for job submission will likely not correspond to what you use in your system. You can either edit the output file "colabfoldrun.bsh" or "jobsetup.py" itself so that the alphafold submission commands have the right syntax. If you do change this, make sure the results are output into a "results" directory, which is important for the analysis command "show_top" to work. The package has only been tested on Colabfold 1.3.0 and therefore its file naming system.

The results can be ranked by iptm score and compiled into a pdf showing all PAEs next to snapshots of the predictions.

## Installation<a name="installation"></a>

* Set up a fresh conda environment with Python >= 3.8: `conda create -n alphascreen python=3.8`

* Activate the environment: `conda activate alphascreen`.

* Install alphascreen: **`pip install alphascreen`**

* Install pymol dependancies: **`conda install -c schrodinger pymol-bundle`**

## Bugs<a name="bugs"></a>

* The pymol package often crashes or hangs during model-snapshot generation.

## Usage<a name="usage"></a>

### Job setup

```
alphascreen --parse filename [options]
```

**Options**

**```--fragment```** *```length```*

Approximate fragment length. Default is 500.

**```--overlap```** *```length```*

Sequence is extended by this amount on either side of slices. Default is 50.

**```--dimerize```** *```uniprot-id```* *or* *```uniprot-ids.txt```*

Dimerize this uniprot ID whenever it is found.

**```--dimerize_all```**

Dimerize all proteins.

**```--dimerize_all_except```** *```uniprot-ids.txt```*

Provide a text file (.txt) with a single column list of uniprot IDs to NOT dimerize. Everything else will be dimerized.

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

Generate summary files for the runs so far. Only iptms above the threshold value provided will be considered (e.g. ```alphascreen --show_top 0.3``` for iptms above 0.3). Optionally, pass the "--overwrite" option to overwrite snapshots that have already been generated, otherwise it will skip those to save time.

```
alphascreen --write_table
```

Output all the results into a table ranked by iptm score.