# alphascreen

1. [Features](#features)
2. [Workflow](#workflow)
3. [Installation](#installation)
4. [Usage](#usage)
    * [Job setup](#jobsetup)
    * [Check runs](#checkruns)
    * [Analysis](#analysis)
5. [License](#license)

## Features<a name="features"></a>

* Fetches the fasta sequences of the uniprot IDs contained in a table of paired proteins.

* Chops up the sequences so that they are of reasonable size for subsequent pairwise predictions.

* Interprets the PAEs of the resulting predictions, only considering the region in the plot relating to the interaction.

* Generates summaries, including a PDF showing the PAE plot next to snapshots of the models.

## Workflow<a name="workflow"></a>

### Setting up the fasta files

This package generates fasta files for a set of interaction partners in order to run Alphafold predictions. The assumption is that your Alphafold implementation takes fasta files as input. The input to ```--parse``` is a table which includes two columns containing uniprot IDs for the interaction partners (headers specified with ```--columnA``` and ```--columnB```). These tables can be generated manually or by a database (e.g. BioGRID). Acceptable extensions are .xlsx and .txt (tab-delimited, e.g. as output by BioGRID).

The command ```--parse``` fetches the sequences from Uniprot and fragments them before generating fasta files, which are stored in the *fastas* folder. Fragmenting the sequences (```--fragment```) helps keep the total sequence length short enough so the jobs don't run out of memory. An overlap is considered (```--overlap```) so that the fragmentation doesn't accidentally cut into an interaction interface. You can dimerize any or all proteins (```--dimerize```, ```--dimerize_all```, ```--dimerize_all_except```) and/or consider a specific stretch of sequence for a protein of interest (```--consider```).

```
alphascreen --parse myproteinpairs.xlsx
```

### Running the predictions

The output is a bash script (*runpredictions.bsh*) with a list of commands to run Alphafold on each of the generated fasta files on your machine/cluster. The syntax is set up for the LMB, and will therefore likely not be directly compatible with your implementation. You can sort this by either editing *runpredictions.bsh* or *jobsetup.py* so that the commands will work. If you do change this, make sure the results are output into the *results* directory, which is important for the analysis command ```--show_top``` to work. The package has only been tested on Colabfold 1.3.0 (when parsing the results, it is assumed that the PDBs and PAE json filenames are those from Colabfold).

Important: The script will run all the jobs and relies on your queuing system to handle them!

```
bash runpredictions.bsh
```

### Analyzing the results

The default behaviour for analysis (```--show_top```) is to go through one result at a time and find the best PAE, only considering the region of the plot corresponding to the protein interaction you are screening for. After doing this for all the predictions, they are ranked relative to each other by again comparing their respective PAEs in the protein interaction region.

```
alphascreen --show_top 0.7
```

If you want instead want to rank by iptm score, you can pass ```--rankby iptm```. This relies on your Alphafold/Colabfold implementation writing the iptm and ptm scores to a *scores.txt* file within the individual results directories. It should have the same number of lines as there are models, each line containing the information in this format:

*iptm:0.09 ptm:0.62*

These values are just an example. Note the presence/absence of spaces and placement of colons.

## Installation<a name="installation"></a>

* Set up a fresh conda environment with Python >= 3.8: `conda create -n alphascreen python=3.8`

* Activate the environment: `conda activate alphascreen`.

* Install alphascreen: **`pip install alphascreen`**

* Install pymol dependancies: **`conda install -c schrodinger pymol-bundle`**

## Usage<a name="usage"></a>

### Job setup<a name="jobsetup"></a>

* If you have two large proteins and want to generate fragmented predictions:

```
alphascreen --parse uniprot-id-1/uniprot-id-2 [options]
```

* If you have a list of potential interactors in two columns within a table (.txt or .xlsx):

```
alphascreen --parse filename [options]
```

**Options**

**```--fragment```** *```length```*

Approximate fragment length. Default is 500.

**```--overlap```** *```length```*

Sequence is extended by this amount on either side of slices. Default is 50.

**```--exhaustive```**

Run every protein in column A against every protein in column B, instead of just row by row.

**```--dimerize```** *```uniprot-id```* or *```uniprot-ids.txt```*

Uniprot ID to dimerize. Alternatively, provide a text file (.txt) with a single column list of uniprot IDs to dimerize.

**```--dimerize_all```**

Dimerize all proteins. You may have to reduce the fragment length if the total sequence length becomes too big.

**```--dimerize_all_except```** *```uniprot-ids.txt```*

Provide a text file (.txt) with a single column list of uniprot IDs to NOT dimerize. Everything else will be dimerized.

**```--consider```** *```uniprot/start/end```*

Uniprot ID and sequence range to consider. Example: *Q86VS8/1/200* only considers amino acids 1-200 for uniprot ID Q86VS8.

**```--alphafold_exec```** *```alphafold-executable```*

Path to script that runs Alphafold for writing the commands. Default is *colabfold2* as per the LMB cluster usage.

**```--columnA```** *```columnA-name```*

Name of column heading for uniprot IDs for the first set of interactors. Default is *SWISS-PROT Accessions Interactor A*, which is just what BioGRID uses.

**```--columnB```** *```columnB-name```*

Name of column heading for uniprot IDs for the second set of interactors. Default is *SWISS-PROT Accessions Interactor B*, which is just what BioGRID uses.

**```--focus```** *```uniprot-id```*

Uniprot ID to focus on. This means that it will the first chain in any predictions that contain it.

**```--dontwrite```**

Don't write out any files and just show the relevant information. This is useful to check how many fastas will be generated from parsing and fragmenting.

### Check runs<a name="checkruns"></a>

* Check how many runs are finished so far and how many remain.

```
alphascreen --check
```

* Check how many runs are finished so far and write out a new bash script with the remaining Alphafold commands (runpredictions-unfinished.bsh). This is useful when the jobs crash.

```
alphascreen --write_unfinished
```

Be careful not to run this while there are predictions in progress.

### Analysis<a name="analysis"></a>

```
alphascreen --show_top threshold [options]
```

Generate summary files for the runs so far. For example, ```alphascreen --show_top 0.7``` will rank predictions by interaction-site PAEs to choose the highest rank, then lists those predictions, ranking by the interaction-site PAE. Only those with scaled PAEs higher than 0.7 are shown. See the ```--rankby``` option below for more information on the scaled PAE. To output all predictions, pass ```--show_top 0```. A table is output (.xlsx and .csv), and the .xlsx can be used as input for a subsequent run of alphascreen if you need to test dimerization or use different alphafold executable on the top hits.

```
alphascreen --write_table [options]
```

Like ```--show_top```, but only outputs the table (.xlsx and .csv). No threshold value is considered; all predictions are ranked and output.

**Options**

**```--rankby```** ```pae``` or ```iptm``` or ```ptm```

Score by which models are ranked (***pae***, ***iptm***, or ***ptm***). The default is *pae*. This is used for both choosing the best model in a prediction as well as ranking those chosen models in the summary files. The option ```pae``` will look for the deepest PAE valleys only in the parts of the plot that are interactions between **different** proteins. The PAE is scaled to be between 0 and 1 where higher values are better predictions (```--show_top 0.7``` is a good starting point). The options ```iptm``` and ```ptm``` rely on a *scores.txt* file in each results directory (see explanation at the top) (in this case ```--show_top 0.3 --rankby iptm``` is a good starting point).

**```--overwrite```**

Overwrite snapshots that have already been generated, otherwise it will skip those to save time. This is only relevant for ``show_top``.

## License<a name="license"></a>

This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/sami-chaaban/alphascreen/blob/main/LICENSE.txt) file for details.