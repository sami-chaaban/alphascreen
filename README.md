# alphascreen

1. [Features](#features)
2. [Workflow Overview](#workflow)
3. [Installation](#installation)
4. [Usage](#usage)
    * [Job setup](#jobsetup)
    * [Check runs](#checkruns)
    * [Analysis](#analysis)
5. [Troubleshooting](#troubleshooting)
6. [License](#license)

## Features<a name="features"></a>

* Fetches the fasta sequences of the uniprot IDs contained in a table of paired proteins.

* Chops up the sequences so that they are of reasonable size for subsequent pairwise predictions.

* Interprets the PAEs of the resulting predictions, only considering the region in the plot relating to the interaction.

* Generates summaries, including a PDF showing the PAE plot next to snapshots of the models.

## Workflow Overview<a name="workflow"></a>

For a step-by-step workflow, see the [Usage](#usage) section.

### Setting up the fasta files

This package generates fasta files for a set of interaction partners in order to run Alphafold predictions. The assumption is that your Alphafold implementation takes fasta files as input. The input to ```--parse``` is a table which includes two columns containing uniprot IDs for the interaction partners (column headers specified with ```--columnA``` and ```--columnB```). These tables can be generated manually or by a database (e.g. BioGRID). Acceptable extensions are .xlsx and .txt (tab-delimited, e.g. as output by BioGRID).

The command ```--parse``` fetches the sequences from Uniprot and fragments them before generating fasta files, which are stored in the *fastas* folder. Fragmenting the sequences (```--fragment```) helps keep the total sequence length short enough so the jobs don't run out of memory. An overlap is considered (```--overlap```) so that the fragmentation doesn't accidentally cut into an interaction interface. You can dimerize any or all proteins (```--dimerize```, ```--dimerize_all```, ```--dimerize_all_except```) and/or consider a specific stretch of sequence for a protein or set of proteins (```--consider```).

```
alphascreen --parse myproteinpairs.xlsx
```

### Running the predictions

One of the files that is output is a bash script (*runpredictions.bsh*) with a list of commands to run Alphafold on each of the generated fasta files on your machine/cluster. The syntax is set up for the LMB, and will therefore likely not be directly compatible with your implementation. You can sort this by either editing *runpredictions.bsh* or *jobsetup.py* so that the commands will work. If you do change this, make sure the results are output into the *results* directory, which is important for the analysis command ```--show_top``` to work. The package has only been tested on Colabfold 1.3.0, 1.4.0, and 1.5.2 (when parsing the results, it is assumed that the PDBs and PAE json filenames are those from Colabfold).

Important: The script will run all the jobs and relies on your queuing system to handle them!

```
bash runpredictions.bsh
```

### Analyzing the results

The default behaviour for analysis (```--show_top```) is to go through one result at a time and find the best PAE, only considering the region of the plot corresponding to the protein interaction you are screening for. After doing this for all the predictions, they are ranked relative to each other by again comparing their respective PAEs in the protein interaction region.

```
alphascreen --show_top 0.7
```

If you want instead want to rank by local interaction score (LIS) (Kim et al., 2024), you can pass ```--rank_by lis```. While iptm and ptm can be used for ranking as well, alphascreen relies on your Alphafold/Colabfold implementation writing the iptm and ptm scores to a *scores.txt* file within the individual results directories. It should have the same number of lines as there are models, each line containing the information in this format:

*iptm:0.09 ptm:0.62*
*iptm:0.02 ptm:0.63*
*...*

These values are just an example. Note the presence/absence of spaces and placement of colons.

## Installation<a name="installation"></a>

* Set up a fresh conda environment with Python >= 3.8: `conda create -n alphascreen python=3.8`

* Activate the environment: `conda activate alphascreen`.

* Install alphascreen: **`pip install alphascreen`**

* Install pymol dependancies: **`conda install -c schrodinger pymol-bundle`**

## Usage<a name="usage"></a>

### Job setup<a name="jobsetup"></a>

If you have a list of potential interactors in two columns within a table (.txt or .xlsx):

```
alphascreen --parse filename [options]
```

For example, the tab delimited Biogrid output *BIOGRID-GENE-111808-4.4.214.tab3.txt* looks like this:

![BiogridExample](https://github.com/sami-chaaban/alphascreen/blob/main/examples/Biogrid.png?raw=true "BiogridExample")

To parse the file, run `alphascreen --parse BIOGRID-GENE-111808-4.4.214.tab3.txt`. The script will generate the following output:

```
>> Reading table...

>> Removed 45 duplicated pairs...

>> Getting sequences...

!! Warning: there was a problem getting the Uniprot data for accessions:

·· -

>> Writing 157 Alphafold commands...

>> Run the Alphafold jobs with "bash runpredictions.bsh"

Done!

```

Sometimes the uniprot ID is invalid, such as the dash ("-") referenced in the warning above. When finished, the script will have generated two folders (***fastas*** and ***results***), some general information files (***log.txt*** and ***uniprots.txt***), and the bash script to run the predictions (***runpredictions.bsh***). See below for additional options, such as dimerization, custom sequences, etc.

If you only have two proteins and want to generate fragmented predictions, you can simply run:

```
alphascreen --parse uniprot-id-1/uniprot-id-2 [options]
```


**Options**

**```--columnA```** *```columnA-name```*

Name of column heading for uniprot IDs for the first set of interactors. Default is *SWISS-PROT Accessions Interactor A*, which is what BioGRID uses.

**```--columnB```** *```columnB-name```*

Name of column heading for uniprot IDs for the second set of interactors. Default is *SWISS-PROT Accessions Interactor B*, which is what BioGRID uses.

**```--fragment```** *```length```*

Approximate fragment length. Default is 500. If you provided two uniprot IDs as input, you have the option to provide two fragment lengths here (*length1/length2*) corresponding to the fragment length of each uniprot ID.

**```--overlap```** *```length```*

Sequence is extended by this amount on either side of slices. Default is 50.

**```--exhaustive```**

Run every protein in column A against every protein in column B, instead of just row by row. This is useful when you want to run everything in a list of proteins against each other, or one protein against many (see example of the latter below). In this case, create a table with two columns both containing the list of proteins, and use the *--exhaustive* option. Any duplicate pairs will be removed during parsing. To ignore proteins being predicted against themselves, use *--ignore_self*.

<img src="https://github.com/sami-chaaban/alphascreen/blob/main/examples/Exhaustive.png"  alt="ExhaustiveExample" width="60%" height="60%">

**```--ignore_self```**

Ignore any instances of a protein being predicted against itself. Useful when using the *--exhaustive* option with two identical columns of proteins to predict everything against everything while ignoring self. This does not include dimerization.

**```--dimerize```** *```uniprot-id```* or *```uniprot-ids.txt```*

Uniprot ID to dimerize (i.e. homodimerize). Alternatively, provide a text file (.txt) with a single column list of uniprot IDs to dimerize.

**```--dimerize_all```**

Dimerize all proteins. You may have to reduce the fragment length if the total sequence length becomes too big.

**```--dimerize_all_except```** *```uniprot-id```* or *```uniprot-ids.txt```*

Dimerize all proteins except the uniprot ID provided here. Alternatively, provide a text file (.txt) with a single column list of uniprot IDs to NOT dimerize. Everything else will be dimerized.

**```--consider```** *```uniprot/start/end```* or *```consider.txt```*

Uniprot ID and sequence range to consider. Example: *Q86VS8/1/200* only considers amino acids 1-200 for uniprot ID Q86VS8. Alternatively, provide a text file with a single column list of sequences to consider with a similar syntax (*id/start/end*).

**```--customids```** *```sequences.fasta```*

Provide a fasta file to define custom IDs that you will use when running *--parse*. This is useful when your sequence is not present in uniprot database. The fasta file can contain multiple entries. Make sure the headers don't contain any dashes and have at most one underscore (working examples: *>ProteinA* or *>ProteinA_Tetrahymena*).

**```--alphafold_exec```** *```alphafold-executable```*

Path to script that runs Alphafold for writing the commands. Default is *colabfold4* as per the LMB cluster usage.

**```--focus```** *```uniprot-id```*

Uniprot ID to focus on. This means that it will the first chain in any predictions that contain it.

**```--dontwrite```**

Don't write out any files and just show the relevant information. This is useful to check how many fastas will be generated from parsing and fragmenting.

### Check runs<a name="checkruns"></a>

```
alphascreen --check
```

Check how many runs are finished so far and how many remain.

```
alphascreen --write_unfinished
```

Check how many runs are finished so far and write out a new bash script with the remaining Alphafold commands (runpredictions-unfinished.bsh). This is useful when the jobs crash.

Be careful not to run this while there are predictions in progress.

### Analysis<a name="analysis"></a>

```
alphascreen --show_top threshold [options]
```

Generate summary files for the runs so far. For example, ```alphascreen --show_top 0.7``` will choose the model with the best interaction-site PAE for each prediction, and then rank all the predictions (also by the interaction-site PAE). Only those with scaled-PAEs higher than the threshold (e.g. 0.7) are output to a pdf. The ```--rank_by``` option below has more information on the scaled-PAE and LIS. To output all predictions, pass ```--show_all```. The PAEs-Models pdf will also contain a snapshot of the prediction (see troubleshooting below if you get a "license" watermark on the models). A table is also output (.xlsx) that can be used as input for a subsequent run of alphascreen if you need to test different parameters on just the top hits.

Here is an example of one page from the PAEs-Models pdf:

![AnalysisExample](https://github.com/sami-chaaban/alphascreen/blob/main/examples/Analysis.png?raw=true "AnalysisExample")

```
alphascreen --write_table [options]
```

Like ```--show_top```, but only outputs the table (.xlsx and .csv). No threshold value is considered; all predictions are ranked and output.

**Options**

**```--rank_by```** ```pae``` or ```lis``` or ```iptm``` or ```ptm```

Score by which models are ranked (***pae***, ***lis***, ***iptm***, or ***ptm***). The default is *pae*. This is used for both choosing the best model in a prediction as well as ranking those chosen models in the summary files. The option ```pae``` will look for the deepest PAE valleys only in the parts of the plot that are interactions between **different** proteins. The PAE is scaled to be between 0 and 1 where higher values are better predictions (```--show_top 0.7``` is a good starting point).

The local interaction score (option ```lis```) (Kim et al., 2024) similarly analyzes those regions of the PAE (see the manuscript for more information). Great LIS scores are above 0.1 but the values can get quite small in some cases (```--rank_by lis --show_all``` is probably best).

The options ```iptm``` and ```ptm``` rely on a *scores.txt* file in each results directory (see explanation at the top) (in this case ```--rank_by iptm --show_top 0.3``` is a good starting point).

**```--overwrite```**

Overwrite .png snapshots that have already been generated, otherwise it will skip those to save time. This is only relevant for ``show_top``.

## Troubleshooting<a name="troubleshooting"></a>

* If a "license expired" watermark appears on the model images in PAEs-Models.pdf, follow the instructions here to update your pymol license: https://pymol.org/dsc/ip/license/. You may have to be on your institute's VPN for the link to direct you to a download page rather than an invoice verification page.

## License<a name="license"></a>

This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/sami-chaaban/alphascreen/blob/main/LICENSE.txt) file for details.