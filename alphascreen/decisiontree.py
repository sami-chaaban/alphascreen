from unipressed import UniprotkbClient
import pandas as pd
import os
from pathlib import Path
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from textwrap import wrap
import cv2 as cv
import sys

import alphascreen
from alphascreen import argparser
from alphascreen import jobsetup
from alphascreen import analyze

def decide():

    ##################################
    """
    This is the master function that is used to determine what the user is asking for from the command-line arguments,
    and set up the variables required to call those functions
    """
    ##################################

    #Get the passed parameters
    params = argparser.argparse()

    table = params['table']
    focus = params['focus']
    fraglen = params['fraglen']
    overlap = params['overlap']
    step = params['step']
    dimerize = params['dimerize']
    dimerize_all = params['dimerize_all']
    dimerize_except = params['dimerize_except']
    consider = params['consider']
    dontwrite = params['dontwrite']
    alphafold_exec = params['alphafold_exec']
    columnA = params['columnA']
    columnB = params['columnB']
    check = params['check']
    check_write = params['check_write']
    threshold = params['threshold']
    overwrite = params['overwrite']
    writetable = params['writetable']
    rankby = params['rank_by']
    exhaustive = params['exhaustive']
    showall = params['showall']
    ignoreself = params['ignoreself']
    customids = params['customids']

    #--overlap and --step are two different ways of relating the fragment length to the spacing
    #between fragments, so only one of them can be in play at a time
    if step != "" and overlap is not None:
        sys.exit("\n>> Error: --overlap and --step cannot be combined. Use --overlap to extend abutting slices on either side, or --step to slide a window of --fragment amino acids along the sequence.\n")

    if overlap is None:
        overlap = 50

    if fraglen == "":
        fraguniprot = 0
        fraglen = 500
    elif len(fraglen.split("/")) == 2:
        if len(table.split("/")) != 2:
            sys.exit("\n>> Error: You provided two fragment lengths but the input was not two uniprot IDs.\n")
        fraguniprot = [table.split("/")[0], table.split("/")[1]]
        fraglen = [int(fraglen.split("/")[0]), int(fraglen.split("/")[1])]
    else:
        fraguniprot = 0
        fraglen = int(fraglen)

    if step == "":
        stepuniprot = 0
        step = None
    elif len(step.split("/")) == 2:
        if len(table.split("/")) != 2:
            sys.exit("\n>> Error: You provided two step lengths but the input was not two uniprot IDs.\n")
        stepuniprot = [table.split("/")[0], table.split("/")[1]]
        step = [int(step.split("/")[0]), int(step.split("/")[1])]
    else:
        stepuniprot = 0
        step = int(step)

    if step is not None:
        for s in (step if type(step) == list else [step]):
            if s < 0:
                sys.exit("\n>> Error: the step length cannot be negative.\n")

        #A step larger than the fragment length skips over residues entirely, which is almost never intended
        if table != "" and type(step) != list and type(fraglen) != list and step > fraglen:
            print("\n!! Warning: the step (" + str(step) + ") is larger than the fragment length (" + str(fraglen) + "), so stretches of sequence between fragments will not be covered.\n")

    #An overlap as large as the fragment length itself makes each slice at least three times the
    #spacing between slices, and beyond that the slices start running off the ends of the sequence
    elif table != "" and overlap >= min(fraglen if type(fraglen) == list else [fraglen]):
        print("\n!! Warning: the overlap (" + str(overlap) + ") is large relative to the fragment length (" + str(fraglen) + "), so the slices will be heavily redundant and clipped at the ends of the sequence. Use --step instead to slide a window of a fixed length along the sequence.\n")

    #log.txt is about to be overwritten with the current version, so the version that set this
    #directory up has to be read before that happens
    if table != "" and os.path.exists('results') and len(os.listdir('results')) > 0:
        jobsetup.warnoldversion(reparsing=True)

    if table != "" and not dontwrite:

        with open("log.txt", 'w') as f:
            f.write("AlphaScreen version: " + str(alphascreen.__version__) + "\n")
            f.write("Input: " + table + "\n")
            f.write("Fragment: " + str(fraglen) + "\n")
            if step is not None:
                f.write("Step: " + str(step) + "\n")
            else:
                f.write("Overlap: " + str(overlap) + "\n")
            if dimerize:
                f.write("Dimerized: " + dimerize + "\n")
            elif dimerize_all:
                f.write("Dimerized all proteins\n")
            elif dimerize_except:
                f.write("Dimerize all except: " + dimerize_except + "\n")
            if consider != "":
                f.write("Consider: " + consider + "\n")
            if customids != "":
                f.write("Custom IDs: " + customids + "\n")
            f.write("Alphafold executable: " + alphafold_exec + "\n")
            if exhaustive:
                f.write("Exhaustive = True.\n")
            if ignoreself:
                f.write("Ignore self = True.\n")


    ##################################
    #Parse input

    towrite = not dontwrite

    #These all read the fragment names that are already on disk rather than generating new ones
    if check or check_write or writetable or showall or threshold != -1:
        jobsetup.warnoldversion()

    if check:
        towrite = False
        jobsetup.findunfinished(alphafold_exec, write=towrite)
        sys.exit()
    elif check_write:
        towrite = True
        jobsetup.findunfinished(alphafold_exec, write=towrite)
        sys.exit()

    if rankby=="pae":
        rankby="scaledPAE"

    if showall:
        threshold = 0

    ##################################

    if writetable:
        print("\n>> Parsing results...")
        df = analyze.getscores(rankby)
        analyze.write_top(df, 0, rankby)
        sys.exit()

    ##################################

    if threshold != -1:
        print("\n>> Parsing results...")
        df = analyze.getscores(rankby)
        
        if df.empty:
            sys.exit("\n>> Error: no results could be found.\n")
        try:
            df[df[rankby]>threshold]
        except TypeError:
            sys.exit("\n>> Error: some values you want to rank by don't seem to exist. Did you try to rank by iptm without those scores existing?\n")

        if df[df[rankby]>threshold].empty:
            sys.exit("\n>> Error: no results could be found with "+rankby+" above " + str(threshold) + ".\n")

        analyze.summarize_pae_pdf(df, threshold, rankby)
        analyze.write_top(df, threshold, rankby)
        analyze.write_modelpngs(df, threshold, rankby, overwrite=overwrite)
        analyze.summarize_paeandmodel_pdf(df, threshold, rankby)
        sys.exit()

    ##################################

    #The start is converted to a zero-based index and the end is left as it is, so that the two
    #are used directly as the ends of a python slice and the range given is inclusive of both
    consideruniprot, considerstart, considerend = [], [], []
    if consider != "":
        if consider[-4:] == ".txt":
            with open(consider) as f:
                considerlist=[line.strip() for line in f.readlines()]
            for c in considerlist:
                considerargs = c.split("/")
                if len(considerargs) != 3:
                    sys.exit("\n>> Error: there is a problem with the line " + c + " in " + consider + ".\n")
                consideruniprot.append(considerargs[0])
                considerstart.append(int(considerargs[1])-1)
                considerend.append(int(considerargs[2]))
        else:
            considerargs = consider.split("/")
            if len(considerargs) != 3:
                sys.exit("\n>> Error: the --consider argument was not passed properly.\n")
            consideruniprot = [considerargs[0]]
            considerstart = [int(considerargs[1])-1]
            considerend = [int(considerargs[2])]

    for c, cs, ce in zip(consideruniprot, considerstart, considerend):
        if cs < 0:
            sys.exit("\n>> Error: the first amino acid to consider for " + c + " is numbered from 1, so it cannot be lower than that.\n")
        if ce <= cs:
            sys.exit("\n>> Error: the last amino acid to consider for " + c + " is not after the first one.\n")

    filetype=""
    if table != "":
        if table[-4:] == ".txt":
            filetype = "table"
        elif table[-4:] == "xlsx":
            filetype = "excel"

        if filetype!="":
            Ainteractors, Binteractors = jobsetup.getinteractors(table, filetype, columnA, columnB, focus, exhaustive)
        else: #if two proteins
            Ainteractors = [table.split("/")[0]]
            Binteractors = [table.split("/")[1]]
            print("\n>> Parsing " + Ainteractors[0] + " and " + Binteractors[0] + "...\n")

        jobsetup.getfastas_writecommands(Ainteractors, Binteractors, consideruniprot, considerstart, considerend, customids, split=True,fraguniprot=fraguniprot,fraglen=fraglen,overlap=overlap,stepuniprot=stepuniprot,step=step,dimerize=dimerize,dimerize_all=dimerize_all,dimerize_except=dimerize_except,write=towrite,alphafold_exec=alphafold_exec, ignoreself=ignoreself)

    ##################################

