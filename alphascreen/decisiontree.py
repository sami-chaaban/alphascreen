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
from pymol import cmd
import cv2 as cv
import sys

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
    rankby = params['rankby']
    exhaustive = params['exhaustive']

    if table != "" and not dontwrite:

        with open("log.txt", 'w') as f:
            f.write("Input: " + table + "\n")
            f.write("Fragment: " + str(fraglen) + "\n")
            f.write("Overlap: " + str(overlap) + "\n")
            if dimerize:
                f.write("Dimerized: " + dimerize + "\n")
            elif dimerize_all:
                f.write("Dimerized all proteins\n")
            elif dimerize_except:
                f.write("Dimerize all except: " + dimerize_except + "\n")
            if consider != "":
                f.write("Consider: " + consider + "\n")
            f.write("Alphafold executable: " + alphafold_exec + "\n")
            if exhaustive:
                f.write("Exhaustive = True.\n")


    ##################################
    #Parse input

    towrite = not dontwrite

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

    consideruniprot, considerstart, considerend = ["", "", ""]
    if consider != "":
        considerargs = consider.split("/")
        if len(considerargs) != 3:
            sys.exit("\n>> Error: the --consider argument was not passed properly.\n")
        consideruniprot = considerargs[0]
        considerstart = int(considerargs[1])-1
        considerend = int(considerargs[2])-1

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

        jobsetup.getfastas_writecommands(Ainteractors, Binteractors, consideruniprot, considerstart, considerend, split=True,fraglen=fraglen,overlap=overlap,dimerize=dimerize,dimerize_all=dimerize_all,dimerize_except=dimerize_except,write=towrite,alphafold_exec=alphafold_exec)

    ##################################

