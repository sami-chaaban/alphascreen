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
    fraglen = params['fraglen']
    overlap = params['overlap']
    dimerize = params['dimerize']
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

    ##################################

    if writetable:
        print("\n>> Parsing results...")
        df = analyze.getscores()
        analyze.write_top(df, 0)
        sys.exit()

    ##################################

    if threshold != -1:
        print("\n>> Parsing results...")
        df = analyze.getscores()
        if df.empty:
            sys.exit("\n>> Error: no results could be found.\n")
        elif df[df['iptm']>threshold].empty:
            sys.exit("\n>> Error: no results could be found with an iptm above " + str(threshold) + ".\n")
        analyze.summarize_pae_pdf(df, threshold)
        analyze.write_top(df, threshold)
        analyze.write_modelpngs(df, threshold, overwrite=overwrite)
        analyze.summarize_paeandmodel_pdf(df, threshold)
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
        if table[-4:] == "txt":
            filetype = "table"
        elif table[-4:] == "xlsx":
            filetype = "excel"
        else:
            sys.exit("\n>> Error: did not recognize filetype.\n")

        Ainteractors, Binteractors = jobsetup.getinteractors(table, filetype, columnA, columnB)

        jobsetup.getfastas_writecommands(Ainteractors, Binteractors, consideruniprot, considerstart, considerend, split=True,fraglen=fraglen,overlap=overlap,dimerize=dimerize,write=towrite,alphafold_exec=alphafold_exec)

    ##################################

