import json
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
from textwrap import wrap
from pymol import cmd
import cv2 as cv
import glob
import sys
import time
import pymol
import math

def getscores(rankby):

    proteinA_lst=[]
    proteinB_lst=[]
    uniprotA_lst=[]
    uniprotB_lst=[]
    iptmhighscore_lst = []
    ptmhighscore_lst = []
    mincrosspae_lst = []
    inv_mincrosspae_lst = []
    pdbname_lst = []
    paepng_lst = []
    paejson_lst = []
    
    existuniprot = False
    #get uniprots
    if os.path.exists("uniprots.txt"):
        existuniprot = True
    uniprotlst = []
    if existuniprot:
        with open("uniprots.txt") as f:
            lines=f.readlines()
        for l in lines:
            uniprotlst.append([l.split(":")[0], l.split(":")[1].replace("\n","")])
        #print(uniprotlst)

    skipped=0

    for resultdir in Path('results').glob('*/'):

        resultdir = str(resultdir)
        resultname = resultdir.split("/")[-1]
            
        try:
            next(Path(resultdir).glob("*rank*.pdb"))
        except StopIteration:
            skipped+=1
            #print("\n>> Warning: skipping " + resultname + " since it is missing PDB files.")
            continue
            
        try:
            next(Path(resultdir).glob("*_scores.json"))
            scorewildcard="_scores.json"
        except StopIteration:
            try:
                next(Path(resultdir).glob("*_pae.json"))
                scorewildcard="_pae.json"
            except StopIteration:
                try:
                    next(Path(resultdir).glob("*predicted_aligned_error_v1.json"))
                    scorewildcard="predicted_aligned_error_v1.json"
                except StopIteration:
                    print("\n>> Warning: skipping " + resultname + " since it is missing the scores json files.")
                    skipped+=1
                    continue

        modelnumlst = []
        pdbpaths = []
        for path in Path(resultdir).glob("*rank*.pdb"):
            pdbpaths.append(str(path))
            modelnum = str(path).split("model_")[-1].split(".pdb")[0].split("_")[0]
            modelnumlst.append(int(modelnum)-1) #-1 is to turn it into index (0-4) since numbers are 1-5

        if os.path.exists(resultdir+"/scores.txt"):
            with open(resultdir+"/scores.txt") as f:
                lines=f.readlines()
                iptms=[float(x.split(' ')[0].split(":")[1]) for x in lines]
                ptms=[float(x.split(' ')[1].split(":")[1]) for x in lines]
        else:
            iptms=["-"] * len(pdbpaths)
            ptms=["-"] * len(pdbpaths)
            if rankby in ["iptm, ptm"]:
                print("\n>> Warning: skipping " + resultname + " since it is missing the scores.txt file with the iptms and ptms.")
                skipped+=1
                continue

        try:
            next(Path(resultdir).glob("*PAE*.png"))
            for path in Path(resultdir).glob("*PAE*.png"):
                paepng_lst.append(str(path))
                break
        except StopIteration:
            paepng_lst.append("-")
        
        ProteinA=resultname.split("-")[0]
        ProteinB=resultname.split("-")[1]

        if len(ProteinA.split("_")) == 4: #Those from uniprot have "NAME_HUMAN_1_100" or similar
            ProteinAmin=ProteinA.split("_")[0]+"_"+ProteinA.split("_")[1]
        else: #not from uniprot
            ProteinAmin=ProteinA.split("_")[0]
        if len(ProteinB.split("_")) == 4: #Those from uniprot have "NAME_HUMAN_1_100" or similar
            ProteinBmin=ProteinB.split("_")[0]+"_"+ProteinB.split("_")[1]
        else: #not from uniprot
            ProteinBmin=ProteinB.split("_")[0]

        proteinA_lst.append(ProteinA)
        proteinB_lst.append(ProteinB)
         
        paejsons=[]
        paenumlst=[]
        if scorewildcard == "_scores.json":
            for path in Path(resultdir).glob("*_model_*_scores.json"):
                paejsons.append(str(path))
                paenumlst.append(int(str(path).split("_scores.json")[0].split("_model_")[-1].split("_seed")[0])-1)
        elif scorewildcard == "_pae.json":
            for path in Path(resultdir).glob("rank_*_model_*_ptm_seed_0_pae.json"):
                paejsons.append(str(path))
                paenumlst.append(int(str(path).split("_ptm_seed_0_pae.json")[0].split("_model_")[-1])-1)
        elif scorewildcard == "predicted_aligned_error_v1.json":
            for path in Path(resultdir).glob("*_rank_*_model_*_seed_000.json"):
                paejsons.append(str(path))
                paenumlst.append(int(str(path).split("_seed_000.json")[0].split("_model_")[-1])-1)
        
        if existuniprot:
            foundA=False
            foundB=False
            for u in uniprotlst:
                if u[1] == ProteinAmin and not foundA:
                    uniprotA_lst.append(u[0])
                    foundA=True
                if u[1] == ProteinBmin and not foundB:
                    uniprotB_lst.append(u[0])
                    foundB=True
                if foundA and foundB:
                    break
            if not foundA:
                print("Could not find uniprot ID for " + ProteinAmin)
            if not foundB:
                print("Could not find uniprot ID for " + ProteinBmin)


        foundflag=False
        if rankby in ["iptm", "ptm"]:
            if rankby == "iptm":
                m = np.argmax(iptms)
            else:
                m = np.argmax(ptms)

        elif rankby=="scaledPAE":
            minpaes=[]
            rankfile = resultdir+"/bestpaerank.txt"
            if os.path.exists(rankfile) and os.stat(rankfile).st_size != 0:
                with open(rankfile) as f:
                    lines = [line for line in f.readlines()]
                if ":" in lines[-1]:
                    m = int(lines[-1].split(":")[-1])-1
                    minpaes = [line.strip() for line in lines[:-1]]
                    minpaes = [float(i.split("=")[-1]) for i in minpaes]                   
                    foundflag=True

            if not foundflag:
                with open(rankfile, 'w') as f:
                    for paenum, paejson in zip(paenumlst, paejsons):
                        minpae = getmincrosspae(paejson)
                        minpaes.append(minpae)
                        f.write("model"+str(paenum+1)+"="+str(minpae)+"\n")
                    m = paenumlst[np.argmin(minpaes)]
                    f.write("bestmodel:"+str(m+1))
            # For the future: if same crosspae, choose best iptm
            # if len([k for k in minpaes if k == min(minpaes)]) > 1:
            #     print("there are multiple solutions for " + str(scorepath))
        
        else:
            sys.exit("\n>> Provide ptm, iptm, or pae as the ranking method.\n")

        iptmhighscore_lst.append(iptms[m])
        ptmhighscore_lst.append(ptms[m])

        paejson_lst.append(paejsons[paenumlst.index(m)])
        if rankby=="scaledPAE":
            mincrosspae_lst.append(minpaes[paenumlst.index(m)])
            inv_mincrosspae_lst.append(1-(minpaes[paenumlst.index(m)]/30))
        else:
            mincrosspae_lst.append("-")
            inv_mincrosspae_lst.append("-")

        pdbname = pdbpaths[modelnumlst.index(m)]
        pdbname_lst.append(pdbname)
                
    if not existuniprot:
        uniprotA_lst = ["-"] * len(proteinA_lst)
        uniprotB_lst = ["-"] * len(proteinB_lst)
            
    #print(len(proteinA_lst), len(proteinB_lst), len(iptmhighscore_lst), len(ptmhighscore_lst), len(pdbname_lst), len(paepng_lst), len(paejson_lst))

    df = pd.DataFrame(
        {'Protein A': proteinA_lst,
         'Protein B': proteinB_lst,
         'iptm': iptmhighscore_lst,
         'ptm': ptmhighscore_lst,
         'minPAE': mincrosspae_lst,
         'scaledPAE': inv_mincrosspae_lst,
         'Model': pdbname_lst,
         'PAE-png': paepng_lst,
         'PAE-json': paejson_lst,
         'SWISS-PROT Accessions Interactor A' : uniprotA_lst,
         'SWISS-PROT Accessions Interactor B' : uniprotB_lst
        })
        
    df.sort_values(by=[rankby], ascending=False, ignore_index=True, inplace=True)

    if skipped > 0:
        print("\n>> Warning: skipped " + str(skipped) + " since they were missing necessary files.")

    return(df)

def make_paeplot(pae, protnames, protlens):
    
    plt.imshow(pae, vmin=0, vmax=30, cmap="bwr")
    currlen=0
    yticklocs = []
    newprotnames = []
    for i, p in enumerate(protnames):
        if len(p.split("_")) == 4: #Those from uniprot have "NAME_HUMAN_1_100" or similar
            newprotnames.append(p.split("_")[0]+"\n"+p.split("_")[2]+"-"+p.split("_")[3])
        else:
            newprotnames.append(p.split("_")[0]+"\n"+p.split("_")[1]+"-"+p.split("_")[2])
        
    for n, l in zip(protnames, protlens):
        
        yticklocs.append(int(l/2)+currlen)
        
        currlen = currlen + l
        
        plt.axvline(x=currlen-1, color='k')
        plt.axhline(y=currlen-1, color='k')
    
    plt.xlim(1, currlen-1)
    plt.ylim(1, currlen-1)
    plt.yticks(yticklocs, newprotnames, rotation='horizontal')
    plt.gca().invert_yaxis()
    
    return(plt)
    

def getinfo_frompaename(paejson):

    scores = json.loads(Path(paejson).read_text())

    pae = scores["pae"]            

    nameA = str(paejson).split("results/")[1].split("/")[0].split("-")[0]
    nameB = str(paejson).split("results/")[1].split("/")[0].split("-")[1]
    fasta=str(paejson).split("results/")[0]+"fastas/"+nameA+"-"+nameB+".fasta"
    with open(fasta, 'r') as f:
        lines = f.readlines()
    protnames = []
    protlens = []
    for i, a in enumerate(lines):
        if a[0]==">":
            protnames.append(a[1:])
            protlens.append(len(lines[i+1]))
    #lenprotA = int(nameA.split("_")[-1]) - int(nameA.split("_")[-2]) + 1
    #lenprotB = int(nameB.split("_")[-1]) - int(nameB.split("_")[-2]) + 1
    #nameprotA = nameA.split("_")[0] + " " + nameA.split("_")[1] + "\n" + nameA.split("_")[2]+ "-" + nameA.split("_")[3]
    #nameprotB = nameB.split("_")[0] + " " + nameB.split("_")[1] + "\n" + nameB.split("_")[2]+ "-" + nameB.split("_")[3]
    
    return(pae, protnames, protlens)

def getinfo_frompaename_legacy(paejson):

    scores = json.loads(Path(paejson).read_text())

    pae = scores[0]["distance"]
    splitby = int(math.sqrt(len(pae)))
    pae = [pae[i:i + splitby] for i in range(0, len(pae), splitby)]

    nameA = str(paejson).split("results/")[1].split("/")[0].split("-")[0]
    nameB = str(paejson).split("results/")[1].split("/")[0].split("-")[1]
    fasta=str(paejson).split("results/")[0]+"fastas/"+nameA+"-"+nameB+".fasta"
    with open(fasta, 'r') as f:
        lines = f.readlines()
    protnames = []
    protlens = []
    for i, a in enumerate(lines):
        if a[0]==">":
            protnames.append(a[1:])
            protlens.append(len(lines[i+1]))

    #lenprotA = int(nameA.split("_")[-1]) - int(nameA.split("_")[-2]) + 1
    #lenprotB = int(nameB.split("_")[-1]) - int(nameB.split("_")[-2]) + 1
    #nameprotA = nameA.split("_")[0] + " " + nameA.split("_")[1] + "\n" + nameA.split("_")[2]+ "-" + nameA.split("_")[3]
    #nameprotB = nameB.split("_")[0] + " " + nameB.split("_")[1] + "\n" + nameB.split("_")[2]+ "-" + nameB.split("_")[3]
    
    return(pae, protnames, protlens)

# def showpae(paejson):

#     pae, protnames, protlens = getinfo_frompaename(paejson)
#     plt = make_paeplot(pae, protnames, protlens)
#     plt.show()
#     return(len(protnames))

def getpae(paejson):
    try:
        return(getinfo_frompaename(paejson))
    except:
        try:
            return(getinfo_frompaename_legacy(paejson))
        except:
            sys.exit("\n>> Could not interpret the PAE json file.\n")

def getmincrosspae(paejson):

    pae, protnames, protlens = getpae(paejson)

    if len(protnames) == 2:

        subgrid = np.array(pae)[(protlens[0]+2):, 0:(protlens[0]-2)]

    elif len(protnames) == 3:

        if protnames[0] == protnames[1]:
            dimerlen = protlens[0]+protlens[1] # which should be the same
            subgrid = np.array(pae)[(dimerlen+2):, 0:(dimerlen-2)]

        else:
            dimerlen = protlens[1]+protlens[2] # which should be the same
            subgrid = np.array(pae)[(protlens[0]+2):, 0:(protlens[0]-2)]

    elif len(protnames) == 4:

        dimerlen = protlens[0]+protlens[1]
        subgrid = np.array(pae)[(dimerlen+2):, 0:(dimerlen-2)]

    return(subgrid.min())

def summarize_pae_pdf(df, threshold, rankby):
    
    if threshold == 0:
        pdfname = "PAEs-rankedby-" + rankby + ".pdf"
    else:
        pdfname = "PAEs-"+rankby+"-above-" + str(threshold).replace(".", "p") + ".pdf"
    
    print("\n>> Writing " + pdfname)

    with PdfPages(pdfname) as pdf:
        for i, result in df[df[rankby]>threshold].iterrows():
            pae, protnames, protlens = getpae(result["PAE-json"])
            plt = make_paeplot(pae, protnames, protlens)
            plt.title("\n".join(wrap(result["Model"], 80)), fontsize=8)
            pdf.savefig()
            plt.close('all')

# def show_pdb(modelname, chains, show_sidechains, show_mainchains=False, color="chain"):
    
#     view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
#     view.addModel(open(modelname,'r').read(),'pdb')

#     if color == "lDDT":
#         view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
#     elif color == "rainbow":
#         view.setStyle({'cartoon': {'color':'spectrum'}})
#     elif color == "chain":
#         #chains = len(queries[0][1]) + 1 if is_complex else 1
#         for n,chain,color in zip(range(chains),list("BCDEFGH"),
#                          ["cyan", "magenta", "orange", "yellow"]): #"["lime","cyan","magenta","yellow","salmon","white","blue","orange"]"
#             view.setStyle({'chain':chain},{'cartoon': {'color':color}})
#     if show_sidechains:
#         BB = ['C','O','N']
#         view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
#                             {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
#         view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
#                             {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
#         view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
#                             {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})  
#     if show_mainchains:
#         BB = ['C','O','N','CA']
#         view.addStyle({'atom':BB},{'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
        
#     nameA = str(modelname).split("results/")[1].split("/")[0].split("-")[0]
#     nameB = str(modelname).split("results/")[1].split("/")[0].split("-")[1]
#     lenprotA = int(nameA.split("_")[-1]) - int(nameA.split("_")[-2]) + 1
#     lenprotB = int(nameB.split("_")[-1]) - int(nameB.split("_")[-2]) + 1
#     nameprotA = nameA.split("_")[0] + " " + nameA.split("_")[1] + "\n" + nameA.split("_")[2]+ "-" + nameA.split("_")[3]
#     nameprotB = nameB.split("_")[0] + " " + nameB.split("_")[1] + "\n" + nameB.split("_")[2]+ "-" + nameB.split("_")[3]
#     #view.addLabel(nameprotA,{'fontOpacity':1, 'fontSize':12, 'fontColor':'black','backgroundOpacity':0.2, 'backgroundColor':'cyan'},{'resi':1})
#     #view.addLabel(nameprotB,{'fontOpacity':1, 'fontSize':12, 'fontColor':'black','backgroundOpacity':0.2, 'backgroundColor':'magenta'},{'resi':lenprotA+10})

#     view.zoomTo()
#     return view

def write_top(df, threshold, rankby):

    if threshold == 0:
        basename = "Results"
    else:
        basename = "Results-"+rankby+"-above-" + str(threshold).replace(".", "p")
        
    excelname = basename + ".xlsx"
    csvname = basename + ".csv"

    print("\n>> Writing " + excelname + " and " + csvname)
    
    with pd.ExcelWriter(excelname) as writer:  
        df[df[rankby]>threshold].to_excel(writer)
        
    df.to_csv(csvname)
        
        
def write_modelpngs(df, threshold, rankby, overwrite=False):

    print("\n>> Writing model snapshots...")

    pymol.finish_launching(['pymol', '-qc']) #-Q will suppress render outputs too if you want

    total=0

    for i,result in df[df[rankby]>threshold].iterrows():
        
        model=result["Model"]
        png1 = model[:-4]+".png"
        png2 = model[:-4]+"-rotated.png"
        
        if not overwrite and os.path.exists(png1) and os.path.exists(png2):
            continue

        if os.path.exists(png1):
            os.rename(png1, model[:-4]+"-backup.png")
        
        cmd.load(model, "current")
        cmd.cartoon("automatic")
        cmd.bg_color("white")
        #cmd.ray(600,600)
        #cmd.draw(300,300,antialias=2)
        cmd.zoom()
        cmd.util.cbc()
        cmd.png(png1, width=600, height=600, dpi=900)

        waitfor(png1)

        if not os.path.exists(png1):
            if os.path.exists(model[:-4]+"-backup.png"):
                os.rename(model[:-4]+"-backup.png", png1)
            sys.exit("\n>> Error: pymol didn't output anything.\nFailed on: "+model+"\n")
        elif os.path.exists(model[:-4]+"-backup.png"):
            os.remove(model[:-4]+"-backup.png")
        
        cmd.rotate(axis='x',angle=90)
        cmd.rotate(axis='y',angle=90)
        cmd.zoom()
        cmd.png(png2, width=600, height=600, dpi=900)

        waitfor(png2)

        cmd.delete("current")
        
        total+=1
        
def waitfor(filename):
    pngexists=False
    tic = time.perf_counter()
    while not pngexists:
        if os.path.exists(filename):
            pngexists=True
        time.sleep(0.25)
        toc = time.perf_counter()
        if toc-tic > 4:
            sys.exit("\n>> Error: pymol didn't output anything in 4 seconds for file:\n"+filename)
    
def summarize_paeandmodel_pdf(df, threshold, rankby):
    
    if threshold == 0:
        pdfname = "PAEs-Models.pdf"
    else:
        pdfname = "PAEs-Models-"+rankby+"-above-"+ str(threshold).replace(".", "p") + ".pdf"
    
    print("\n>> Writing " + pdfname + "\n")

    with PdfPages(pdfname) as pdf:
        for i, result in df[df[rankby]>threshold].iterrows():
            
            if result["iptm"] == "-":
                iptmstring = "-"
                ptmstring = "-"
            else:
                iptmstring = str("{:.2f}".format(result["iptm"]))
                ptmstring = str("{:.2f}".format(result["ptm"]))
            
            if rankby=="scaledPAE":
                plottitle = [result["Model"],"iptm: "+ iptmstring + ", ptm: "+ ptmstring,"minimum-pae: "+str("{:.2f}".format(result["minPAE"])) + ", scaled-pae: "+str("{:.2f}".format(result["scaledPAE"])), result["Protein A"]+" ("+result["SWISS-PROT Accessions Interactor A"]+")",result["Protein B"]+" ("+result["SWISS-PROT Accessions Interactor B"]+")"]
            else:
                plottitle = [result["Model"],"iptm: "+ iptmstring + ", ptm: "+ iptmstring, result["Protein A"]+" ("+result["SWISS-PROT Accessions Interactor A"]+")",result["Protein B"]+" ("+result["SWISS-PROT Accessions Interactor B"]+")"]
            plottitle = "\n".join(plottitle)
            
            try:
                pae, protnames, protlens = getinfo_frompaename(result["PAE-json"])
            except:
                try:
                    pae, protnames, protlens = getinfo_frompaename_legacy(result["PAE-json"])
                except:
                    sys.exit("\n>> Could not interpret the PAE json file.\n")

            fig = plt.figure(figsize=(60,20))
            fig.suptitle(plottitle, fontsize=30)
            plt.subplot(1, 3, 1)
            make_paeplot(pae, protnames, protlens)
            plt.xticks(fontsize=36)
            plt.yticks(fontsize=36)

            png1 = result["Model"][:-4]+".png"
            png2 = result["Model"][:-4]+"-rotated.png"
            
            if not os.path.exists(png1):
                sys.exit("\n>> Error: The model snapshots could not be found.\n")
            elif not os.path.exists(png2):
                sys.exit("\n>> Error: The model snapshots could not be found.\n")
            
            plt.subplot(1, 3, 2)
            img1=plt.imread(png1)
            plt.axis('off')
            #img1=cv.imread(png1)
            plt.imshow(img1)
            
            plt.subplot(1, 3, 3)
            img2=plt.imread(png2)
            plt.axis('off')
            plt.imshow(img2)
            
            pdf.savefig()
            plt.close('all')

    print(">> If there is an error below, ignore it." + "\n")
