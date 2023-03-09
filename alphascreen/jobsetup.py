from unipressed import UniprotkbClient
import pandas as pd
import os
from pathlib import Path
import sys
import numpy as np
from Bio import SeqIO

def getuni(ACCESSION):
    ENTRY = UniprotkbClient.fetch_one(ACCESSION)
    SEQ = ENTRY['sequence']['value']
    NAME = ENTRY['uniProtkbId']
    return(NAME, SEQ)

def getinteractors(file,filetype, columnA, columnB, focus, exhaustive):
    
    print("\n>> Reading table...\n")
    if filetype == "table":
        table = pd.read_table(file, usecols=[columnA, columnB])
    elif filetype == "excel":
        table = pd.read_excel(file, usecols=[columnA, columnB])
    
    Ainteractors = table[columnA].tolist()
    Binteractors = table[columnB].tolist()

    #cleanup
    if exhaustive:
        Ainteractors = [x for x in Ainteractors if x == x] #remove nan float which result from empty cells
        Binteractors = [x for x in Binteractors if x == x]
    else:
        Ainteractors = ["empty-cell" if x != x else x for x in Ainteractors] #replace nan-float with 0, which result from empty cells
        Binteractors = ["empty-cell" if x != x else x for x in Binteractors]
    Ainteractors = [s.strip() for s in Ainteractors]
    Binteractors = [s.strip() for s in Binteractors]

    if exhaustive:
        Ainteractors = list(set(Ainteractors))
        Binteractors = list(set(Binteractors))
        Ainteractors_exh = []
        Binteractors_exh = []
        for A in Ainteractors:
            for B in Binteractors:
                Ainteractors_exh.append(A)
                Binteractors_exh.append(B)
        Ainteractors = Ainteractors_exh
        Binteractors = Binteractors_exh

    if focus != "":
        if focus not in Ainteractors and focus not in Binteractors:
            sys.exit("!! Error: the uniprot ID " + focus + " to focus on was not found.\n")
        Ainteractors_fixed = []
        Binteractors_fixed = []
        for A, B in zip(Ainteractors,Binteractors):
            if B == focus:
                Ainteractors_fixed.append(B)
                Binteractors_fixed.append(A)
            else:
                Ainteractors_fixed.append(A)
                Binteractors_fixed.append(B)
        Ainteractors = Ainteractors_fixed
        Binteractors = Binteractors_fixed

    return(Ainteractors, Binteractors)

def getfastas_writecommands(Ainteractors, Binteractors, consideruniprot,considerstart,considerend, customids, split=True,fraguniprot=0,fraglen=500,overlap=50,dimerize="",dimerize_all=False,dimerize_except="",write=True,alphafold_exec="colabfold2",ignoreself=False):
    
    dimerizelst = []
    dontdimerizelst = []

    if dimerize_all:
        print(">> Dimerizing all proteins...\n")

    elif dimerize_except:
        if dimerize_except[-4:] == ".txt":
            print(">> Dimerizing all proteins except those in " + dimerize_except + "...\n")
            with open(dimerize_except) as f:
                dontdimerizelst=[line.strip() for line in f.readlines()]
        else:
            dontdimerizelst = [dimerize_except]
            print(">> Dimerizing all proteins except " + dimerize_except + "...\n")

    elif dimerize != "":
        if dimerize[-4:] == ".txt":
            print(">> Dimerizing uniprot IDs in " + dimerize + "...\n")
            with open(dimerize) as f:
                dimerizelst=[line.strip() for line in f.readlines()]
        else:
            dimerizelst = [dimerize]
            print(">> Dimerizing uniprot ID " + dimerize + "...\n")

    #remove exact duplicates
    rawtotal = len(Ainteractors)
    Ainteractors, Binteractors = zip(*list(dict.fromkeys(zip(Ainteractors, Binteractors))))
    
    #remove reversed duplicates
    repeating_indices = []
    for i, pair in enumerate(zip(Ainteractors, Binteractors)):
        for j, pair_same in enumerate(zip(Ainteractors, Binteractors)):
            if pair[0] == pair_same[1] and pair[1] == pair_same[0] and j > i:
                repeating_indices.append(j)
                #print(pair[0], pair[1], pair_same[0], pair_same[1])
    Ainteractors = [j for i, j in enumerate(Ainteractors) if i not in repeating_indices]
    Binteractors = [j for i, j in enumerate(Binteractors) if i not in repeating_indices]

    newtotal = len(Ainteractors)
    numremoved = rawtotal - newtotal
    if len(Ainteractors) > 1 and len(Binteractors) > 1:
        print(">> Removed " + str(numremoved) + " duplicated pairs...\n")

    #remove self
    if ignoreself:  
        self_indices = []
        for i, pair in enumerate(zip(Ainteractors, Binteractors)):
            if pair[0] == pair[1]:
                self_indices.append(i)
        Ainteractors = [j for i, j in enumerate(Ainteractors) if i not in self_indices]
        Binteractors = [j for i, j in enumerate(Binteractors) if i not in self_indices]
        numremoved = newtotal - len(Ainteractors)
        if len(Ainteractors) > 1 and len(Binteractors) > 1:
            print(">> Removed " + str(numremoved) + " pairs where the proteins are the same...\n")

    if write:
        os.makedirs('fastas', exist_ok=True)
        os.makedirs('results', exist_ok=True)

    for fname in os.listdir('fastas'):
        if fname.endswith('.fasta'):
            sys.exit(">> Error: it looks like fasta files already exist in the fastas folder.\n")

    colabfoldcommand = []
    unilst = []
    failed = []

    if len(dontdimerizelst) != 0:
        alluniprots = list(set(list(Ainteractors+Binteractors)))
        dimerizelst = [i for i in alluniprots if i not in dontdimerizelst]

    print(">> Getting sequences...\n")

    customid_names = []
    customid_seqs = []
    if customids != "":
        with open(customids, mode='r') as handle:
             for record in SeqIO.parse(handle, 'fasta'):
                    customid_names.append(record.id)
                    customid_seqs.append(str(record.seq))

    Unis_found = []
    Names_found = []
    Seqs_found = []

    notconsidered = consideruniprot.copy() #fully populated and will be removed one by one to check that all were found
    for A, B in zip(Ainteractors, Binteractors):

        try:
            if A in Unis_found:
                Aindex = Unis_found.index(A)
                Aname = Names_found[Aindex]
                Aseq = Seqs_found[Aindex]
            else:
                if A in customid_names:
                    Aindex = customid_names.index(A)
                    Aname = customid_names[Aindex]
                    Aseq = customid_seqs[Aindex]
                else:
                    Aname,Aseq=getuni(A)
                Unis_found.append(A)
                Names_found.append(Aname)
                Seqs_found.append(Aseq)
        except:
            failed.append(A)
            continue
        try:
            if B in Unis_found:
                Bindex = Unis_found.index(B)
                Bname = Names_found[Bindex]
                Bseq = Seqs_found[Bindex]
            else:
                if B in customid_names:
                    Bindex = customid_names.index(B)
                    Bname = customid_names[Bindex]
                    Bseq = customid_seqs[Bindex]
                else:
                    Bname,Bseq=getuni(B)
                Unis_found.append(B)
                Names_found.append(Bname)
                Seqs_found.append(Bseq)
        except:
            failed.append(B)
            continue
            
        unilst.append(A + ":" + Aname)
        unilst.append(B + ":" + Bname)

        addmeA=0
        addmeB=0
        if len(consideruniprot) > 0:

            if A in consideruniprot:

                a = consideruniprot.index(A)
                if considerend[a] > len(Aseq):
                    sys.exit(">> Error: last amino acid to consider is outside the protein length for " + consideruniprot[a] + ". Maximum is " + str(len(Aseq)+1) + ".\n")
                Aseq = Aseq[considerstart[a]:considerend[a]]
                addmeA = considerstart[a]

                if A in notconsidered:
                    notconsidered.remove(A) 

            if B in consideruniprot:

                b = consideruniprot.index(B)
                if considerend[b] > len(Bseq):
                    sys.exit(">> Error: last amino acid to consider is outside the protein length for " + consideruniprot[b] + ". Maximum is " + str(len(Bseq)+1) + ".\n")
                Bseq = Bseq[considerstart[b]:considerend[b]]
                addmeB = considerstart[b]

                if B in notconsidered:
                    notconsidered.remove(B) 
        
        if split:

            if type(fraglen) == list:
                if A in fraguniprot:
                    fraglen_cur=fraglen[fraguniprot.index(A)]
                else:
                    fraglen_cur = 500
            else:
                fraglen_cur = fraglen

            Anamelst, Aseqlst = splitfasta(Aname,Aseq,fraglen_cur,overlap,addmeA)

            if type(fraglen) == list:
                if B in fraguniprot:
                    fraglen_cur=fraglen[fraguniprot.index(B)]
                else:
                    fraglen_cur = 500
            else:
                fraglen_cur = fraglen

            Bnamelst, Bseqlst = splitfasta(Bname,Bseq,fraglen_cur,overlap,addmeB)

        else:

            Anamelst = [Aname]
            Aseqlst = [Aseq]
            Bnamelst = [Bname]
            Bseqlst = [Bseq]
            
        
        for An, As in zip(Anamelst, Aseqlst):
            
            for Bn, Bs in zip(Bnamelst, Bseqlst):

                fastaname = 'fastas/'+An+'-'+Bn+'.fasta'
                fastaname_minimal = 'fastas/'+An+'-'+Bn+'.fasta'
                
                if write:
                    
                    with open(fastaname, 'w') as f:
                        
                        f.write(">"+An+"\n"+As+"\n")
                        
                        if A in dimerizelst or dimerize_all:
                            f.write(">"+An+"\n"+As+"\n")
                            
                        f.write(">"+Bn+"\n"+Bs+"\n")
                        
                        if B in dimerizelst or dimerize_all:
                            f.write(">"+Bn+"\n"+Bs+"\n")

                colabfoldcommand.append(alphafold_exec + " " + fastaname_minimal + " --output=results/" + An + "-" + Bn)

    if len(failed) > 0:
        failed = list(set(failed))
        print("!! Warning: there was a problem getting the Uniprot data for accessions:\n")
        for f in failed:
            print("·· "+ str(f) + "\n")

    if len(notconsidered) > 0:
        for c in notconsidered:
            print("!! Warning: you asked to consider a region within " + c + ", but this ID wasn't in your input.\n")

    colabfoldcommand=list(dict.fromkeys(colabfoldcommand)) 

    if write:
        print(">> Writing " + str(len(colabfoldcommand)) + " Alphafold commands...\n")
        with open("runpredictions.bsh", 'w') as f:
            for c in colabfoldcommand:
                f.write(c+"\n")
                
        with open("uniprots.txt", 'w') as f:
            for u in list(set(unilst)):
                f.write(u+"\n")

        print(">> Run the Alphafold jobs with \"bash runpredictions.bsh\"\n")
    else:
        print(">> There are " + str(len(colabfoldcommand)) + " Alphafold commands...\n")
            
    print("Done!\n")

def splitfasta(name,seq,fraglen,overlap,addme):

    seqlen = len(seq)
    #print(name, seqlen)
    
    splitinto = round(seqlen/fraglen)
    #print(splitinto)
    
    test=[]
    if splitinto in [0, 1]:
        index_withOverlap = [[0,seqlen]]
    else:
        currfraglen = int(seqlen/splitinto)
        middlefragindices=[]
        for n in range(1,splitinto):
            if n == 1 and splitinto > 2:
                middlefragindices.append(currfraglen*n+int(overlap/2)) #account for overlap
            elif n == splitinto-1 and splitinto > 2:
                middlefragindices.append(currfraglen*n-int(overlap/2)) #account for overlap
            else:
                middlefragindices.append(currfraglen*n) 
        
        #print(middlefragindices)
    
        index_withOverlap=[]
        for i,m in enumerate(middlefragindices):
            if i == 0:
                index_withOverlap.append([0, m+overlap])
            else:
                index_withOverlap.append([middlefragindices[i-1]-overlap, m+overlap])
                
        index_withOverlap.append([middlefragindices[-1]-overlap, seqlen])
                    
        #print(index_withOverlap)
        
    splitlens_withOverlap = []
    for i in index_withOverlap:
        splitlens_withOverlap.append(i[1]-i[0])
    #print("New lengths", splitlens_withOverlap)
        
    splitseqs_withOverlap = []
    newnames = []
    
    for i in index_withOverlap:
        
        splitseqs_withOverlap.append(seq[i[0]:i[1]])
        newnames.append(name+"_"+str(i[0]+1+addme)+"_"+str(i[1]+1+addme))
        
    #print(newnames)
        
    return(newnames, splitseqs_withOverlap)

def findunfinished(alphafold_exec, write=True):
    
    colabfoldcommand = []
    numfound=0

    warned=False    
    for fastapath in Path('fastas/').rglob('**/*.fasta'):
        
        resultdir = "results" + str(fastapath).split("fastas")[1][:-6]

        if not os.path.exists(resultdir) and not warned:
            print("\n!! Warning: could not match one or more fasta files with its result folder.\n!! Wait for any running jobs to finish, then run alphascreen --write_unfinished.")
            warned=True
        
        found = False
        if os.path.exists(resultdir):
            for fname in os.listdir(resultdir):
                if fname.endswith('.pdb') and len(fname.split("rank")) == 2:
                    found = True
                    numfound+=1
                    break
        if not found:
            colabfoldcommand.append(alphafold_exec + " fastas/" + str(fastapath).split("fastas/")[1] + " --output=results" + str(resultdir).split("results")[1] + "\n")

                
    print("\n>> Finished runs: " + str(numfound) + "\n>> Unfinished runs: " + str(len(colabfoldcommand))+"\n")
    
    colabfoldcommand=list(dict.fromkeys(colabfoldcommand)) 

    if write:
        if len(colabfoldcommand) == 0:
            sys.exit("\n>> There are no jobs left.\n")
        with open("runpredictions-unfinished.bsh", 'w') as f:
            for c in colabfoldcommand:
                f.write(c)
        print(">> Wrote runpredictions-unfinished.bsh\n")
