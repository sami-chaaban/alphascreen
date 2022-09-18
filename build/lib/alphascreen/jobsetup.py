from unipressed import UniprotkbClient
import pandas as pd
import os
from pathlib import Path
import sys

def getuni(ACCESSION):
    ENTRY = UniprotkbClient.fetch_one(ACCESSION)
    SEQ = ENTRY['sequence']['value']
    NAME = ENTRY['uniProtkbId']
    return(NAME, SEQ)

def getinteractors(file,filetype, columnA, columnB):
    
    if filetype == "table":
        table = pd.read_table(file)
    elif filetype == "excel":
        table = pd.read_excel(file)
    
    Ainteractors = table[columnA].tolist()
    Binteractors = table[columnB].tolist()
    
    return(Ainteractors, Binteractors)

def getfastas_writecommands(Ainteractors, Binteractors, consideruniprot,considerstart,considerend, split=True,fraglen=500,overlap=50,dimerize="",write=True,alphafold_exec="colabfold2"):
    
    print("\nStarting...\n")
    
    if dimerize != "":
        print(">> Dimerizing uniprot " + dimerize + "...\n")

    if consideruniprot != "":
        print(">> Only considering " + str(considerstart+1) + "-" + str(considerend+1) + " for uniprot ID"+consideruniprot+"...\n")
    
    #remove duplicates
    repeating_indices = []
    for i, pair in enumerate(zip(Ainteractors, Binteractors)):
        for j, pair_same in enumerate(zip(Ainteractors, Binteractors)):
            if pair[0] == pair_same[0] and pair[1] == pair_same[1] and i != j:
                repeating_indices.append(j)
                #print(pair[0], pair[1], pair_same[0], pair_same[1])
            if pair[0] == pair_same[1] and pair[1] == pair_same[0]:
                repeating_indices.append(j)
                #print(pair[0], pair[1], pair_same[0], pair_same[1])
    Ainteractors = [j for i, j in enumerate(Ainteractors) if i not in repeating_indices]
    Binteractors = [j for i, j in enumerate(Binteractors) if i not in repeating_indices]
    
    print(">> Removed " + str(len(repeating_indices)) + " duplicates...\n")
    
    os.makedirs('fastas', exist_ok=True)
    os.makedirs('results', exist_ok=True)

    colabfoldcommand = []
    unilst = []
    Anamefail = []
    Bnamefail = []

    for A, B in zip(Ainteractors, Binteractors):

        try:
            Aname,Aseq=getuni(A)
        except:
            Anamefail.append(Aname)
            continue
        try:
            Bname,Bseq=getuni(B)
        except:
            Anamefail.append(Aname)
            continue
            
        unilst.append(A + ":" + Aname)
        unilst.append(B + ":" + Bname)

        addmeA=0
        addmeB=0
        if consideruniprot != "":
            if consideruniprot == A:
                if considerend > len(Aseq):
                    sys.exit(">> Error: last amino acid to consider is outside the protein length. Maximum is " + str(len(Aseq)+1) + "\n")
                Aseq = Aseq[considerstart:considerend]
                addmeA = considerstart
            if consideruniprot == B:
                if considerend > len(Bseq):
                    sys.exit(">> Error: last amino acid to consider is outside the protein length. Maximum is " + str(len(Bseq)+1) + "\n")
                Bseq = Bseq[considerstart:considerend]
                addmeB = considerstart

        todimer = ""
        if dimerize != "":
            if dimerize == A:
                todimer = "A"
            elif dimerize == B:
                todimer = "B"
            with open("dimerized-"+A+"-"+Aname, 'w') as f:
                f.write("")
        
        if split:
            #print(Aname)
            Anamelst, Aseqlst = splitfasta(Aname,Aseq,fraglen,overlap,addmeA)
            #print(Bname)
            Bnamelst, Bseqlst = splitfasta(Bname,Bseq,fraglen,overlap,addmeB)
        else:
            Anamelst = [Aname]
            Aseqlst = [Aseq]
            Bnamelst = [Bname]
            Bseqlst = [Bseq]
            
        #print(Anamelst, Aseqlst)
        
        for An, As in zip(Anamelst, Aseqlst):
            
            for Bn, Bs in zip(Bnamelst, Bseqlst):

                fastaname = 'fastas/'+An+'-'+Bn+'.fasta'
                fastaname_minimal = 'fastas/'+An+'-'+Bn+'.fasta'
                
                if write:
                    
                    with open(fastaname, 'w') as f:
                        
                        f.write(">"+An+"\n"+As+"\n")
                        
                        if todimer == "A":
                            f.write(">"+An+"\n"+As+"\n")
                            
                        f.write(">"+Bn+"\n"+Bs+"\n")
                        
                        if todimer == "B":
                            f.write(">"+Bn+"\n"+Bs+"\n")

                colabfoldcommand.append(alphafold_exec + " " + fastaname_minimal + " --output=results/" + An + "-" + Bn)
    
    failed = Anamefail + Bnamefail
    if len(failed) > 0:
        print(">> Warning: there was a problem getting the Uniprot data for accessions: ")
        for f in failed:
            print(f)

    colabfoldcommand=list(dict.fromkeys(colabfoldcommand)) 

    if write:
        print(">> Writing " + str(len(colabfoldcommand)) + " colabfold commands...\n")
        with open("colabfoldrun.bsh", 'w') as f:
            for c in colabfoldcommand:
                f.write(c+"\n")
                
        with open("uniprots.txt", 'w') as f:
            for u in list(set(unilst)):
                f.write(u+"\n")

        print(">> Run the alphafold jobs with bash \"colabfoldrun.bsh\"")
    else:
        print(">> There are " + str(len(colabfoldcommand)) + " colabfold commands...\n")
        
            
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
            print("\n>> Warning: could not find one of the directories. Make sure this the right directory where the colabfold jobs were run.")
            warned=True
        
        found = False
        if os.path.exists(resultdir):
            for fname in os.listdir(resultdir):
                if fname.endswith('.pdb'):
                    found = True
                    numfound+=1
                    break
        if not found:
            colabfoldcommand.append(alphafold_exec + " fastas/" + str(fastapath).split("fastas/")[1] + " --output=results" + str(resultdir).split("results")[1] + "\n")

                
    print("\n>> Finished runs: " + str(numfound) + "\n>> Unfinished runs: " + str(len(colabfoldcommand))+"\n")
    
    colabfoldcommand=list(dict.fromkeys(colabfoldcommand)) 

    if write:
        with open("colabfoldrun-unfinished.bsh", 'w') as f:
            for c in colabfoldcommand:
                f.write(c)
