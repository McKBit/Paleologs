#Plugins we are using
import shutil, glob, os, csv, numpy, scipy.stats, sys, math

#Get list file name input. Example command: python3 Paleolog.Sort.py "My.CDS.List.File.txt"
file = sys.argv[1]
f = open(f"{file}", 'r' ).read()
f = f.split("\n")

#################################################################
############Do the analysis for each file in the list############
#################################################################
for item in f:
    print(f"------------------------------------------{item}----------------------------------------------")
    #Get Names of files
    CDS = f"{item}"
    Candidates = f"{item}.Can.csv"
    Syntenny = f"{item}.Syn.txt"
    Peak = f"{item}.Peak.csv"

    #Make lists we will be using later on for this species
    CDSList = []
    SYNList = []
    KSList = []
    CANList = []
    PEAKList= []
    BOOKList = []

    #Featch the peak list file and make a list with each peak and it's info
    print("Reading All Peaks in Peak File...")
    f = open(f"Peak.Set/{Peak}", 'r' ).read()
    f = f.split("\n")
    for i in f:
        i = i.split(",")
        PEAKList.append(i)
    del PEAKList[-1]


    #Make a list of the genes in the CDS file and assign an index number to it (ex: >GeneABCD1 is now "1", GeneABCD2 is now "2")
    print("Making list of genes in CDS file...")
    f = open(f"CDS.Files/{CDS}", 'r' ).read()
    f1 = f.split(">")
    x=0
    for i in f1:
        if i == "": # this makes sure we always start at the first gene, sometimes the first item is blank!!!!
            pass
        else:
            x=x+1
            i = i.split("\n")
            i = i[0]
            i = (i, x)
            CDSList.append(i)

    #Make a list of paleologs from Synteny file.
    print("Reading Synteny File...")
    f3 = open(f"CoGe.Synteny.Files/{Syntenny}", 'r' ).read()
    f3 = f3.split("\t")
    f32 = []
    SYNList=[]
    for item in f3:
        if "||" in item:
            item = item.split("||")
            abc = len(item)
            if abc > 8:
                z = item[3]
                f32.append(z)
            else:
                pass

    # Removes duplicated items!!!
    SYNListA = list(set(f32))

    #Translate synteny gene names to CDS index numbers
    a = int(len(SYNListA))
    track=0
    for i in SYNListA:
        track=track+1
        tracker = math.ceil((track/a)*100)
        print(f'Translating Synteny Gene Names... {tracker}% complete', end='\r')
        for b in CDSList:
            if i in b[0]:
                SYNList.append(int(b[1]))

    #Total number of syntenic paleologs
    SynPal = len(SYNList)

    #Print error if Synteny genes are missing or not found in CDSList. This is is a big step to check if the gene annotaiton versions are the same!!!
    warn = (len(SYNListA) - SynPal)
    if warn > 0:
        print(f"!!!!!!!!!!! WARNING {warn} Sytenic Genes Missing From Analysis !!!!!!!!!!!!!  ")

    #Read in the DupPipe KS output file and make a list for it, translating it to the propper code.
    print("Making list of KS values...")
    f2 = open(f"DupPipe.Ks.Files/pamloutput_{CDS}", 'r' ).read()
    f2 = f2.split("\n")
    for line in f2[1:]:
        l = len(line)
        if l >= 2:
            line = line.split("\t")
            KSList.append(line)

    #Make a list of all the gene numbers to be searched based on the canaididate gene lists
    f3 = open(f"Candidate.Lists/{Candidates}", 'r' ).read()
    f3 = f3.split("\n")
    warn1 = 0
    warn2 = 0
    for i in f3:
        warn1=warn1 + 1
        for g in CDSList:
            if i in g[0]:
                warn2 = warn2 + 1
                temp = (i,g[1])
                CANList.append(temp)
                break
            else:
                pass

    del CANList[-1] #there is a blank candidate in last position that needs removing

    #Print error if Candidate genes are missing or not found in CDSList. This is a frequent error if you use other people's gene lists, as annotation versions change and lead ot missing candidates. Either find the origional annotation or move forward knowing you are missing a candidate!!!
    warn = warn2 - warn1
    if warn > 0:
        print(f"!!!!!!!!!!! WARNING {warn} Candidate Genes Missing From Analysis !!!!!!!!!!!!!  ")

    #Make the final lists to be printed to files.
    CandidateList = []
    PeakOutput = []
    names = ["Peak","Max","Min","Total Genes","Total Candidates","Synteny Paleologs","Ks Paleologs","Method","Paleologs From Method","Non-Paleologs","Candidate Paleologs","Candidate Non-Paleologs","Chi.Sqr","P.val"]
    PeakOutput.append(names)

    #####################################################
    ### For Each Peak Begin Sorting all the Paleologs ###
    #####################################################
    x=0
    for thing in PEAKList:
        x=x+1
        max =float(thing[2])
        min = float(thing[3])
        if len(thing) > 4:
            prior = thing[4]
            print(prior)
        else:
            prior = "no"
        print(f"----------------------------- Begin Analysis Peak {x} -----------------------------")

        #########Find all the gene pairs inside the KS peak#########
        print(f"Queering all genes inside KS peak {x}...")
        peak = []
        for i in KSList:
            if i[4] == "":
                i[4] = 0
            Ks = float(i[4])
            if Ks <= max and Ks >= min:
                    peak.append(int(i[0]))
                    peak.append(int(i[1]))
            else:
                pass
        #Remove Duplicates
        peakKS = list(set(peak))
        #Total number of KS based paleologs
        KsPal = len(peakKS)


        ##########Remove those genes which are not in the synteny list aswell#########
        peak = []
        y=0
        track = len(peakKS)
        for i in peakKS:
            y=y+1
            tracker = math.ceil((y/track)*100)
            print(f'Comparing Synteny list to KS list for peak {x}... {tracker}% complete', end='\r')
            for z in SYNList:
                if i == z:
                    peak.append(i)
                else:
                    pass
        print("")
        peak = list(set(peak))


        ##########Find all the genes lower than the KS peak#########
        print(f"Queering all genes bellow KS peak {x}...")
        strict = []
        for i in KSList:
            if i[4] == "":
                i[4] = 0
            Ks = float(i[4])
            if Ks < min:
                    strict.append(int(i[0]))
                    strict.append(int(i[1]))
            else:
                pass

        #Remove duplicates
        PaleologsAStrict = list(set(strict))


        ##########Find all the candidate genes that are inside the ks only peak#########
        print(f"Finding all candidate genes in the KS only peak {x}...")
        CanPaleologs= []
        for gene in CANList:
            genie = int(gene[1])
            for line in peakKS:
                if genie == line:
                        CanPaleologs.append(genie)
                else:
                    pass


        ##########Find all the candidate genes that are inside the peak#########
        print(f"Finding all candidate genes inside peak {x}...")
        CanPaleologsRelaxed= []
        for gene in CANList:
            genie = int(gene[1])
            for line in peak:
                if genie == line:
                        CanPaleologsRelaxed.append(genie)
                else:
                    pass

        ##########Find all the candidate genes that are bellow the peak#########
        print(f"Finding all candidate genes bellow the peak {x}...")
        CanPaleologsAStrict= []
        for gene in CANList:
            genie = int(gene[1])
            for line in strict:
                if genie == line:
                        CanPaleologsAStrict.append(genie)
                else:
                    pass


        ##########Remove Paleologs that have earlier duplications (strict list)#########
        PaleologsStrict = []
        SynDupStrict = []
        for gene in peak:
            if gene not in PaleologsAStrict:
                PaleologsStrict.append(gene)
            else:
                SynDupStrict.append(gene)


        ##########Remove Paleologs that have earlier duplications (strict list)#########
        KSPaleologsStrict = []
        KsDupStrict = []
        for gene in peakKS:
            if gene not in PaleologsAStrict:
                KSPaleologsStrict.append(gene)
            else:
                KsDupStrict.append(gene)

        ##########Remove Candidate Paleologs that have earlier duplications (strict list)#########
        CanPaleologsStrict = []
        SynCanDupStrict = []
        for gene in CanPaleologsRelaxed:
            if gene not in CanPaleologsAStrict:
                CanPaleologsStrict.append(gene)
            else:
                SynCanDupStrict.append(gene)


        ##########Remove Candidate Paleologs that have earlier duplications (strict list)#########
        KSCanPaleologsStrict = []
        KsCanDupStrict = []
        for gene in CanPaleologs:
            if gene not in CanPaleologsAStrict:
                KSCanPaleologsStrict.append(gene)
            else:
                KsCanDupStrict.append(gene)


        ##################Book Keeping Duplicated Ks Only###########################
        #Total Genes
        RDG = int(len(CDSList))
        #Number of Paleologs
        RDP = int(len(KsDupStrict))
        #Number of non-Paleologs
        RDNP = int(RDG - RDP)
        #Total candidates
        RDC = int(len(CANList))
        #Number of candidate Paleologs
        RDCP = int(len(KsCanDupStrict))
        #Number of candidate non-paleologs
        RDCNP = int(len(CANList) - RDCP)

        ##################Book Keeping Dupuplicated Ks + Syn ###########################
        #Total Genes
        SDG = int(len(CDSList))
        #Number of Paleologs
        SDP = int(len(SynDupStrict))
        #Number of non-Paleologs
        SDNP = int(SDG - SDP)
        #Total candidates
        SDC = int(len(CANList))
        #Number of candidate Paleologs
        SDCP = int(len(SynCanDupStrict))
        #Number of candidate non-paleologs
        SDCNP = int(len(CANList) - SDCP)


        ##################Book Keeping KS only Relaxed###########################
        #Number of Paleologs
        P = KsPal
        #Number of non-Paleologs
        NP = int(len(CDSList) - P)
        #Number of candidate Paleologs
        CP = int(len(CanPaleologs))
        #Number of candidate non-paleologs
        CNP = int(len(CANList) - CP)

        ##################Book Keeping KS only Strict###########################

        #Number of  Paleologs
        KSP = int(len(KSPaleologsStrict))
        #Number of non-paleologs
        KSNP = int(len(CDSList) - KSP)
        #Number of candidate Paleologs
        KSCP = int(len(KSCanPaleologsStrict))
        #Number of candidate non-paleologs
        KSCNP = int(len(CANList) - KSCP)

        ##################Book Keeping KS + Synteny Relaxed###########################
        #Total Genes
        RG = int(len(CDSList))
        #Number of Paleologs
        RP = int(len(peak))
        #Number of non-Paleologs
        RNP = int(len(CDSList) - RP)
        #Total candidates
        RC = int(len(CANList))
        #Number of candidate Paleologs
        RCP = int(len(CanPaleologsRelaxed))
        #Number of candidate non-paleologs
        RCNP = int(len(CANList) - RCP)

        ##################Book Keeping KS + Synteny Strict###########################
        #Total Genes
        SG = int(len(CDSList))
        #Number of Paleologs
        SP = int(len(PaleologsStrict))
        #Number of non-Paleologs
        SNP = int(SG - SP)
        #Total candidates
        SC = int(len(CANList))
        #Number of candidate Paleologs
        SCP = int(len(CanPaleologsStrict))
        #Number of candidate non-paleologs
        SCNP = int(len(CANList) - SCP)



        ####################Compute Chisqr for all types###########################


    #Fisher's Exact Test if any value is less than five - KS only Relaxed
        obs = numpy.array([[CP,P], [CNP,NP]])
        chi2, p = scipy.stats.fisher_exact(obs, 'greater')
        chi2 = f"{chi2}"
        p = f"{p}"


    #Fisher's Exact Test if any value is less than five - KS only Strict
        obs = numpy.array([[KSCP,KSP], [KSCNP,KSNP]])
        KSchi2, KSp = scipy.stats.fisher_exact(obs, 'greater')
        KSchi2 = f"{KSchi2}"
        Ksp = f"{KSp}"


    #Fisher's Exact Test if any value is less than five - KS + Synteny Relaxed
        obs = numpy.array([[RCP,RP], [RCNP,RNP]])
        Rchi2, Rp = scipy.stats.fisher_exact(obs, 'greater')
        Rchi2 = f"{Rchi2}"

    #Fisher's Exact Test if any value is less than five - KS + Synteny Strict
        obs = numpy.array([[SCP,SP], [SCNP,SNP]])
        Schi2, Sp = scipy.stats.fisher_exact(obs, 'greater')
        Schi2 = f"{Schi2}"
        sp = f"{Sp}"

    #Fisher's Exact Test if any value is less than five - Dup Relaxed
        obs = numpy.array([[RDCP,RDP], [RDCNP,RDNP]])
        Schi2, Sp = scipy.stats.fisher_exact(obs, 'greater')
        RDchi2 = f"{Schi2}"
        RDp = f"{Sp}"

    #Fisher's Exact Test if any value is less than five - Dup Strict
        obs = numpy.array([[SDCP,SDP], [SDCNP,SDNP]])
        Schi2, Sp = scipy.stats.fisher_exact(obs, 'greater')
        SDchi2 = f"{Schi2}"
        SDp = f"{Sp}"


        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),SG,SC,SynPal,KsPal,"KS Only Relaxed", P,NP,CP,CNP,chi2,p]
        PeakOutput.append(Line)

        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),SG,SC,SynPal,KsPal,"KS Only Strict Paleologs.", KSP,KSNP,KSCP,KSCNP,KSchi2,Ksp]
        PeakOutput.append(Line)

        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),SG,SC,SynPal,KsPal,"KS Only Duplicated Paleologs", RDP,RDNP,RDCP,RDCNP,RDchi2,RDp]
        PeakOutput.append(Line)

        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),RG,RC,SynPal,KsPal,"KS + Syn Relaxed",RP,RNP,RCP,RCNP,Rchi2,Rp]
        PeakOutput.append(Line)

        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),SG,SC,SynPal,KsPal,"KS + Syn Strict Paleologs", SP,SNP,SCP,SCNP,Schi2,sp]
        PeakOutput.append(Line)

        #Make Lists for output files
        Z = str(thing[1])
        Line = [f"{Z}",str(max),str(min),SG,SC,SynPal,KsPal,"KS + Syn Duplicated Paleologs", SDP,SDNP,SDCP,SDCNP,SDchi2,SDp]
        PeakOutput.append(Line)


        if prior == "*":
            #####################Make a list of each paralog###########################################
            print(f"Making lists of Paleologs in primary peak (takes a while)...")
            A=[]
            B=[]
            C=[]
            D=[]
            for i in peakKS:
                for z in CDSList:
                    if i == z[1]:
                        A.append(z[0])
            for i in KSPaleologsStrict:
                for z in CDSList:
                    if i == z[1]:
                        B.append(z[0])
            for i in peak:
                for z in CDSList:
                    if i == z[1]:
                        C.append(z[0])
            for i in PaleologsStrict:
                for z in CDSList:
                    if i == z[1]:
                        D.append(z[0])

            ParalogsSummary = []
            x = 0
            for i in A:
                line = []
                line.append(i)
                if x < KSP:
                    line.append(B[x])
                if x < RP:
                    line.append(C[x])
                if x < SP:
                    line.append(D[x])
                x=x+1
                ParalogsSummary.append(line)

            #Write to a file
            if prior == "*":
                with open(f"Output/{CDS}.Paleologs.csv", "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerows(ParalogsSummary)




    #Write the Results to a file
    with open(f"Output/{CDS}.Summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(PeakOutput)
