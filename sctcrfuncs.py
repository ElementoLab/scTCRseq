# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 15:24:39 2015

Functions for single-cell TCRseq analysis for paired-end fastq data
@author:David Redmond (email david.redmond@outlook.com)
"""
import sys, os, commands, csv, commands, operator, tempfile, subprocess, numpy
from itertools import groupby, count
from collections import Counter
from Bio.Blast import NCBIXML

def analysis_seq_vidjil(outlabel,outdir,speciesVidjilRef,vidjildir):
    vidjil_cline=vidjildir+"vidjil-2015.10.2_x86_64 -c clones -r 1 -g "+speciesVidjilRef+" "+outdir+"/"+outlabel+"/"+outlabel+".gapfilled.final.fa -o "+outdir
    os.system(vidjil_cline)    
    
def create_gapFiller_libraries(myFastq1,myFastq2,outName,label,insertSize):
    f=open(outName+".gapfiller.libraries.txt", 'w+')
    gfFastq1=os.path.abspath(myFastq1)
    gfFastq2=os.path.abspath(myFastq2)
    print >> f, label+" bwa "+gfFastq1+" "+gfFastq2+" "+str(insertSize)+" 0.5 FR"
    f.close()
    
def run_gapFiller(label,minCov,gapFillerDir,outdir):
    myPrevDir=os.getcwd()    
    os.chdir(outdir)
    gapFiller_cline = "perl "+gapFillerDir+"GapFiller.pl -l "+label+".gapfiller.libraries.txt"+" -s "+label+".Gapped.Targets.fa"+" -m 20 -o "+str(minCov/2)+" -r 0.7 -n 5 -d 50 -t 0 -g 2 -T 1 -i 3 -b "+label
    os.system(gapFiller_cline)
    os.chdir(myPrevDir)

def print_tcr_summary_log(myGene,geneType):
    if is_empty(myGene):
        return "No "+geneType+" detected at coverage level"    
    if type(myGene) is tuple:
        return geneType+" detected:\n"+myGene[0]+"\nSequence:\n"+myGene[1][0]+"\nCoverage:\n"+str(myGene[1][1])+"\n"
    else:
        myText=""
        for i in range(0, len(myGene)):
            myText+=geneType+" detected:\n"+list(myGene)[i][0]+"\nSequence:\n"+list(myGene)[i][1][0]+"\nCoverage:\n"+str(list(myGene)[i][1][1])+"\n"
        return myText

def is_empty(any_structure):
    if any_structure:
        return False
    else:
        return True

def create_gap_fill_to_be(geneV,geneC,trimSize,gapSize):
    junctions=[]
    myGap="N"*gapSize 
    if type(geneV) is tuple:
        vFrag=geneV[1][0][:-trimSize]    
        cFrag=geneC[1][0][trimSize:]
        junctions.append((geneV[0]+"."+geneC[0],vFrag+myGap+cFrag))
    else:
        for i in range(0, len(geneV)):   
            vFrag=list(geneV)[i][1][0][:-trimSize]    
            cFrag=geneC[1][0][trimSize:]
            junctions.append((list(geneV)[i][0]+"."+geneC[0],vFrag+myGap+cFrag))
    return junctions

#choose region with longest seq then highest coverage
def choose_coverage_region(myRegion):
    result=[]
    maxLen=0
    maxCov=0    
    if type(myRegion) is tuple:
        return myRegion
    for i in range(0, len(myRegion)):
        if len(myRegion[i][1][0]) > maxLen:
            maxLen=len(myRegion[i][1][0])
            result=myRegion[i]
            maxCov=myRegion[i][1][1]
        if(len(myRegion[i][1][0])) == maxLen:
            if myRegion[i][1][1] > maxCov:
               result=myRegion[i]
               maxCov=myRegion[i][1][1] 
    return result        

#count reads in fastq files
def count_total_reads(myFastq1,myFastq2):
        result=int(commands.getoutput("zcat %s | wc -l" % myFastq1))
        result+=int(commands.getoutput("zcat %s | wc -l" % myFastq2))
        return result/4

def return_fastq_counts(myFastq1,myFastq2,outfile):
    f1=open(outfile, 'w+')
    print >>f1, count_total_reads(myFastq1,myFastq2)
    f1.close()

#return median read legth of fasta file
def return_fastq_median_read_lengths(myFastq1,outfile,lengthScript):
   cmd="zcat "+myFastq1+" | perl "+lengthScript+" - > "+outfile
   os.system(cmd)
   cmd="zcat "+myFastq1+" | perl "+lengthScript+" -"
   return(int(os.popen(cmd).read()))    

# gunzip fastq files
def gunzip_fastq(myFastq):
    command="gunzip %s" % myFastq
    os.system(command)
    
# gzip fastq files
def gzip_fastq(myFastq):
    command="gzip -1 %s" % myFastq
    os.system(command)
    
# Prepare fq files in format for blast mapping
def blast_fq_format(myFastq,outFastq):
    #command="sed '3~4d;4~4d;s/@/>/g' "+myFastq+" > "+outFastq
    command="zcat "+myFastq+" | sed '3~4d;4~4d;s/@/>/g' > "+outFastq    
    os.system(command)

#split fasta file into temporary files of 10k lines
def tempfile_split(filename, temp_dir, chunk=10**4):
    fns={}
    with open(filename, 'r') as datafile:
        groups = groupby(datafile, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=temp_dir,prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
    return fns        

#blast fastqs against ref tcr databases
def blastall_v_regions(myFastq1,myFastq2,myRef,outputfile,eVal,blastallDir):
    fns={}
    chunk=10**4
    with open(myFastq1, 'r') as datafile1:
        groups = groupby(datafile1, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=tempfile.mkdtemp(),prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
            blastn_cline = blastallDir+"blastall -p blastn -o "+str(outfile.name)+".blast.out -i "+str(outfile.name)+" -d "+myRef+" -e "+str(eVal)+" -m 8 -b 1"    
            os.system(blastn_cline+" > /dev/null 2>&1")
            os.system("cat "+str(outfile.name)+".blast.out >> "+outputfile)
            os.remove(str(outfile.name)+".blast.out")
            os.remove(str(outfile.name))
            testvar=commands.getstatusoutput("dirname "+str(outfile.name))
            os.system("rm -r "+testvar[1])
    fns={}
    with open(myFastq2, 'r') as datafile2:
        groups = groupby(datafile2, key=lambda k, line=count(): next(line) // chunk)
        for k, group in groups:
            with tempfile.NamedTemporaryFile(delete=False,
                           dir=tempfile.mkdtemp(),prefix='{}_'.format(str(k))) as outfile:
                outfile.write(''.join(group))
                fns[k]=outfile.name   
            blastn_cline = blastallDir+"blastall -p blastn -o "+str(outfile.name)+".blast.out -i "+str(outfile.name)+" -d "+myRef+" -e "+str(eVal)+" -m 8 -b 1"    
            os.system(blastn_cline+" > /dev/null 2>&1")
            os.system("cat "+str(outfile.name)+".blast.out >> "+outputfile)
            os.remove(str(outfile.name)+".blast.out")
            os.remove(str(outfile.name))
            testvar=commands.getstatusoutput("dirname "+str(outfile.name))
            os.system("rm -r "+testvar[1])

def listToStringWithoutBrackets(list1):
    return str(list1).replace('[','').replace(']','')

def run_seqtk(inputList,inputFastq1,inputFastq2,outputFastq,seqTkDir):
    command1=seqTkDir+"seqtk subseq "+inputFastq1+" "+inputList+" >> "+outputFastq
    command2=seqTkDir+"seqtk subseq "+inputFastq2+" "+inputList+" >> "+outputFastq
    os.system(command1)
    os.system(command2)

#main ftn to return V or C region counts and alignments for gapfiller processing
def return_counts_and_alignment(blastHitsFile,outName,fastq1,fastq2,seqTkDir,threshold,minBlastAlignedLength):
    myHits=[]
    with open(blastHitsFile) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
                if(int(row[3])>=minBlastAlignedLength):
                        myHits.append(row[1])
    gene_table={}
    gene_table=Counter(myHits)
    perc_table={}
    for gene in gene_table:
        perc_table[gene]=float(gene_table[gene])/float(sum(gene_table.values()))
    gene_table_output = { k: [ gene_table[k], perc_table[k] ] for k in gene_table }
    sorted_gto = sorted(gene_table_output.items(), key=operator.itemgetter(1),reverse=True)
    f1=open(outName+".counts.txt", 'w+')
    for item in sorted_gto:
        print >>f1, item[0],",",listToStringWithoutBrackets(item[1])
    f1.close() 
    with open(blastHitsFile) as f:    
       reader = csv.reader(f, delimiter="\t")
       d = list(reader)   
    myReads = [myReads[0:2] for myReads in d]
    f2=open(blastHitsFile+".candidates.txt", 'w+')    
    for i in range(0,len(sorted_gto)):
        if sorted_gto[i][1][1] > threshold:
            for j in range(0,len(myReads)):
                if myReads[j][1]==sorted_gto[i][0]:
                        print >> f2, myReads[j][0]
        else:
            next
    f2.close()
    run_seqtk(blastHitsFile+".candidates.txt",fastq1,fastq2,outName+".local.cands.fq",seqTkDir)


#blastall the major V and return XML ouput
def blast_single_v_region(myFastq,myRef,outfile,eVal,blastallDir):
    tempFasta = tempfile_split(myFastq, tempfile.mkdtemp(), chunk=10**4)
    for key, value in tempFasta.iteritems():
        blastn_cline = blastallDir+"blastall -p blastn -o "+value+".blast.out -i "+value+" -d "+myRef+" -e "+str(eVal)+" -m 7 -b 1"
        os.system(blastn_cline+" > /dev/null 2>&1")
        os.system("cat "+value+".blast.out >> "+outfile)
        os.remove(value+".blast.out")
        os.remove(value)

def perform_targeted_alignment(candidateFile, outAlignment, candidateGene):
    result=open(candidateFile,"r")
    f=open(outAlignment, 'w+')     
    records=NCBIXML.parse(result)
    for item in records:
        for alignment in item.alignments:
            if alignment.accession == candidateGene:
                for hsp in alignment.hsps:
                    myAlignment="_"*(hsp.sbjct_start-1)
                    myAlignment=myAlignment+hsp.query[0:190]
                    myAlignment=myAlignment+"_"*(alignment.length-len(myAlignment))
                    print >> f, (myAlignment)
    f.close()
    result.close()

#create pileup of reads mapping to particular V or C
def return_consensus(alignmentFile,minCov):
    with open(alignmentFile,"rt") as infile:
        matrix = [list(line.strip()) for line in infile.readlines()]
    transpose=[list(x) for x in zip(*matrix)]
    vGeneA=[]
    vGeneC=[]
    vGeneG=[]
    vGeneT=[]
    for i in range(0,len(transpose)):
        vGeneA.append(transpose[i].count("A"))
        vGeneC.append(transpose[i].count("C"))
        vGeneG.append(transpose[i].count("G"))
        vGeneT.append(transpose[i].count("T"))
    vGeneCoverage=[]
    vGeneCoverage=[a + b + c + d for a, b, c, d in zip(vGeneA, vGeneC, vGeneG, vGeneT)]
    vGeneAlignment={}
    vGeneAlignment={"A":vGeneA,"C":vGeneC,"G":vGeneG,"T":vGeneT}
    vGeneConsensus=[]
    vGeneConsensusCoverage=[]
    for i in range(0,len(vGeneAlignment["A"])):
        consensusCount=0    
        consensusBase="N"
	for base in vGeneAlignment:
            if(vGeneAlignment[base][i] > consensusCount):
                consensusCount = vGeneAlignment[base][i] 
                consensusBase = base
        vGeneConsensus.append(consensusBase)
        vGeneConsensusCoverage.append(consensusCount)
    oldsubseq = []
    newsubseq = []
    for i in range(0,len(vGeneConsensus)):
        if (vGeneConsensusCoverage[i] > minCov):
            newsubseq.append(vGeneConsensus[i])
        else:
            if (len(newsubseq) > len(oldsubseq)):
                oldsubseq = newsubseq
                newsubseq = []
    #vGeneTarget=max(oldsubseq,newsubseq)
    if len(oldsubseq) >= len(newsubseq):
        vGeneTarget=oldsubseq
    else:
        vGeneTarget=newsubseq    
    return ("".join(vGeneTarget),numpy.mean(vGeneConsensusCoverage))
   
def compare_alignments(candGene):
    if len(candGene) == 2:    
        return(compare_agg(candGene))
    else:
        temp1=map(candGene.__getitem__, (0,1))
        temp2=map(candGene.__getitem__, (0,2))
        temp3=map(candGene.__getitem__, (1,2))
        temp1=compare_agg(temp1)
        temp2=compare_agg(temp2)
        temp3=compare_agg(temp3)
        ###ADDED IN
        if type(temp1) is tuple:
            temp1=list("")
        if type(temp2) is tuple:
            temp2=list("")
        if type(temp3) is tuple:
            temp3=list("")
        #END OF ADDIN    
        return(set(temp1+temp2+temp3))
                
def compare_agg(candGene):
    if candGene[1][1][0] in candGene[0][1][0]:
        return candGene[0]
    elif candGene[0][1][0] in candGene[1][1][0]:
        return candGene[1]
    else:
        return candGene    
        
def get_variable_regions(myFastq1,myFastq2,myRef,outName,minCov,eVal,blastallDir,seqTkDir,threshold,minBlastAlignedLength):
    blastall_v_regions(myFastq1,myFastq2,myRef,outName+".matches.txt",eVal,blastallDir)
    return_counts_and_alignment(outName+".matches.txt",outName,myFastq1,myFastq2,seqTkDir,threshold,minBlastAlignedLength)
    blast_single_v_region(outName+".local.cands.fq",myRef,outName+".matches.xml",eVal,blastallDir)
    candGene=[]
    with open(outName+".counts.txt", 'r') as csvfile:
        reader = csv.reader(csvfile)
        table = [[e for e in r] for r in reader]
    for i in range(0, len(table)):
        if float(table[i][2]) > threshold:
            myCand=table[i][0]           
            candGene.append(myCand.strip())
    for i in range(0, len(candGene)):
        perform_targeted_alignment(outName+".matches.xml", outName+"."+str(candGene[i]).replace("/","")+".aln", candGene[i])
        try:
            candGene[i]=(candGene[i],return_consensus(outName+"."+str(candGene[i]).replace("/","")+".aln",minCov))
        except Exception:
            pass 
    
    for i in candGene[:]:
        if len(i) != 2:
            candGene.remove(i)
    for i in candGene:
        if numpy.isnan(float(i[1][1])):
            candGene.remove(i)
    if len(candGene) > 1:
        candGene=compare_alignments(candGene)
    return candGene
        
            
