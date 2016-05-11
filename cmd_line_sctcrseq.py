"""
Created on Fri Sep  4 12:36:46 2015
single-cell TCRseq analysis for paired-end fastq data
@author:David Redmond (email david.redmond@outlook.com)
"""
import sys, getopt, os, commands, csv, commands, operator, tempfile, subprocess, numpy
from optparse import OptionParser
from itertools import groupby, count
from collections import Counter
from Bio.Blast import NCBIXML
import sys
#insert path of scTCRseq program here - NEEDS TO BE MANUALLY CHANGED
sys.path.insert(0, '/path/to/sctcrseq/')
import sctcrfuncs

#directories for programs - NEEDS TO BE MANUALLY CHANGED
seqTkDir="/path/to/seqtk/"
blastallDir="/path/to/blastall/"
gapFillerDir="/path/to/GapFiller_v1-10_linux-x86_64/"
vidjildir="/path/to/vidjil/"
lengthScript="/path/to/calc.median.read.length.pl"

#location for FASTA BLAST reference sequences downloadable from imgt.org - NEEDS TO BE MANUALLY CHANGED
humanTRAVblast="/path/to/TRAV.human.fa"
humanTRBVblast="/path/to/TRBV.human.fa"
humanTRACblast="/path/to/TRAC.human.fa"
humanTRBCblast="/path/to/TRBC.human.fa"
mouseTRAVblast="/path/to/TRAV.mouse.fa"
mouseTRBVblast="/path/to/TRBV.mouse.fa"
mouseTRACblast="/path/to/TRAC.mouse.fa"
mouseTRBCblast="/path/to/TRBC.mouse.fa"

#location for Vidjil BLAST reference sequences in vidjil program -  NEEDS TO BE MANUALLY CHANGED
humanVidjilRef="/path/to/tr_germline/human"
mouseVidjilRef="/path/to/tr_germline/mouse"


# Input and check variables from command line
# Read command line args

parser = OptionParser()
usage = "usage: %prog [options] --fastq1 FASTQ1 --fastq2 FASTQ2 --species human/mouse --outdir OUTPUT DIRECTORY --label OUTPUT LABEL"
parser = OptionParser(usage=usage)
parser.add_option("--fastq1", dest="myFastq1",action="store",type="string",
                  help="enter R1 reads as FASTQ1 in fastq.gz format", metavar="FASTQ1")
parser.add_option("--fastq2", dest="myFastq2",action="store",type="string",
                  help="enter R1 reads as FASTQ2 in fastq.gz format", metavar="FASTQ2")                  
parser.add_option("-s","--species", dest="species",action="store",type="string",
                  help="enter SPECIES either human or mouse (deault human)",default="human", metavar="SPECIES")
parser.add_option("-e","--eval", dest="eVal",action="store",type="float",default=1e-10,
                  help="enter BLAST E-VALUE threshold (default 10e-10)", metavar="E-VALUE")
parser.add_option("-c","--cov", dest="minCov",action="store",type="float",default=5,
                  help="enter minimum coverage (default 5)", metavar="MIN COVERAGE")
parser.add_option("-i","--insertsize", dest="insertSize",action="store",type="float",default=300,
                  help="enter insert size (default 300)", metavar="INSERT SIZE")
parser.add_option("-o","--outdir", dest="outdir",action="store",type="string",
                  help="enter OUTDIR of output if required",default="", metavar="OUTDIR")
parser.add_option("-l","--label", dest="outlabel",action="store",type="string",
                  help="enter LABEL of output", metavar="LABEL")                  
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

#standard error handling in cmd line
if not (str(options.myFastq1).endswith(".fastq.gz") or str(options.myFastq1).endswith(".fq.gz")):
    parser.error("--fastq1 must be in .fastq.gz or .fq.gz format")
if not (str(options.myFastq2).endswith(".fastq.gz") or str(options.myFastq2).endswith(".fq.gz")):
    parser.error("--fastq2 must be in .fastq.gz or .fq.gz format")
if not (str(options.species)=="human" or str(options.species)=="mouse"):
    parser.error("must select human or mouse for --species")
if str(options.outlabel)=="None":
    parser.error("must select output label for --label")
#print parameters used    
print "running on PE FQ files %s and %s" % (options.myFastq1, options.myFastq2)
print "species=%s" % options.species
print "BLAST eval=%.3g" % options.eVal
print "MIN COVERAGE=%.3g" % options.minCov
print "INSERT SIZE=%.3g" % options.insertSize
print "outdir=%s" % options.outdir
print "label=%s" %options.outlabel
#set ref databases
if options.species=="human":
    blastRef=(humanTRAVblast, humanTRACblast, humanTRBVblast, humanTRBCblast)
    speciesVidjilRef=humanVidjilRef
else:
    blastRef=(mouseTRAVblast, mouseTRACblast, mouseTRBVblast, mouseTRBCblast)
    speciesVidjilRef=mouseVidjilRef

myFastq1=options.myFastq1
myFastq2=options.myFastq2
if sctcrfuncs.is_empty(options.outdir):
    outName=options.outlabel
else:
    outName=options.outdir+"/"+options.outlabel
minCov=options.minCov
eVal=options.eVal
insertSize=options.insertSize
species=options.species
#threshold to consider V regions
threshold=0.15

#gunzip fq files, format and rezip them
sctcrfuncs.blast_fq_format(myFastq1,outName+".formatted.1.fq")
sctcrfuncs.blast_fq_format(myFastq2,outName+".formatted.2.fq")

#count reads
sctcrfuncs.return_fastq_counts(myFastq1,myFastq2,outName+".readcounts.txt")
#get median read lengths
medianLengths=sctcrfuncs.return_fastq_median_read_lengths(myFastq1,outName+".medianreadlength.txt",lengthScript)
minBlastAlignedLength=max(medianLengths/3,20)

#get TCR V and C genes
alphaVgene=sctcrfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[0],outName+".alpha.V",minCov,eVal,blastallDir,seqTkDir,threshold,minBlastAlignedLength)
alphaCgene=sctcrfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[1],outName+".alpha.C",minCov,eVal,blastallDir,seqTkDir,threshold,minBlastAlignedLength)
betaVgene=sctcrfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[2],outName+".beta.V",minCov,eVal,blastallDir,seqTkDir,threshold,minBlastAlignedLength)
betaCgene=sctcrfuncs.get_variable_regions(outName+".formatted.1.fq",outName+".formatted.2.fq",blastRef[3],outName+".beta.C",minCov,eVal,blastallDir,seqTkDir,threshold,minBlastAlignedLength)

#choose constant regions
alphaCgene=sctcrfuncs.choose_coverage_region(alphaCgene)
betaCgene=sctcrfuncs.choose_coverage_region(betaCgene)

#write to output files and build receptors to be gap-filled
gapped_junctions=[]
if not (sctcrfuncs.is_empty(alphaVgene) or sctcrfuncs.is_empty(alphaCgene)):
    gapped_junctions.append(sctcrfuncs.create_gap_fill_to_be(alphaVgene,alphaCgene,3,66))
if not (sctcrfuncs.is_empty(betaVgene) or sctcrfuncs.is_empty(betaCgene)):
    gapped_junctions.append(sctcrfuncs.create_gap_fill_to_be(betaVgene,betaCgene,3,66))

#gapped targets
f=open(outName+".Gapped.Targets.fa", 'w+')    
for i in range(0, len(gapped_junctions)):
    if sctcrfuncs.is_empty(gapped_junctions[i]):
        next
    else:
        for j in range(0,len(gapped_junctions[i])):
          print >> f, ">"+gapped_junctions[i][j][0]
          print >> f, gapped_junctions[i][j][1]
f.close()

#print summary of targets to output.log
f=open(outName+".Output.log", 'w+')
print >> f, sctcrfuncs.print_tcr_summary_log(alphaVgene,"TCR Alpha Variable Region")
print >> f, sctcrfuncs.print_tcr_summary_log(alphaCgene,"TCR Alpha Constant Region")
print >> f, sctcrfuncs.print_tcr_summary_log(betaVgene,"TCR Beta Variable Region")
print >> f, sctcrfuncs.print_tcr_summary_log(betaCgene,"TCR Beta Constant Region")
if ((sctcrfuncs.is_empty(alphaVgene) or sctcrfuncs.is_empty(alphaCgene)) and (sctcrfuncs.is_empty(betaVgene) or sctcrfuncs.is_empty(betaCgene))):
    print >> f, "Not sufficient information to run gapFiller on either alpha or beta chains. Process exited."
    f.close()    
    print "Not sufficient information to run gapFiller on either alpha or beta chains. Process exited."
    sys.exit()
f.close()

#create gapfiller library files
sctcrfuncs.create_gapFiller_libraries(myFastq1,myFastq2,outName,options.outlabel,insertSize)
#run gapfiller
sctcrfuncs.run_gapFiller(options.outlabel,minCov,gapFillerDir,options.outdir)
#analysis output joined tcr regions in vidjil to process junction names quickly
sctcrfuncs.analysis_seq_vidjil(options.outlabel,options.outdir,speciesVidjilRef,vidjildir)





