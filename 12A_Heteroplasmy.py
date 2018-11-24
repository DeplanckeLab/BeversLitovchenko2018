#Heteroplasmy.py

#------------------------------------------------------------------------------
# STEP 12A: PYTHON SCRIPT TO DETECT HETEROPLASMY IN MITOCHONDRIAL NGS
# SEQUENCING DATA. THE SCRIPT WAS TAKEN FROM Rensch, Thomas, et al.
# "Mitochondrial heteroplasmy in vertebrates using ChIP-sequencing data."
# Genome biology 17.1 (2016): 139 and put here to ensure the reproducibility of
# the calculations.
#------------------------------------------------------------------------------

#Running this script in the folder containing pre-processed bam files will
#generate heteroplasmy files for each bam containing at least one heteroplasmic
#position.
#Calling the runAnalysis() function will then print
#the analysis results (to regenerate figures uncomment the function calls).
#Feel free to modify this script at will.

import os
import glob
import StringIO
import subprocess
import sys
import collections
import pysam as PS
import pandas as PD
import matplotlib.pyplot as PLT
import scipy.stats as ST
import numpy as NP
from matplotlib.backends.backend_pdf import PdfPages


#Main detection parameters:
#Minimum coverage
MINCOV=20
#Minimum number of bases per strand with minor allele
MINSTRAND=2
#Minimum quality of base pair
MINQUAL=23
#Minimum quality for neighbouring bases
NMINQUAL=15
#Range of neighbouring bases
NR=5
#Minimum minor allele ratio
RATIO=0.15
#Bam file mtDNA identifier (MT or chrM)
MTID="chrM"
#Use an LSF type cluster submission
CLUSTER=False
#Cluster job submission command
CLUSTERCOMMAND="bsub -o /dev/null"
#Path to Heteroplasmy.py
FILEPATH="~/Documents/MitoSeq_RPJB_ML_Aug2017/8_chrM_dm6/13_Heteroplasmy.py"

#Without arg -> checks all bams in folder for heteroplasmies
#With file as arg -> checks for heteroplasmies in file only
def main(argv):
	if len(argv)<1:
		deleteHetFilesCurrentFolder()
		files=getBams()
		for file in files:
			print file
			if(CLUSTER):
				command="%s \"python %s %s\""%(CLUSTERCOMMAND, FILEPATH, file)
				os.system(command)
			else:
				main([file])
	else:
		bamFile=argv[0]
		heteroplasmies=detectHeteroplasmies(bamFile)
		if heteroplasmies is not None:
			file=bamFile.split(".")[0]
			heteroplasmies.to_csv(file+".het", index=False)

####################################UTILITYFUNCTIONS################################
#Deletes .het files
def deleteHetFilesCurrentFolder():
	hetFiles=glob.glob("*.het")
	for h in hetFiles:
		os.remove(h)

#Returns list of bams in current folder
def getBams():
	files=glob.glob("*.bam")
	return sorted(files)

#Retrieves the length of the MT
def getMTLengthPysam(bamObj):
	return int(bamObj.lengths[bamObj.references.index(MTID)])

#Returns true if the basequality is higher than MINQUAL
def checkQualChar(char, qual):
	return (ord(char)-33)>=qual

#check neighbouring bases quality
def checkNeighbouringBases(posRR,bases,quals):
	neighbourQual=True
	for i in xrange(max(0,posRR-NR),min(len(bases)-1,posRR+NR+1),1):
		if not checkQualChar(quals[i],NMINQUAL):
			neighbourQual=False
			#print "dumping read QS"
			break	
	return neighbourQual

#returns ratio, major allele and minor allele
def getRatioMajorMinor(basesDict):
	bases=[('A','a'),('T','t'),('C','c'),('G','g')]
	basesCoverage=[]
	for tup in bases:
		X,y=tup
		cov=basesDict.get(X,0) + basesDict.get(y,0)
		basesCoverage.append((X,cov))
	basesCoverage.sort(key=lambda x: x[1], reverse=True)
	ratio=basesCoverage[1][1]/float(sum(basesDict.values()))
	return (ratio,basesCoverage[0][0],basesCoverage[1][0])

#Returns DF containing all heteroplasmies
def readHetsCSV():
	files=glob.glob("*.het")
	tabs=[]
	for file in files:
		tabs.append(PD.read_csv(file))
	df=PD.concat(tabs)
	return speciesSort(df)

#Creates a global csv file with all heteroplasmies
def writeHetsCSV():
	df=readHetsCSV()
	df.to_csv("allHets.csv",index=False)

#Creates a BED file to upload positions into third-party programs (e.g. Jalview)
def writeHetsBED():
	df=readHetsCSV()
	df["type"]="Heteroplasmy"
	df["num"]=-1
	df["col"]="red"
	df["start"]=df["Position"]-10
	df["end"]=df["Position"]+10
	bed=df[["type","Species","num","start","end", "col"]]
	bed.to_csv("allHets.bed",index=False,sep="\t")

#Returns sorted df according to species tree
def speciesSort(df):
	returnDF=df
	orderDict={"Hsap":0,"Mmul":1,"Csab":2,"Cjac":3,"Ogar":4,"Mmus":5,\
	"Rnor":6,"Hgla":7,"Ocun":8,"Btau":9,"Ddel":10,"Sscr":11,\
	"Cfam":12,"Mfur":13,"Shar":14,"Ggal":15}
	returnDF["tmp"]=returnDF['Species'].replace(orderDict)
	returnDF.sort(["tmp","Individual"],inplace=True)
	returnDF.reset_index(inplace=True, drop=True)
	returnDF=returnDF.drop("tmp",1)
	return returnDF

#Returns df with coverage data
def getCoverageDataAllBams(v=False, covThreshold=True):
	files=getBams()
	coverageInfo=[]
	for file in files:
		if v: print f
		mtLenCom="samtools view -H %s | grep 'SN:MT' | cut -f3"%(file)	
		
		mtLength=int(subprocess.check_output(mtLenCom, shell=True).split(":")[1])
		covData=getCoverageData(file)
		zeroCovBases=mtLength-len(covData)
		
		#Add the correct number of 0 coverage bases
		if zeroCovBases>0:
			v1=[MTID]*zeroCovBases
			v2=[0]*zeroCovBases
			zeroCovDF=PD.DataFrame(zip(v1,v2,v2), columns=["mt","pos","cov"])
			covData=PD.concat([covData,zeroCovDF])
		if covThreshold:
			covData=covData[covData["cov"]>=MINCOV]
		
		covRatio=(len(covData)/float(mtLength))*100
		covMean=covData["cov"].mean()
		covSD=covData["cov"].std()
		individual=file.split('.')[0]
		species=individual[0:4]
		coverageInfo.append({"Individual":individual,"Mean":covMean,\
		"SD":covSD,"Ratio":covRatio,"Species":species})
	return speciesSort(PD.DataFrame(coverageInfo))

def getCoverageData(file):
	bamObj=PS.Samfile(file, "rb")	
	mtLen=getMTLengthPysam(bamObj)
	df=PD.DataFrame(range(1,mtLen+1),columns=["pos"])
	df["cov"]=0
	df=df.set_index(["pos"], drop=False)
	depthCom="samtools depth %s"%(file)
	depth=subprocess.check_output(depthCom, shell=True)
	cd=PD.read_csv(StringIO.StringIO(depth),sep='\t', header=None, names=["mt","pos","cov"])
	cd=cd.set_index(["pos"], drop=False)
	df.update(cd)
	return df

#returns dictionary of data if position is heteroplasmic or None
def checkPosition(pos0, posReads, v=False):
	basesDict=filterQuality(pos0,posReads,v)	
	coverage=sum(basesDict.values())
	if coverage<MINCOV:
		return None
	ratio,major,minor=getRatioMajorMinor(basesDict)
	if ratio<RATIO:
		return None
	#Strand
	tab=basesDict.values()
	tab.sort(reverse=True)
	if len(tab)>=4:
		if tab[3]>=MINSTRAND:
			
			return {'Position':pos0+1,'Coverage':coverage,\
			'Ratio':ratio,'Bases':basesDict,'Major allele':major,\
			'Minor allele':minor}

#Returns a dictionary with the bases from filtered reads
def filterQuality(pos0,reads, v=False):
	basesDict=collections.defaultdict(int)
	for read in reads:
		positions=read.get_reference_positions(full_length=True)
		try: posRR=positions.index(pos0)
		except ValueError:
			if v: print "Warning skipped read at %d since position deleted in read."%(pos0)
			continue #continue loop with next read

		bases=read.seq
		quals=read.qual
		#check main base quality
		if checkQualChar(quals[posRR],MINQUAL):
			#check neighbouring bases quality
			if checkNeighbouringBases(posRR, bases, quals):
				char=bases[posRR]
				if read.is_reverse: char=char.lower()
				basesDict[char]+=1
	return basesDict
#Get colour from range of 4
def getColour(num):
	TAB=[(20,67,153),(164,25,130),(240,104,8),(4,126,52),(102,188,83),(207,225,151)]	
	r,g,b=TAB[num]
	return (r/255.,g/255.,b/255.)
#Get colour from range of 20
def getColourRange(num):
	TAB=[(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),\
	(44, 160, 44),(152, 223, 138), (214, 39, 40), (255, 152, 150),(148, 103, 189),\
	(197, 176, 213), (140, 86, 75), (196, 156, 148),(227, 119, 194),\
	(247, 182, 210), (127, 127, 127), (199, 199, 199),(188, 189, 34), (219, 219, 141),\
	(23, 190, 207), (158, 218, 229)]  	
	r,g,b=TAB[num]
	return (r/255.,g/255.,b/255.)

def printMitoMapStats(hetData):
	mitomap=getMitoMapDB()
	mitomapS=getCounts(mitomap)
	counter=collections.Counter(hetData["Major allele"]+hetData["Minor allele"])
	totalCount=sum(counter.values())
	spectrum=collections.OrderedDict()
	forward=['AT','AC','CT','AG','GT','GC']
	reverse=['TA','TG','GA','TC','CA','CG']
	for f,r in zip(forward,reverse):
		spectrum[f]=counter.get(f,0)+counter.get(r,0)
	tvtsHets=collections.OrderedDict()
	tvtsHets["TV"]=spectrum["AT"]+spectrum["AC"]+spectrum["GT"]+spectrum["GC"]
	tvtsHets["TS"]=spectrum["AG"]+spectrum["CT"]
	tvtsmito=collections.OrderedDict()
	tvtsmito["TV"]=mitomapS["AT"]+mitomapS["AC"]+mitomapS["GT"]+mitomapS["GC"]
	tvtsmito["TS"]=mitomapS["AG"]+mitomapS["CT"]
	print "TS/TV ratio in our data: %.2f"%(tvtsHets['TS']/float(tvtsHets['TV']))
	print "TS/TV ratio in MITOMAP: %.2f"%(tvtsmito['TS']/float(tvtsmito['TV']))
	dfMito=PD.DataFrame([PD.Series(tvtsHets),PD.Series(tvtsmito)])
	o2,pMito=ST.fisher_exact(dfMito)
	print "Fisher's exact test between the TS/TV datasets: %.2f"%(pMito)
	heteroplasmies=PD.Series(spectrum,name="Hets")
	mitoMapMutations=PD.Series(mitomapS,name="Mito")
	chiDF=PD.concat([mitoMapMutations,heteroplasmies], axis=1)
	chiDF['Mito107']=chiDF['Mito']/sum(chiDF['Mito'])*107
	testDF=PD.concat([chiDF['Hets'],chiDF['Mito107'].round()], axis=1)
	chi2,p,ddof,expected=ST.chi2_contingency(testDF.transpose().as_matrix())
	print "Chi-squared test between the mutational spectrums: %.2f"%(p)
	diseaseVsLi=[[5,33],[4,34]]
	o3,pDisease1=ST.fisher_exact(diseaseVsLi)
	print "Fisher's exact test between human disease positions (our data and Li 2010 data): %.2f"%(pDisease1)

	mitoDisease=getMitoMapDiseaseDB()
	diseaseOtherSpeciesVsMito=[[2,43],[len(mitoDisease),len(mitomap)]]
	#print diseaseOtherSpeciesVsMito
	#print NP.array(diseaseOtherSpeciesVsMito)
	o4,pDisease2=ST.fisher_exact(NP.array(diseaseOtherSpeciesVsMito))
	print "Fisher's exact test between non-human disease positions and mitomap: %.2f"%(pDisease2)
	


	return chiDF


def getMitoMapDB():
	a=PD.read_csv("dataFiles/PolymorphismsCoding.csv")
	b=PD.read_csv("dataFiles/PolymorphismsControl.csv")
	return PD.concat([a,b])

def getMitoMapDiseaseDB():
	a=PD.read_csv("dataFiles/MutationsCodingControlDisease.csv")
	b=PD.read_csv("dataFiles/MutationsRNADisease.csv")
	return PD.concat([a,b])

def getCounts(mitoDB):
	c=mitoDB["Nucleotide Change"].value_counts()
	d=collections.OrderedDict()
	d["AT"]=c["A-T"]+c["T-A"]
	d["AC"]=c["A-C"]+c["T-G"]
	d["CT"]=c["C-T"]+c["G-A"]
	d["AG"]=c["A-G"]+c["T-C"]
	d["GT"]=c["G-T"]+c["C-A"]
	d["GC"]=c["G-C"]+c["C-G"]
	
	return PD.Series(d, name=["Mitomap Spectrum"])

#returns the hetData DF with the additional annotations
def addGeneAnnotation(hetData):
	annotationDF=PD.read_csv("dataFiles/hetAnnotation.csv")
	if not len(annotationDF)==len(hetData):
		print "Warning: Gene annotations and hetData dont concur"
		return annotationDF
	ret=PD.merge(hetData,annotationDF)
	if not len(ret)==len(annotationDF):
		print "Warning: Gene annotations and hetData dont concur"
	if not len(ret)==len(hetData):
		print "Warning: Gene annotations and hetData dont concur"
	return ret


def printStats(hetData,covData):
	nbInds=len(covData)
	nbSpecies=len(covData["Species"].unique())
	nbHets=len(hetData)
	nbHetInds=len(hetData["Individual"].unique())
	nbHetSpecies=len(hetData["Species"].unique())
	print "Results:"
	print "Initial dataset: %d individuals from %d species"%(nbInds, nbSpecies)
	print "%d heteroplasmies detected in %d individuals (from %d species)."%(nbHets,nbHetInds, nbHetSpecies)
	print "Average coverage (SD) per heteroplasmic position: %.2f (%.2f)"%(hetData["Coverage"].mean(),hetData["Coverage"].std())
	covData["Heteroplasmies"]=covData["Individual"].map(lambda ind: len(hetData[hetData["Individual"]==ind]))
	a,b=ST.stats.pearsonr(covData["Mean"],covData["Heteroplasmies"])
	print "Correlation between number of heteroplasmies and average coverage per individual: %.2f"%(a)
	
	cd40=covData[covData["Mean"]>40]
	a2,b2=ST.stats.pearsonr(cd40["Mean"],cd40["Heteroplasmies"])
	print "Correlation between number of heteroplasmies and average coverage per individual (for individuals with average coverage >40): %.2f"%(a2)
	printMitoMapStats(hetData)	

####################################MAINFUNCTION####################################
#Returns dataframe of heteroplasmies in file
def detectHeteroplasmies(bamFile, v=False):
	bamObj=PS.Samfile(bamFile, "rb")	
	mtLen=getMTLengthPysam(bamObj)
	data=[]
	for pos in xrange(1,mtLen+1):
		pos0=pos-1
		posReads=list(bamObj.fetch(MTID,pos0,pos))	
		if len(posReads)>=MINCOV:
			posDict=checkPosition(pos0, posReads, v)
			if posDict is not None:
				data.append(posDict)
	bamObj.close()
	
	if len(data)>0:
		hetDF=PD.DataFrame(data)
		hetDF["Individual"]=bamFile.split(".")[0]
		hetDF["Species"]=bamFile[0:4]
		cols=["Species","Individual","Position",\
		"Major allele","Minor allele","Ratio","Coverage","Bases"]
		return hetDF[cols]

def runAnalysis():
	print "Running analysis..."
	cd=getCoverageDataAllBams()
	print "Coverage data uploaded."
	hd=readHetsCSV()
	print "Heteroplasmy data uploaded."
	printStats(hd,cd)

	#createCoveragePlot()
	createFigure5()
	#createSupFigureCoverageAll()
	#createSupHetCoverage()

####################################FIGUREFUNCTIONS#################################
#Creates Figure 2
def createCoveragePlot():
	covData=getCoverageDataAllBams()
	figure,axes=PLT.subplots(2,sharex=True)
	colour=0	
	for species in covData["Species"].unique():
		subDF=covData[covData["Species"] == species]
		c=getColourRange(colour)
		colour+=1
		axes[0].bar(subDF.index,subDF["Mean"],align="center",color=c, edgecolor="none")
		axes[0].errorbar(subDF.index,subDF["Mean"],yerr=subDF["SD"].values,linestyle="None", color=c)
		axes[1].bar(subDF.index,subDF["Ratio"], align="center",color=c, edgecolor="none")
	
	
	PLT.setp(axes,xticks=covData.index,xticklabels=covData["Individual"])
	#axes[0].set_title("Average read coverage per mitochondrial genome")
	axes[0].set_ylabel("Mean coverage (SD) [rds per bp]")
	axes[0].set_ylim((0,120))
	axes[0].set_xlim((-1,len(covData)))
	axes[0].spines["top"].set_visible(False)
	axes[0].spines["right"].set_visible(False)
	axes[0].tick_params(axis="x",which="both",bottom="off",top="off")
	axes[0].get_yaxis().tick_left()
	#axes[1].set_title("% bases with > 20X coverage")
	axes[1].spines["top"].set_visible(False)
	axes[1].spines["right"].set_visible(False)
	axes[1].tick_params(axis="x",which="both",bottom="off",top="off")
	axes[1].get_yaxis().tick_left()
	axes[1].set_ylabel("Coverage ratio [% of mtDNA]")
	PLT.setp(axes[0].get_xticklabels(), rotation="vertical", visible=True)
	PLT.setp(axes[1].get_xticklabels(), rotation="vertical")
	figure.tight_layout()
	PLT.show()

#Creates Figure 5
def createFigure5():
	hetData=readHetsCSV()
	covData=getCoverageDataAllBams()
	hetData=addGeneAnnotation(hetData)
	figure, axes=PLT.subplots(2,2)

	#Figure A
	bins=[0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5]
	axes[0,0].hist(hetData["Ratio"].values, bins, color=getColourRange(14))
	axes[0,0].set_xlim([0.15,0.5])
	#axes[0,0].set_yticks([0,10,20,30,40])
	axes[0,0].set_xlabel("Minor allele fraction", fontsize=18)
	axes[0,0].set_ylabel("Heteroplasmies", fontsize=18)
	axes[0,0].spines["top"].set_visible(False)
	axes[0,0].spines["right"].set_visible(False)
	axes[0,0].get_xaxis().tick_bottom()  
	axes[0,0].get_yaxis().tick_left()
	
	#Figure B
	covData["Heteroplasmies"]=covData["Individual"].map(lambda ind: len(hetData[hetData["Individual"]==ind]))
	axes[0,1].scatter(covData["Mean"],covData["Heteroplasmies"], s=50.0,edgecolor="black", color=getColourRange(14))
	axes[0,1].set_xlabel("Coverage mean [rds per bps]", fontsize=18)
	axes[0,1].set_ylabel("Number of Heteroplasmies", fontsize=18)
	axes[0,1].set_ylim([-0.5,11])
	axes[0,1].spines["top"].set_visible(False)
	axes[0,1].spines["right"].set_visible(False)
	axes[0,1].get_xaxis().tick_bottom()  
	axes[0,1].get_yaxis().tick_left()

	#Figure C
	counter=collections.Counter(hetData["Major allele"]+hetData["Minor allele"])
	totalCount=sum(counter.values())
	spectrum=collections.OrderedDict()
	forward=['AT','AC','CT','AG','GT','GC']
	reverse=['TA','TG','GA','TC','CA','CG']
	for f,r in zip(forward,reverse):
		spectrum[f]=(counter.get(f,0)+counter.get(r,0))/float(totalCount)
	plotRange=range(len(spectrum))
	axes[1,0].bar(plotRange,spectrum.values(),align="center", edgecolor="black", color=getColourRange(14))
	axes[1,0].set_xticks(plotRange)
	axes[1,0].set_xticklabels(spectrum.keys())
	axes[1,0].set_xlabel("Mutation type", fontsize=18)
	axes[1,0].set_ylabel("Fraction", fontsize=18)
	axes[1,0].set_yticks([0.1,0.2,0.3,0.4])
	axes[1,0].spines["top"].set_visible(False)
	axes[1,0].spines["right"].set_visible(False)	
	axes[1,0].get_xaxis().tick_bottom()  
	axes[1,0].get_yaxis().tick_left()

	#Figure D
	hets=addGeneAnnotation(hetData)
	counts=PD.DataFrame(hetData["Annotation"].value_counts(), columns=["Counts"])
	counts["Heteroplasmies"]=counts["Counts"].map(lambda c: c/float(sum(counts["Counts"])))
	humanMT=PD.DataFrame([1190,2511,1486,11382],index=["NC","rRNA","tRNA","Protein"], columns=["hCounts"])
	humanMT["mtDNA"]=humanMT["hCounts"].map(lambda c: c/float(sum(humanMT["hCounts"])))
	plotDF=counts.add(humanMT, fill_value=0).fillna(0).transpose()
	plotDF=plotDF.drop(["Counts","hCounts"])
	plotDF=plotDF.reindex(["mtDNA","Heteroplasmies"])
	plotRange=[0,0.65]
	nc=axes[1,1].bar(plotRange,plotDF["NC"],align="center",label="Non coding",\
	edgecolor="none",color=getColour(0),width=0.5)
	rrna=axes[1,1].bar(plotRange,plotDF["rRNA"],align="center",label="rRNA",\
	bottom=plotDF["NC"], edgecolor="none",color=getColour(1),width=0.5)
	trna=axes[1,1].bar(plotRange,plotDF["tRNA"],align="center",label="tRNA",\
	bottom=plotDF["NC"]+plotDF["rRNA"], edgecolor="none",color=getColour(2),width=0.5)
	prot=axes[1,1].bar(plotRange,plotDF["Protein"],align="center",label="Protein coding",\
	bottom=plotDF["NC"]+plotDF["rRNA"]+plotDF["tRNA"], edgecolor="none",color=getColour(3),width=0.5)
	syn=axes[1,1].bar(plotRange,plotDF["SYN"],align="center",label="Synonymous",\
	bottom=plotDF["NC"]+plotDF["rRNA"]+plotDF["tRNA"]+plotDF["Protein"],\
	edgecolor="none",color=getColour(4),width=0.5)
	nsyn=axes[1,1].bar(plotRange,plotDF["NonSYN"],align="center",label="Non synonymous",\
	bottom=plotDF["NC"]+plotDF["rRNA"]+plotDF["tRNA"]+plotDF["Protein"]+plotDF["SYN"],\
	edgecolor="none",color=getColour(5),width=0.5)
	axes[1,1].set_xticks(plotRange)
	axes[1,1].set_xticklabels(plotDF.index, fontsize=18)
	axes[1,1].spines["top"].set_visible(False)
	axes[1,1].spines["right"].set_visible(False)	
	axes[1,1].get_xaxis().tick_bottom()  
	axes[1,1].get_yaxis().tick_left()
	axes[1,1].set_xlim([-0.5,2])
	axes[1,1].set_ylabel("Fraction", fontsize=18)
	axes[1,1].legend(handles=[nsyn,syn,prot,trna,rrna,nc],loc="upper right", frameon=0)
	PLT.show()

#Figure S1 (exports it to a pdf file)
def createSupFigureCoverageAll():
	files=getBams()
	with PdfPages("figureS1.pdf") as pdfP:
		for file in files:
			covDat=getCoverageData(file)
			ind=file.split('.')[0]
			PLT.plot(covDat["pos"],covDat["cov"])
			PLT.title(ind)
			PLT.xlabel("mtDNA genome positions", fontsize=16)
			PLT.ylabel("Read coverage [reads per bp]", fontsize=16)
			PLT.fill_between(covDat["pos"],0,covDat["cov"])
			pdfP.savefig()
		PLT.close()

#Creates Figure S2
def createSupHetCoverage():
	hetData=readHetsCSV()
	bins=[20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
	hd=hetData.sort(["Coverage"])
	PLT.hist(hd["Coverage"].values, bins,color=getColourRange(14))
	#PLT.set_xlim([0.15,0.5])
	#axes[0,0].set_yticks([0,10,20,30,40])
	PLT.title("Figure S2: Read coverage of heteroplasmies", fontsize=20)
	PLT.xlabel("Heteroplasmy read coverage", fontsize=16)
	PLT.ylabel("Number of heteroplasmies", fontsize=16)
	PLT.gca().spines["top"].set_visible(False)
	PLT.gca().spines["right"].set_visible(False)
	PLT.gca().get_xaxis().tick_bottom()  
	PLT.gca().get_yaxis().tick_left()
	PLT.show()

#Creates Figure S3
#To create this figure the files need to be separated in three folders
#and coverageData with getCoverageDataAllBams() obtained separately.
def createSupFigureBoxplots(cdInput,cdProteins,cdHistones, pdf=True):
	cdInput.fillna(0,inplace=True)
	cdProteins.fillna(0,inplace=True)
	cdHistones.fillna(0,inplace=True)
	labels=["Input","Transcription factors","Histone marks"]
	bp=PLT.boxplot([cdInput["Mean"],cdProteins["Mean"],cdHistones["Mean"]], labels=labels,\
	notch=True, whis=[5,95])
	for box in bp['boxes']:
		box.set(color='black', linewidth=1)
	for whis in bp["whiskers"]:
		whis.set(color='black', linewidth=2)
	for med in bp["medians"]:
		med.set(color="black", linewidth=4)
	for flier in bp["fliers"]:
	    flier.set(marker='o', color='#e7298a', alpha=0.8)
	PLT.ylim([0,110])
	PLT.ylabel("Mean mtDNA coverage")
	PLT.title("Average coverage of different ChIP-seq files")
	if pdf:
		PLT.savefig("figureS3.pdf")
		PLT.close()
	else:
		PLT.show()


if  __name__ == "__main__":
	main(sys.argv[1:])
