#
## this program does deep seq data analysis
#
### import necessary libraries
import re,os,sys,logging,time,datetime,pysam;
from distutils.version import LooseVersion
#
### checking out the python version
if sys.version_info < (2,6) or sys.version_info >= (3,0):
  print ("Python Version error: must use phthon 2.6 or greater (not supporting python 3 yet)");
  sys.exit(-1);
#else:
#  print "Passed version test";
#  sys.exit(0);
import commands;
#
###MATS version
#
MATS_ver="3.2.5" 
#
############################ parameter variables
### seven required values
gtf=''; 	## gtf file 
bIndex='';	## bowtie index
outDir='';	## output directory
readType = '';   ## read type, single or paired
readLength = 0; ## read length
analysis='U'; ## P for paired analysis, U for unpaired analysis. default is U
libType='fr-unstranded'; ## it coulbe fr-firststrand or fr-secondstrand
novelSS=0; ## do not find novel splice sites by default
#
### input is either fastq or bam
## fastq
s1=''; 		## sample_1  
s2=''; 		## sample_2 
## bam
b1='';	## bam_1
b2='';  ## bam_2
#
### with default values
insertLength=[];
sigma = [];
insertLength2=[];
sigma2 = [];
anchorLength=1;
tophatAnchor=6; ## default
c=0.0001; 
bamFile=0; ## by default, no bam file
MacOS = 0; ## by default, not a MacOS
if os.uname()[0]=="Darwin": ## it is Mac OS
  MacOS=1; ## MacOS is true
keepTemp = 0; ## by default, do not keep temp
lite=0; ## by default, do not do the lite

### checking out the argument names
#validArgList=['-s1','-b1','-s2','-b2','-gtf','-bi','-o','-t','-len','-a','-r1','-r2','-sd1','-sd2','-c','-analysis','-expressionChange','-keepTemp'];
validArgList=['-s1','-b1','-s2','-b2','-gtf','-bi','-o','-t','-len','-a','-r1','-r2','-sd1','-sd2','-c','-analysis','-keepTemp','-libType','-lite', '-novelSS'];
for argIndex in range(1,len(sys.argv)): ## going through the all parameters 
  if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
    print ('Not valid argument: %s' % sys.argv[argIndex]);
    print ('Please provide valid arguments.');
    sys.exit();

for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
  if(sys.argv[paramIndex] == '-s1' or sys.argv[paramIndex] == '-b1'):  ## sample_1
    if (sys.argv[paramIndex] == '-b1'): ## bam file here
      bamFile=1;
      bIndex='dummy'; 
    paramIndex += 1;  ## increase index
    s1 = sys.argv[paramIndex];    
  elif (sys.argv[paramIndex] == '-s2' or sys.argv[paramIndex] == '-b2'):  ## sample_2
    if (sys.argv[paramIndex] == '-b2'): ## bam file here
      bamFile=1;
      bIndex='dummy'; 
    paramIndex += 1;  ## increase index
    s2 = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-gtf'):  ## gtf file
    paramIndex += 1;  ## increase index
    gtf = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-bi'):  ## bowtie index base
    paramIndex += 1;  ## increase index
    bIndex = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-o'):  ## output folder
    paramIndex += 1;  ## increase index
    outDir = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-t'):  ## readtype, single or paired
    paramIndex += 1;  ## increase index
    readType = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-len'):  ## read length
    paramIndex += 1;  ## increase index
    readLength = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-r1'):  ## insert length average, sample 1
    paramIndex += 1;  ## increase index
    insertLength = [int(float(kk)) for kk in sys.argv[paramIndex].split(',')];
  elif (sys.argv[paramIndex] == '-r2'):  ## insert length average, sample 1
    paramIndex += 1;  ## increase index
    insertLength2 = [int(float(kk)) for kk in sys.argv[paramIndex].split(',')];
  elif (sys.argv[paramIndex] == '-sd1'):  ## standard deviation, sample 1
    paramIndex += 1;  ## increase index
    sigma = [int(float(kk)) for kk in sys.argv[paramIndex].split(',')];
  elif (sys.argv[paramIndex] == '-sd2'):  ## standard deviation, sample 2
    paramIndex += 1;  ## increase index
    sigma2 = [int(float(kk)) for kk in sys.argv[paramIndex].split(',')];
  elif (sys.argv[paramIndex] == '-c'):  ## MATS c option, deltaPSI
    paramIndex += 1;  ## increase index
    c = float(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-a'):  ## anchor length for aligner, not junction reads
    paramIndex += 1;  ## increase index
    tophatAnchor = int(sys.argv[paramIndex]);
  elif (sys.argv[paramIndex] == '-analysis'):  ## MATS t option, paired data or not, it's either U or P
    paramIndex += 1;  ## increase index
    analysis = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-libType'):  ## library type, strand-specific or not
    paramIndex += 1;  ## increase index
    libType = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-keepTemp'):  ## keep temp files, no value needed
    keepTemp = 1;  ## keep temp files
  elif (sys.argv[paramIndex] == '-lite'):  ## keep temp files, no value needed
    lite = 1;  ## keep temp files
  elif (sys.argv[paramIndex] == '-novelSS'):  ## find novel splice sites (0 no, 1 yes)
    paramIndex += 1;  ## increase index
    novelSS = int(sys.argv[paramIndex]);

#  else: ### not valid param.. exit
#    print("Not a valid param detected: %s" % sys.argv[paramIndex]);
#    sys.exit();

### checking out the required arguments
if lite==0:
  if (s1=='' or  s2=='' or  gtf=='' or bIndex=='' or outDir=='' or readLength==0 or readType==''): ### at least one required param is missing
    print ('Not enough arguments!!');
    print ('Usage (with fastq files):\n\tpython RNASeq-MATS.py -s1 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -s2 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -gtf gtfFile -bi bowtieIndexBase -o outDir -t readType -len readLength [options]*');
    print ('Example\n\tpython RNASeq-MATS.py -s1 sample_1.rep_1.R1.fastq:sample_1.rep_1.R2.fastq,sample_1.rep_2.R1.fastq:sample_1.rep_2.R2.fastq -s2 sample_2.rep_1.R1.fastq:sample_2.rep_1.R2.fastq,sample_2.rep_2.R1.fastq:sample_2.rep_2.R2.fastq -gtf gtf/Homo_sapiens.Ensembl.GRCh37.65.gtf -bi ~/bowtieIndexes/hg19 -o out_test -t paired -len 50 -a 8 -r1 72,75 -sd1 40,35 -r2 70,65 -sd2 48,45 8 -c 0.05 -analysis U \n');
    print ('Usage (with bam files):\n\tpython RNASeq-MATS.py -b1 s1_rep1.bam[,s1_rep2.bam]* -b2 s2.rep1.bam[,s2.rep2.bam]* -gtf gtfFile -o outDir -t readType -len readLength [options]*');
    print ('Example\n\tpython RNASeq-MATS.py  -b1 ESRP_1.bam,ESRP_2.bam -b2 EV_1.bam,EV_2.bam -gtf gtf/Homo_sapiens.Ensembl.GRCh37.65.gtf -o out_test -t paired -len 50 -a 8 -r1 72,75 -sd1 40,35 -r2 70,65 -sd2 48,45 8 -c 0.05 -analysis P');
    sys.exit();
elif lite==1:
  if ((s1=='' and  s2=='') or gtf=='' or bIndex=='' or outDir=='' or readLength==0 or readType==''): ### at least one required param is missing
    print ('Not enough arguments!!');
    print ('Usage (with fastq files):\n\tpython RNASeq-MATS.py -s1 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -s2 rep1_1[:rep1_2][,rep2_1[:rep2_2]]* -gtf gtfFile -bi bowtieIndexBase -o outDir -t readType -len readLength [options]*');
    print ('Example\n\tpython RNASeq-MATS.py -s1 sample_1.rep_1.R1.fastq:sample_1.rep_1.R2.fastq,sample_1.rep_2.R1.fastq:sample_1.rep_2.R2.fastq -s2 sample_2.rep_1.R1.fastq:sample_2.rep_1.R2.fastq,sample_2.rep_2.R1.fastq:sample_2.rep_2.R2.fastq -gtf gtf/Homo_sapiens.Ensembl.GRCh37.65.gtf -bi ~/bowtieIndexes/hg19 -o out_test -t paired -len 50 -a 8 -r1 72,75 -sd1 40,35 -r2 70,65 -sd2 48,45 8 -c 0.05 -analysis U\n');
    print ('Usage (with bam files):\n\tpython RNASeq-MATS.py -b1 s1_rep1.bam[,s1_rep2.bam]* -b2 s2.rep1.bam[,s2.rep2.bam]* -gtf gtfFile -o outDir -t readType -len readLength [options]*');
    print ('Example\n\tpython RNASeq-MATS.py  -b1 ESRP_1.bam,ESRP_2.bam -b2 EV_1.bam,EV_2.bam -gtf gtf/Homo_sapiens.Ensembl.GRCh37.65.gtf -o out_test -t paired -len 50 -a 8 -r1 72,75 -sd1 40,35 -r2 70,65 -sd2 48,45 8 -c 0.05 -analysis U');
    sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

### process input params ####
#
### fastq files or bam file
#
sample_1=s1.split(','); ## each end of a pair is separated by :
sample_2=s2.split(','); ## each end of a pair is separated by :
if lite==1:
  sample_2=[];
#
SEPE = 'SE'; ## single-end or paired
#
if analysis=='P': ## check if we got the same number of replicates here
  if (len(sample_1)!=len(sample_2)) or len(sample_1)<3: ## different number of replicates per sample.. wrong!!
    print("PAIRED analysis requires the same number of replicates per sample.\nThe number of replicate must be greater than 2 for paired analysis.");
    sys.exit();
#
if readType=="paired": ## single-end
  SEPE='PE';

if lite==0:
  if bamFile==1 and ( ((sample_1[0].split('.'))[-1].strip()).upper() !='BAM' or ((sample_2[0].split('.'))[-1].strip()).upper() !='BAM'):
    print "Incorrect file type. Need to provide bam file for -b1 and -b2";
    sys.exit();
elif lite==1:
  if bamFile==1 and ( ((sample_1[0].split('.'))[-1].strip()).upper() !='BAM'):
    print "Incorrect file type. Need to provide bam file for -b1";
    sys.exit();

##### checking GTF format #####
#
tempGTF_file = open(gtf); ## open gtf file
for line in tempGTF_file:
  if line.strip()[0]=='#': ## comments, skip this line
    continue;
  gtfEle = line.strip().split('\t');
  if len(gtfEle)<9: ## may be incorrect gtf format
    print "Incorrect GTF file format. GTF file must have 9 tab-delimited columns.";
    sys.exit();
  break;  ## just check the first non-comment column


####### checking readLength ##########
if bamFile==0: ## fastq file was provided
  for fq in sample_1: ## examine each file
    fpfp=[];
    if readType=="paired":
      fpfp = fq.split(':');
    else:
      fpfp.append(fq);
    for yyy in fpfp: ## for each fastq file
      fff=open(yyy);
      line1=fff.readline();
      line2=fff.readline().strip();
      if len(line2)!=readLength: ### different readLength
        print "Incorrect readLength. %s has a read length of %d, while readLength param is %d" % (yyy,len(line2),readLength);
        fff.close();
        sys.exit();
      fff.close();
  if lite==0:
    for fq in sample_2: ## examine each file
      fpfp=[];
      if readType=="paired":
        fpfp = fq.split(':');
      else:
        fpfp.append(fq);
      for yyy in fpfp: ## for each fastq file
        fff=open(yyy);
        line1=fff.readline();
        line2=fff.readline().strip();
        if len(line2)!=readLength: ### different readLength
          print "Incorrect readLength. %s has a read length of %d, while readLength param is %d" % (yyy,len(line2),readLength);
          fff.close();
          sys.exit();
        fff.close();

else: ## bam file was provided
  for fq in sample_1: ## examine each file

    myCmd = 'samtools view '+fq+' | head -n 1';
    status,output=commands.getstatusoutput(myCmd);
    myLen=len(output.strip().split('\t')[9]);
    if myLen !=readLength: ### different readLength
      print "Incorrect readLength. %s has a read length of %d, while readLength param is %d" % (fq,myLen,readLength);
      sys.exit();

## getting readLength and junctionLength here
#
junctionLength = 2*(readLength-anchorLength);

### assigning default insertLength and sigma

if len(insertLength)==0: ## not provide, assign the default
  insertLength = [100 for kk in sample_1];
if len(sigma)==0: ## assign default
  sigma = [300 for kk in sample_1];
if len(insertLength2)==0: ## not provide, assign the default
  insertLength2 = [100 for kk in sample_2];
if len(sigma2)==0: ## assign default
  sigma2 = [300 for kk in sample_2];

## check if the insertLength and sigma values are assigned correctly (correct number of values)
if len(insertLength)!=len(sample_1) or len(sigma) != len(sample_1): ## not correct number of values
  print("The number of -r1 and -sd1 values must match the number of replicates for sample_1.");    
  print("Please check -r1 and -sd1 parameters");    
  sys.exit();

if len(insertLength2)!=len(sample_2) or len(sigma2) != len(sample_2): ## not correct number of values
  print("The number of -r2 and -sd2 values must match the number of replicates for sample_2.");    
  print("Please check -r2 and -sd2 parameters");    
  sys.exit();

if libType not in ['fr-unstranded','fr-firststrand','fr-secondstrand'] : ## strand specific information is not correct
  print("The libType must be fr-unstranded,fr-firststrand, or fr-secondstrand");
  print libType;
  sys.exit();

os.system('mkdir -p '+ outDir);
oFile = open(outDir+'/commands.txt', 'a'); ## file that will contain list of commands excuted here

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir+'/log.RNASeq-MATS.'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w');

##### Getting Start Time ######
logging.debug('rMATS version: %s' % MATS_ver);
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

scriptPath = os.path.abspath(os.path.dirname(__file__));  ## absolute script path
binPath = scriptPath + '/bin';  ## absolute bin path
outPath = os.path.abspath(outDir); ## absolute output path

s1Path = outPath + '/SAMPLE_1';
os.system('mkdir -p '+ s1Path);
s2Path = outPath + '/SAMPLE_2';
os.system('mkdir -p '+ s2Path);

## making folders for replicates ##
s1rPath = s1Path+'/REP_';
s2rPath = s2Path+'/REP_';
for rr in range(0,len(sample_1)): ## sample_1
  os.system('mkdir -p '+ s1rPath+str(rr+1));
for rr in range(0,len(sample_2)): ## sample_2
  os.system('mkdir -p '+ s2rPath+str(rr+1));

finalPath = outPath+'/MATS_output';
os.system('mkdir -p '+ finalPath);

tempPath = outPath + '/temp';
os.system('mkdir -p '+ tempPath);

#geneExprFile = tempPath+"/geneExprOutput.txt";

#
### putting keys in log file
#
logging.debug("################### folder names and associated input files #############");
for fki in range(0,len(sample_1)): ## for each replicate of sample_1
  repTempFolder = "SAMPLE_1\REP_"+str(fki+1);
  associatedFile = sample_1[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);

for fki in range(0,len(sample_2)): ## for each replicate of sample_2
  repTempFolder = "SAMPLE_2\REP_"+str(fki+1);
  associatedFile = sample_2[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);
  
logging.debug("#########################################################################\n");

########## functions here... ############

def doSTARMapping(): ## do STAR mapping
  logging.debug("mapping the first sample");

  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
    cmd += '--alignSJDBoverhangMin ' + str(tophatAnchor) + ' --alignIntronMax 300000 --genomeDir '+bIndex+ ' --sjdbGTFfile ' + gtf; 
    cmd += ' --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn ';
    if SEPE=='PE': ## paired-end
      cmd += sample_1[rr].split(':')[0]+' '+sample_1[rr].split(':')[1];
    else: ## single-end
      cmd += sample_1[rr];
    oFile.write('######  running STAR for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("mapping sample_1, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_1, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  if lite==1:
    return;

  logging.debug("mapping the second sample");
  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
    cmd += '--alignSJDBoverhangMin ' + str(tophatAnchor) + ' --alignIntronMax 300000 --genomeDir '+bIndex+ ' --sjdbGTFfile ' + gtf; 
    cmd += ' --outFileNamePrefix ' + rTempFolder + '/ --readFilesIn ';
    if SEPE=='PE': ## paired-end
      cmd += sample_2[rr].split(':')[0]+' '+sample_2[rr].split(':')[1];
    else: ## single-end
      cmd += sample_2[rr];
    oFile.write('######  running STAR for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("mapping sample_2, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_2, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  return;
##### end of doSTARMapping ####


def doTophatMapping(): ## do tophat mapping
  logging.debug("mapping the first sample");

  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    cmd = 'tophat -a '+str(tophatAnchor)+' -m 0 -I 300000 -p 4 -g 20 --library-type ' + libType + ' --no-novel-indels ';
    cmd += ' -N 3 --segment-mismatches 2 -G '+gtf+' -o '+rTempFolder;
    if SEPE=='PE': ## paired-end
      cmd += ' -r '+str(insertLength[rr]) + ' --mate-std-dev ' + str(sigma[rr])+' ' + bIndex +' '+sample_1[rr].split(':')[0]+' '+sample_1[rr].split(':')[1];
    else: ## single-end
      cmd += ' ' + bIndex +' '+sample_1[rr];
    oFile.write('######  running tophat for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("mapping sample_1, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_1, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  if lite==1:
    return;
 
  logging.debug("mapping the second sample");
  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    cmd = 'tophat -a '+str(tophatAnchor)+' -m 0 -I 300000 -p 4 -g 20 --library-type '+ libType +' --no-novel-indels ';
    cmd += ' -N 3 --segment-mismatches 2 -G '+gtf+' -o '+rTempFolder;
    if SEPE=='PE': ## paired-end
      cmd += ' -r '+str(insertLength2[rr]) + ' --mate-std-dev ' + str(sigma2[rr])+' ' + bIndex +' '+sample_2[rr].split(':')[0]+' '+sample_2[rr].split(':')[1];
    else: ## single-end
      cmd += ' ' + bIndex +' '+sample_2[rr];
    oFile.write('######  running tophat for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("mapping sample_2, rep_"+str(rr+1)+" is done with status %s" % status);
    if (int(status)!=0): ## it did not go well
      logging.debug("error in mapping sample_2, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      raise Exception();
    logging.debug(output);

  return;
##### end of doTophatMapping ####

def indexBamFile(): ## indexing bam files to use pysam
  logging.debug("getting unique SAM function..");

  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    bam_fn='';
    if bamFile==0:  ## we know the location of the bam file
      bam_fn = rTempFolder+'/Aligned.sortedByCoord.out.bam';
    else: ## bam file is provided
      bam_fn = sample_1[rr];

    if LooseVersion(pysam.version.__samtools_version__) < LooseVersion('1.3'):
      pysam.sort(bam_fn, rTempFolder+'/aligned.sorted'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file
    else:
      pysam.sort(bam_fn, '-o', rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file

  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    rTempFolder = s2rPath+str(rr+1);
    bam_fn='';
    if bamFile==0:  ## we know the location of the bam file
      bam_fn = rTempFolder+'/Aligned.sortedByCoord.out.bam';
    else: ## bam file is provided
      bam_fn = sample_2[rr];
      
    if LooseVersion(pysam.version.__samtools_version__) < LooseVersion('1.3'):
      pysam.sort(bam_fn, rTempFolder+'/aligned.sorted'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file
    else:
      pysam.sort(bam_fn, '-o', rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam file
      pysam.index(rTempFolder+'/aligned.sorted.bam'); ## it will make aligned.sorted.bam.bai file

        

def getUniqueSAM(): ## getting uniquely mapped reads or pairs, DEPRECATED!!
  logging.debug("getting unique SAM function..");

  for rr in range(0,len(sample_1)): ## for each replicate of sample_1
    rTempFolder = s1rPath+str(rr+1);
    cmd = '';
    if SEPE=='PE': ## paired-end
      if bamFile==0:
        cmd += 'samtools view -h '+rTempFolder+'/Aligned.sortedByCoord.out.bam';
      else: ## bam file is provided
        cmd += 'samtools view -h '+sample_1[rr];
      cmd += ' | awk -F"\\t" \'(($0 ~ "NH:i:1[^0-9]"||$0 ~ "NH:i:1$") && ((and($2,0x2))&&($6=="'+str(readLength)+'M"||($6~"N"&&$6!~"D"&&$6!~"I"))) )|| NF<7\' > ' +rTempFolder+'/unique.S1.sam'; ## genome, junction reads and header
    else: ## single-end
      if bamFile==0:
        cmd += 'samtools view -h '+rTempFolder+'/Aligned.sortedByCoord.out.bam';
      else:
        cmd += 'samtools view -h '+sample_1[rr];
      cmd += ' | awk -F"\\t" \'(($0 ~ "NH:i:1[^0-9]"||$0 ~ "NH:i:1$") && ($6=="'+str(readLength)+'M"||($6~"N"&&$6!~"D"&&$6!~"I")) )|| NF<7\' > ' +rTempFolder+'/unique.S1.sam'; ## genome, junction reads and header
    oFile.write('######  getting unique reads or pairs for sample_1, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("getting uniquely mapped reads or pairs for sample_1, rep_%d is done with status %s" % ((rr+1),status));
    if (int(status)!=0): ## it did not go well
      logging.debug("error in getting uniquely mapped reads from smaple_1, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      logging.debug("Retry up to 3 more times..");
      for rp in range(0,3): ## try up to 3 more times
        logging.debug("Retry: " + str(rp+1));
        status,output=commands.getstatusoutput(cmd);
        logging.debug("getting uniquely mapped reads or pairs for sample_1, rep_%d is done with status %s" % ((rr+1),status));
        if (int(status)==0): ## worked okay
          break; ## break for loop
        else: ### error
          logging.debug("error in getting uniquely mapped reads from smaple_1, rep_%d: %s" % ((rr+1),status));
          logging.debug("error detail: %s" % output);
      if (int(status)!=0): ## didn't go well in all retries
        raise Exception();
    logging.debug(output);

  for rr in range(0,len(sample_2)): ## for each replicate of sample_2
    if lite==1:
      break;
    rTempFolder = s2rPath+str(rr+1);
    cmd = '';
    if SEPE=='PE': ## paired-end
      if bamFile==0:
        cmd += 'samtools view -h '+rTempFolder+'/Aligned.sortedByCoord.out.bam';
      else: ## bam file is provided
        cmd += 'samtools view -h '+sample_2[rr];
      cmd += ' | awk -F"\\t" \'(($0 ~ "NH:i:1[^0-9]"||$0 ~ "NH:i:1$") && ((and($2,0x2))&&($6=="'+str(readLength)+'M"||($6~"N"&&$6!~"D"&&$6!~"I"))) )|| NF<7\' > ' +rTempFolder+'/unique.S2.sam'; ## genome, junction reads and header
    else: ## single-end
      if bamFile==0:
        cmd += 'samtools view -h '+rTempFolder+'/Aligned.sortedByCoord.out.bam';
      else:
        cmd += 'samtools view -h '+sample_2[rr];
      cmd += ' | awk -F"\\t" \'(($0 ~ "NH:i:1[^0-9]"||$0 ~ "NH:i:1$") && ($6=="'+str(readLength)+'M"||($6~"N"&&$6!~"D"&&$6!~"I")) )|| NF<7\' > ' +rTempFolder+'/unique.S2.sam'; ## genome, junction reads and header
    oFile.write('######  getting unique reads or pairs for sample_2, replicate_'+ str(rr+1)+'#####\n'+cmd+'\n#\n');
    oFile.flush();
    status,output=commands.getstatusoutput(cmd);
    logging.debug("getting uniquely mapped reads or pairs for sample_2, rep_%d is done with status %s" % ((rr+1),status));
    if (int(status)!=0): ## it did not go well
      logging.debug("error in getting uniquely mapped reads from smaple_2, rep_%d: %s" % ((rr+1),status));
      logging.debug("error detail: %s" % output);
      logging.debug("Retry up to 3 more times..");
      for rp in range(0,3): ## try up to 3 more times
        logging.debug("Retry: " + str(rp+1));
        status,output=commands.getstatusoutput(cmd);
        logging.debug("getting uniquely mapped reads or pairs for sample_2, rep_%d is done with status %s" % ((rr+1),status));
        if (int(status)==0): ## worked okay
          break; ## break for loop
        else: ### error
          logging.debug("error in getting uniquely mapped reads from smaple_2, rep_%d: %s" % ((rr+1),status));
          logging.debug("error detail: %s" % output);
      if (int(status)!=0): ## didn't go well in all retries             
        raise Exception();
    logging.debug(output);

  return;
########## end of getUniqueSAM ####


def uniqSamNames(sampleNum): ## getting unique sam file names for sample_1 or sample_2

  rValue=[];

  if sampleNum==1: ## for the sample_1
    for rr in range(0,len(sample_1)): ## for each replicate of sample_1
      rTempFolder = s1rPath+str(rr+1);
      #rValue.append(rTempFolder+'/unique.S1.sam');
      rValue.append(rTempFolder+'/aligned.sorted.bam');
  elif sampleNum==2: ## for the sample_2
    for rr in range(0,len(sample_2)): ## for each replicate of sample_2
      rTempFolder = s2rPath+str(rr+1);
      #rValue.append(rTempFolder+'/unique.S2.sam');
      rValue.append(rTempFolder+'/aligned.sorted.bam');
  else: ## something is wrong here..
    logging.debug("error in uniqSamNames function. Incorrect sample number %d" % sampleNum);
    raise Exception();
  return rValue;
############### end of uniqSamNames #################

def getASEvents(asDir,): ## get AS events from GTF and SAM files
  logging.debug("getting AS events function..");

  samNames1=uniqSamNames(1);
  samNames2=uniqSamNames(2);
  finalSamNames = ','.join(samNames1)+','+','.join(samNames2);  
  if lite==1:
    finalSamNames = ','.join(samNames1); ## for the lite version

  #cmd = 'python '+binPath+'/processGTF.SAMs.py '+gtf+' '+asDir+'/fromGTF'+' '+finalSamNames+' '+libType+' '+tempPath;
  cmd = 'python '+binPath+'/processGTF.BAMs.py '+gtf+' '+asDir+'/fromGTF'+' '+finalSamNames+' '+libType+' '+SEPE+' '+str(novelSS) +' '+tempPath;
  #oFile.write('###### getting AS events from GTF and SAM files #####\n'+cmd+'\n#\n');
  oFile.write('###### getting AS events from GTF and BAM files #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("getting AS events is done with status %s" % status);
  if (int(status)!=0): ## it did not go well
    logging.debug("error in getting AS events %s" % status);
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  return;
############ end of getASEvents #####

def makingMATSInput(configFile): ## making MATS input 
  logging.debug("making MATS input function..");
  #cmd = 'python '+binPath+'/MATS.processsUnique.sam.py ' + configFile;
  cmd = 'python '+binPath+'/MATS.processsUnique.bam.py ' + configFile;
  oFile.write('###### making MATS input from AS events and BAM files #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("making MATS input is done with status %s" % status);
  if (int(status)!=0): ## it did not go well
    logging.debug("error in making MATS input %s" % status);
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  return;
############ end of makingMATSInput #####

def runningMATS(asType): ## running MATS here
  logging.debug("running MATS for " + asType + ". Using Junction Counts only");
  allInput = tempPath+"/JC.RNASeq."+asType+".MATS.input.txt"; ## input for all possible events
  bothIso = tempPath + "/"+asType + ".JC.input.txt";  ## events with both isoforms detected
  cmd = "awk '{ split($2,ic_1,\",\"); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,\",\"); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,\",\"); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,\",\"); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' " + allInput+" > " + bothIso;
  cmd += ";"+scriptPath+"/MATS/rMATS.sh -d "+bothIso+" -o "+tempPath+"/"+asType+"out_JC/ -c "+str(c)+" -p 4 -t "+analysis;
  oFile.write('###### running MATS input for ' + asType + ' using Junction Counts only #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("running MATS for %s using JC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in running MATS for %s with JCs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  logging.debug("running MATS for " + asType + ". Using Junction Counts and Reads on target Exon Counts");
  allInput = tempPath+"/JCEC.RNASeq."+asType+".MATS.input.txt"; ## input for all possible events
  bothIso = tempPath + "/"+asType + ".JCEC.input.txt";  ## events with both isoforms detected
  cmd = "awk '{ split($2,ic_1,\",\"); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,\",\"); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,\",\"); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,\",\"); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' " + allInput+" > " + bothIso;
  cmd += ";"+scriptPath+"/MATS/rMATS.sh -d "+bothIso+" -o "+tempPath+"/"+asType+"out_JCEC/ -c "+str(c)+" -p 8 -t "+analysis;
  oFile.write('###### running MATS input for ' + asType + ' using Junction Counts and Reads on target Exon Counts #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("running MATS for %s using JCEC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in running MATS for %s with JCECs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  return;
########## end of runningMATS ##########



def joiningMATS(asType): ## joining MATS output here
  allEvents = outPath+"/ASEvents/fromGTF."+asType+".txt"; ## all possible events info

  logging.debug("joining MATS for " + asType + ". Using Junction Counts only");
  outMATS = tempPath+"/"+asType+"out_JC/rMATS_Result.txt";
  tMATS = tempPath+"/"+asType+".MATS.JunctionCountOnly.txt";
  finalMATS = finalPath+"/"+asType+".MATS.JunctionCountOnly.txt";
  cmd = "python " + binPath+"/joinFiles.py " + allEvents + " " + outMATS + " " + " 0 0 " + tMATS; ## unsorted
  cmd += " ; awk '{print $(NF-3),$0}' "+tMATS+" | sort -g | cut -d' ' -f2- > "+finalMATS;
  oFile.write('###### joining MATS results for ' + asType + ' using Junction Counts only #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("joining MATS for %s using JC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in joining MATS for %s with JCs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  logging.debug("joining MATS for " + asType + ". Using Junction Counts and Reads on target Exon Counts");
  outMATS = tempPath+"/"+asType+"out_JCEC/rMATS_Result.txt";
  tMATS = tempPath+"/"+asType+".MATS.JunctionCountOnly.txt";
  finalMATS = finalPath+"/"+asType+".MATS.ReadsOnTargetAndJunctionCounts.txt";
  cmd = "python " + binPath+"/joinFiles.py " + allEvents + " " + outMATS + " " + " 0 0 " + tMATS; ## unsorted
  cmd += " ; awk '{print $(NF-3),$0}' "+tMATS+" | sort -g | cut -d' ' -f2- > "+finalMATS;
  oFile.write('###### joining MATS results for ' + asType + ' using Junction Counts and Reads on target Exon Counts #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("joining MATS for %s using JCEC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in joining MATS for %s with JCECs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  return;
################ end of joiningMATS ############



def liteJoiningMATS(asType): ## joining MATS lite here (for future release)
  allEvents = outPath+"/ASEvents/fromGTF."+asType+".txt"; ## all possible events info

  logging.debug("joining lite MATS for " + asType + ". Using Junction Counts only");
  eMATS = tempPath+"/JC.RNASeq."+asType+".MATS.input.txt";
  tMATS = tempPath+"/noS2.JC.RNASeq."+asType+".MATS.input.txt";
  finalMATS = finalPath+"/lite."+asType+".MATS.JunctionCountOnly.txt";
  cmd = "cut -f1,2,3,6,7 " + eMATS + " > " + tMATS; 
  cmd += ";paste " + allEvents + " " + tMATS + " > " + finalMATS;
  oFile.write('###### joining MATS lite results for ' + asType + ' using Junction Counts only #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("joining MATS lite for %s using JC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in joining MATS lite for %s with JCs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  logging.debug("joining lite MATS for " + asType + ". Using Junction Counts and Reads on target Exon Counts");
  eMATS = tempPath+"/JCEC.RNASeq."+asType+".MATS.input.txt";
  tMATS = tempPath+"/noS2.JCEC.RNASeq."+asType+".MATS.input.txt";
  finalMATS = finalPath+"/lite."+asType+".MATS.ReadsOnTargetAndJunctionCounts.txt";
  cmd = "cut -f1,2,3,6,7 " + eMATS + " > " + tMATS; 
  cmd += ";paste " + allEvents + " " + tMATS + " > " + finalMATS;
  oFile.write('###### joining MATS lite results for ' + asType + ' using Junction Counts and Reads on target Exon Counts #####\n'+cmd+'\n#\n');
  oFile.flush();
  status,output=commands.getstatusoutput(cmd);
  logging.debug("joining MATS lite for %s using JCEC is done with status %s" % (asType, status));
  if (int(status)!=0): ## it did not go well
    logging.debug("error in joining MATS for %s with JCECs with status %s" % (asType,status));
    logging.debug("error detail: %s" % output);
    raise Exception();
  logging.debug(output);

  return;
################ end of liteJoiningMATS ############


def printStats(asType, sF): ## printing stats here

  logging.debug("getting stats for " + asType + ". Using Junction Counts only.");
  finalMATS = finalPath+"/"+asType+".MATS.JunctionCountOnly.txt";
  cmd_0 = "wc -l " + finalMATS; ## total # of events
  cmd_1 = "awk '$(NF-3)<0.05 && $NF>0' "+finalMATS+" | wc -l"; ## Sample_1 has higher inclusion level
  cmd_2 = "awk '$(NF-3)<0.05 && $NF<0' "+finalMATS+" | wc -l"; ## Sample_2 has higher inclusion level
  status0,output0=commands.getstatusoutput(cmd_0);
  status1,output1=commands.getstatusoutput(cmd_1);
  status2,output2=commands.getstatusoutput(cmd_2);
  logging.debug("getting stats for %s using JC is done with status %s,%s and %s" % (asType, status0, status1, status2));
  if (int(status0)!=0): ## it did not go well
    logging.debug("error in getting status0 for %s with JCs. status: %s" % (asType,status0));
    logging.debug("error detail: %s" % output0);
    raise Exception();
  logging.debug("numLines: "+output0);
  if (int(status1)!=0): ## it did not go well
    logging.debug("error in getting status1 for %s with JCs. status: %s" % (asType,status1));
    logging.debug("error detail: %s" % output1);
    raise Exception();
  logging.debug("Upregulated: " + output1);
  if (int(status2)!=0): ## it did not go well
    logging.debug("error in getting status2 for %s with JCs. status: %s" % (asType,status2));
    logging.debug("error detail: %s" % output2);
    raise Exception();
  logging.debug("Downregulated: " + output2); 

  print('========== ' + asType + ' ========');
  #sF.write('========== ' + asType + ' ========\n');
  sFString = asType;
  logging.debug('========== ' + asType + ' ========');
  out0=str(int(output0.split()[0])-1);
  out1=output1.split()[0];
  out2=output2.split()[0];
  print('Junction Counts Only: There are %d AS events. Of these, %d events are statistically significant' % (int(out0),int(out1)+int(out2))); 
  print ('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2' % (int(out1),int(out2)));   
  #sF.write('Junction Counts Only: There are %d AS events. Of these, %d events are statistically significant\n' % (int(out0),int(out1)+int(out2))); 
  #sF.write('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2\n' % (int(out1),int(out2)));   
  sFString += '\t'+out0+'\t'+str(int(out1)+int(out2))+' ('+out1+':'+out2+')';
  logging.debug('Junction Counts Only: There are %d AS events. Of these, %d events are statistically significant' % (int(out0),int(out1)+int(out2))); 
  logging.debug('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2' % (int(out1),int(out2)));   


  logging.debug("getting stats for " + asType + ". Using Junction Counts and Reads on target Exon Counts");
  finalMATS = finalPath+"/"+asType+".MATS.ReadsOnTargetAndJunctionCounts.txt";
  cmd_0 = "wc -l " + finalMATS; ## total # of events
  cmd_1 = "awk '$(NF-3)<0.05 && $NF>0' "+finalMATS+" | wc -l"; ## Sample_1 has higher inclusion level
  cmd_2 = "awk '$(NF-3)<0.05 && $NF<0' "+finalMATS+" | wc -l"; ## Sample_2 has higher inclusion level
  status0,output0=commands.getstatusoutput(cmd_0);
  status1,output1=commands.getstatusoutput(cmd_1);
  status2,output2=commands.getstatusoutput(cmd_2);
  logging.debug("getting stats for %s using JCEC is done with status %s,%s and %s" % (asType, status0, status1, status2));
  if (int(status0)!=0): ## it did not go well
    logging.debug("error in getting status0 for %s with JCECs. status: %s" % (asType,status0));
    logging.debug("error detail: %s" % output0);
    raise Exception();
  logging.debug("numLines: " +output0);
  if (int(status1)!=0): ## it did not go well
    logging.debug("error in getting status1 for %s with JCECs. status: %s" % (asType,status1));
    logging.debug("error detail: %s" % output1);
    raise Exception();
  logging.debug("Upregulated: "+output1);
  if (int(status2)!=0): ## it did not go well
    logging.debug("error in getting status2 for %s with JCECs. status: %s" % (asType,status2));
    logging.debug("error detail: %s" % output2);
    raise Exception();
  logging.debug("Downregulated: "+output2);

  out0=str(int(output0.split()[0])-1);
  out1=output1.split()[0];
  out2=output2.split()[0];

  print('Junction Counts and Reads on target Exon Counts: There are %d AS events. Of these, %d events are statistically significant' % (int(out0),int(out1)+int(out2))); 
  print ('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2' % (int(out1),int(out2)));   
  #sF.write('Junction Counts and Reads on target Exon Counts: There are %d AS events. Of these, %d events are statistically significant\n' % (int(out0),int(out1)+int(out2))); 
  #sF.write('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2\n' % (int(out1),int(out2)));   
  sFString += '\t'+out0+'\t'+str(int(out1)+int(out2))+' ('+out1+':'+out2+')';
  logging.debug('Junction Counts and Reads on target Exon Counts: There are %d AS events. Of these, %d events are statistically significant' % (int(out0),int(out1)+int(out2))); 
  logging.debug('%d significant events have higher inclusion level for SAMPLE_1 and %d events for SAMPLE_2' % (int(out1),int(out2)));   

  sF.write(sFString+'\n');

  return;
###################### end of  printStats #############

######## end of functions ##############


################## actual process ##############

####
#### 1. tophat mapping (we now moved to STAR)
####
#logging.debug("start mapping..")
#try:
#  if bamFile==0: ## no bam file, start mapping
#    doTophatMapping();
#    pass;
#  else: ## bam file is provided
#    logging.debug("bam files are provided. skip mapping..");
#except:
#  logging.debug("There is an exception in mapping");
#  logging.debug("Exception: %s" % sys.exc_info()[0]);
#  logging.debug("Detail: %s" % sys.exc_info()[1]);
#  sys.exit(-1);
#logging.debug("done mapping..");

####
#### 1. STAR mapping
####
logging.debug("start mapping..")
try:
  if bamFile==0: ## no bam file, start mapping
    doSTARMapping();
    pass;
  else: ## bam file is provided
    logging.debug("bam files are provided. skip mapping..");
except:
  logging.debug("There is an exception in mapping");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  sys.exit(-1);
logging.debug("done mapping..");


####
#### 2. index bam files ###get uniquely mapped reads or pairs
####
#logging.debug("getting uniquely mapped reads or pairs");
logging.debug("indexing bam files to use pysam");

try:
  #getUniqueSAM();
  indexBamFile();
  pass;
except:
  logging.debug("There is an exception in indexing bam files");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  sys.exit(-2);
logging.debug("done indexing bam files..");


####
#### 3. detecting AS events
####
asFolder = outPath+'/ASEvents';
os.system('mkdir -p '+ asFolder);

logging.debug("start getting AS events from GTF and BAM files");
try:
  getASEvents(asFolder);
  pass;
except:
  logging.debug("There is an exception in getting AS events");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  sys.exit(-3);
logging.debug("done getting AS events..");

####
#### 4. making MATS input for each AS event
####
conFile = outPath + '/config.txt';
os.system('cp ' + scriptPath + '/data/config.sample.txt ' + conFile); ## getting config file
#
logging.debug("Setting proper string");
#
##### setting proper string
#
### setting sed string to support MacOS
#
sedString = "sed -i";
if MacOS==1: ## different sed syntax for MacOS
  sedString = "sed -i ''";
## READLENGTH, JUNCTIONLENGTH 
os.system(sedString+" 's/READLENGTH/"+str(readLength)+"/g' " + conFile);
os.system(sedString+" 's/JUNCTIONLENGTH/"+str(junctionLength)+"/g' " + conFile);
#
## SEPATH, MXEPATH, A5SS, A3SS, RI
os.system(sedString+" 's/SEPATH/"+asFolder.replace('/','\/') +"\\/fromGTF.SE.txt/g' " + conFile);
os.system(sedString+" 's/MXEPATH/"+asFolder.replace('/','\/') +"\\/fromGTF.MXE.txt/g' " + conFile);
os.system(sedString+" 's/A5SSPATH/"+asFolder.replace('/','\/') +"\\/fromGTF.A5SS.txt/g' " + conFile);
os.system(sedString+" 's/A3SSPATH/"+asFolder.replace('/','\/') +"\\/fromGTF.A3SS.txt/g' " + conFile);
os.system(sedString+" 's/RIPATH/"+asFolder.replace('/','\/') +"\\/fromGTF.RI.txt/g' " + conFile);
#
## SAMFOLDER
os.system(sedString+" 's/SAMFOLDER//g' " + conFile);
#
## OUTFOLDER
os.system(sedString+" 's/OUTFOLDER/"+tempPath.replace('/','\/')+"/g' " + conFile);
#
## LIBTYPE
os.system(sedString+" 's/LIBTYPE/"+libType+"/g' " + conFile);
#
## SEPE
if SEPE=='SE':
  os.system(sedString+" 's/SEPE/single/g' " + conFile);
else:
  os.system(sedString+" 's/SEPE/paired/g' " + conFile);
#
## INPUT_1, INPUT_2
#
sams_1=uniqSamNames(1); ## sam files for smaple_1
sams_2=uniqSamNames(2); ## sam files for sample_2

inString = ','.join(sams_1); ## sam files for smaple_1
os.system(sedString+" 's/INPUT_1/"+inString.replace('/','\/') +"/g' " + conFile);
inString = ','.join(sams_2); ## sam files for sample_2
os.system(sedString+" 's/INPUT_2/"+inString.replace('/','\/') +"/g' " + conFile);
#
### 
#
logging.debug("start making MATS input files from AS events and SAM files");
try:
  makingMATSInput(conFile);
  pass;
except:
  logging.debug("There is an exception in making MATS input");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  sys.exit(-4);
logging.debug("done making MATS input..");
#
#

########
######## 4.5 do lite operation here 
########
if lite==1:
  logging.debug("start joining lite results for each AS event");
  try:
    liteJoiningMATS("SE");         ## skipped exon
    liteJoiningMATS("MXE");        ## mutually exclusive
    liteJoiningMATS("A5SS");       ## alt-5
    liteJoiningMATS("A3SS");       ## alt-3
    liteJoiningMATS("RI");         ## retained intron
    pass;
  except:
    logging.debug("There is an exception in lite joining MATS");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-45);
  logging.debug("done lite joining MATS results..");
#


####
#### 5. running rMATS
####
#
if lite==0:
  logging.debug("start running MATS for each AS event");
  try:
    runningMATS("SE");         ## skipped exon
    runningMATS("MXE");        ## mutually exclusive 
    runningMATS("A5SS");       ## alt-5
    runningMATS("A3SS");       ## alt-3
    runningMATS("RI");         ## retained intron
    pass;
  except:
    logging.debug("There is an exception in running MATS");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-7);
  logging.debug("done running MATS for all AS event types..");
#


####
#### 6. joining rMATS input and output files
####
if lite==0:
  logging.debug("start joining MATS results for each AS event");
  try:
    joiningMATS("SE");         ## skipped exon
    joiningMATS("MXE");        ## mutually exclusive
    joiningMATS("A5SS");       ## alt-5
    joiningMATS("A3SS");       ## alt-3
    joiningMATS("RI");         ## retained intron
    pass;
  except:
    logging.debug("There is an exception in joining MATS");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-8);
  logging.debug("done joining MATS results..");
#


####
#### 7. print out stats
####
if lite==0:
  summaryFile = open(outPath+'/summary.txt', 'w');
  summaryFile.write("\n################### folder names and associated input files #############\n");
  for fki in range(0,len(sample_1)): ## for each replicate of sample_1
    repTempFolder = "SAMPLE_1\REP_"+str(fki+1);
    associatedFile = sample_1[fki];
    summaryFile.write(repTempFolder+"\t"+associatedFile+"\n");
  
  for fki in range(0,len(sample_2)): ## for each replicate of sample_2
    repTempFolder = "SAMPLE_2\REP_"+str(fki+1);
    associatedFile = sample_2[fki];
    summaryFile.write(repTempFolder+"\t"+associatedFile+"\n");
  
  summaryFile.write("#########################################################################\n");
  
  summaryFile.write("\n############################# MATS Report #############################\n");
  summaryFile.write("="*70+'\n');
  sumFileHeader = "EventType\tNumEvents.JC.only\tSigEvents.JC.only\tNumEvents.JC+readsOnTarget\tSigEvents.JC+readsOnTarget";
  summaryFile.write(sumFileHeader+"\n");
  summaryFile.write("="*70+'\n');
  logging.debug("======================= Final Report =============");
  try:
    printStats("SE", summaryFile);         ## skipped exon
    printStats("MXE", summaryFile);        ## mutually exclusive
    printStats("A5SS", summaryFile);       ## alt-5
    printStats("A3SS", summaryFile);       ## alt-3
    printStats("RI", summaryFile);         ## retained intron
  except:
    logging.debug("There is an exception in printingStats");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-9);
  logging.debug("done printing out stats..");
  #
  summaryFile.write("="*70+'\n');
  #summaryFile.write("#######################################################################\n");
  summaryFile.write("################# Report Legend ################\nEventType: Type of AS event\n\t");
  summaryFile.write("SE: Skipped exon\n\tMXE: Mutually exclusive exon\n\tA5SS: Alternative 5' splice site\n\t");
  summaryFile.write("A3SS: Alternative 3' splice site\n\tRI: Retained intron\n");
  summaryFile.write("NumEvents.JC.only: total number of events detected using Junction Counts only\n");
  summaryFile.write("SigEvents.JC.only: number of significant events detected using Junction Counts only\n\t");
  summaryFile.write("The numbers in the parentheses (n1:n2) indicate the number of significant events that have higher inclusion level for SAMPLE_1 (n1) or for SAMPLE_2 (n2)\n");   
  summaryFile.write("NumEvents.JC+readsOnTarget: total number of events detected using both Junction Counts and reads on target\n");
  summaryFile.write("SigEvents.JC+readsOnTarget: number of significant events detected using both Junction Counts and reads on target\n\t");
  summaryFile.write("The numbers in the parentheses (n1:n2) indicate the number of significant events that have higher inclusion level for SAMPLE_1 (n1) or for SAMPLE_2 (n2)\n");   
  summaryFile.write("################################################\n");
  
  summaryFile.close();
  
  oFile.close();

### clean up temp folder if necessary ###
if keepTemp==0: ## delete temp folder, by default
  os.system('rm -rf '+ tempPath);
  logging.debug("Temp folder is deleted..");  
#
#### message to suggest sam file deletion
#
##print "\n#####\nPlease delete sam files in "+outPath+"/SAMPLE_*/REP_*/ folder to save storage space!\n#####";
#
#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
