#
## this program trims fastq file to the given length
#

### import necessary libraries
import re,os,sys,logging,time,datetime;

### checking out the number of arguments
if (len(sys.argv)<4): 
  print('Not enough arguments!!');
  print ('It takes at least 3 arguments.');
  print ('Usage:\n\tProgramName inputFastq outputFastq desiredLength');
  print ('Example\n\tProgramName in.fastq out.fastq 32');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename='log.trimFastq.v2.'+  str(datetime.datetime.now())+'.txt',
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

###
MAX_READ_LENGTH=200;

iFile = open(sys.argv[1]); ## input fastq file
oFile = open(sys.argv[2], 'w'); ## out fastq file
seqLen = int(sys.argv[3]); ## desired read length, put 'N' with quality '#' if the given read is short

reads={}; ## dictioinary for stats
for i in range(MAX_READ_LENGTH):
  reads[i]=0; ## initialization;
frontClip=6; ## bps clipped from the beginning of the reads if a read>desiredLength;
c=0;

for l in iFile: ## for each line
  c+=1;
  
  line=l.strip();
  if(c%2==0):  ## sequence or quality
    s='';
    sLen=len(line);
    if sLen<seqLen: ## need to fill it up
      if (c%4==0): ## quality
        s=line[0:sLen]+'#'*(seqLen-sLen);
      else: ## sequence
        s=line[0:sLen]+'N'*(seqLen-sLen);
        reads[sLen]+=1;
    else: ## remove frontClip
      s=line[min(sLen-seqLen, frontClip): seqLen+(min(sLen-seqLen, frontClip))];
      if (c%4==2):  
        reads[sLen]+=1;

    oFile.write(s+'\n');
  else: ## IDs
    oFile.write(line+'\n');

iFile.close();
oFile.close();

logging.debug("Done processing %d reads" % int(c/4));
logging.debug("Read distribution");
logging.debug("length\tfrequency");
for i in range(MAX_READ_LENGTH):
  logging.debug("%3d\t%d" % (i,reads[i]));

#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
