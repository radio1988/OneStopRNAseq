mkdir -p chr19 mkdir fastq

#conda activate bio # samtools, bioconda::bedtools
# workdir: where bam are located
# outputdir: ./chr19 ./fastq
# function: get chr19 from bam, export R1 R2 as fastq.gz, for example-test data
# why chr19: shortest in mm10, pretty short in hs38 too

for f in *bam
do 
samtools view -h -@ 4 $f 19 | samtools sort -n -@ 4 -o chr19/$f 
bamToFastq -i chr19/$f -fq fastq/${f/bam/R1.fastq.gz} -fq2 fastq/${f/bam/R2.fastq.gz}
done
