dir=mccb/green/minggang/2008_minggang_RNASeq/snakemake/
echo "uploading files to $dir"
dropbox_uploader.sh upload fastq $dir
dropbox_uploader.sh upload sorted_reads $dir
dropbox_uploader.sh upload bigWig $dir
dropbox_uploader.sh upload gsea_compressed $dir

