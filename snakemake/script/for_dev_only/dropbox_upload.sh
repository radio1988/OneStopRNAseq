# upload larger files to dropbox dropbox_pathectly 
# so that smaller files can be downloaded via rsync_download.sh or website
# todo: running from HPC job submission can't find dropbox_uploader and shasum
local_path=/home/rl44w/mount/onestoprnaseq_test/Sunil_GSE121103_05.02.47-08.17.2020/analysis_1/
dropbox_path=mccb/green/sunil/GSE121103_RNAseq/snakemake/
echo "uploading files from: $local_path"
echo "to: $dropbox_path"
echo "should run from local"
cd $local_path
dropbox_uploader.sh upload fastq $dropbox_path
dropbox_uploader.sh upload sorted_reads $dropbox_path
dropbox_uploader.sh upload bigWig $dropbox_path
dropbox_uploader.sh upload gsea_compressed $dropbox_path

