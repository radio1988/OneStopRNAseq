# for internal use only
# for smaller result files
echo "should run from personal computer viewing the files"
server_path=/home/rl44w/mount/onestoprnaseq_test/Sunil_GSE121103_05.02.47-08.17.2020/analysis_1/
local_path="./" # can be the dropbox local folder, don't run with dropbox_uploader.sh at the same time
rsync --max-size=100M --exclude=".snake*" --exclude="gsea/*" -aP rl44w@ghpcc06.umassrc.org:$server_path $local_path
