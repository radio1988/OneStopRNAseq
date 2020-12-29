rsync -aP --delete --exclude='envs*' ~/develop/ ~/mount/onestoprnaseq_test/develop_test_run2/
cd  $HOME/mount/onestoprnaseq_test/develop_test_run2/ 
ln -s $HOME/develop/envs
echo "Done"

