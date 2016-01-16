declare -a arr=(2 4 8 16 24 48 72 96 120 144 168)
## now loop through the above array
for nproc in "${arr[@]}"
do
   dir="lev9_10_test$nproc"
   mkdir $dir
   cd $dir
   #Call python script that generates output of new run.pbs
   python ../getscript.py $nproc > run.pbs
   qsub run.pbs
   cd .. 
done
