#declare -a arr=(1  24 48 )
declare -a arr=(1  24 48 120 240 360 480 720 960 1200 1440 1680  1920 2400) 
for nproc in "${arr[@]}"
do
   dir="test$nproc"
   cp -r ../ramses8 $dir 
   #mkdir $dir
   cd $dir
   #Call python script that generates output of new run.pbs
   /global/homes/d/dorislee/anaconda/bin/python ../getscript.py $nproc > run.pbs
   echo "Submitting $nproc"
   sbatch run.pbs
   cd .. 
done
