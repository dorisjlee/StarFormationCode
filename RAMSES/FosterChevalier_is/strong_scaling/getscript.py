import sys
for arg in sys.argv: 
    k=arg
k = int(k)
print "#!/bin/bash -l"
print "#SBATCH -p regular"
print "#SBATCH -n {}".format(k)
print "#SBATCH -o %j.out"
print "#SBATCH -e %j.err"
print "#SBATCH -t 03:00:00"
print "#SBATCH -J nproc{}_run".format(k)
print "cd  ~/project/strong_scaling/test{}/trunk/ramses/bin/".format(k) 
print "srun -n {0} --unbuffered ./ramses3d ../patch/hydro/isothermal_sphere/fc.nml".format(k)
