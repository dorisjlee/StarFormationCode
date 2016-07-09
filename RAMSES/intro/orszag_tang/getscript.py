import sys

for arg in sys.argv: 
    nproc=arg
nproc = int(nproc)
print "#NPROC {} TEST".format(nproc)
print "#PBS -q regular"
print "#PBS -l mppwidth={}".format(((nproc/24)+1)*24)
print "#PBS -l walltime=03:00:00"
print "#PBS -N {}_test".format(nproc)
print "#PBS -e test.$PBS_JOBID.err"
print "#PBS -o test.$PBS_JOBID.out"
print "#PBS -A m2218"
print "module swap PrgEnv-gnu PrgEnv-intel"
print "cd /scratch/scratchdirs/dorislee/test{}".format(nproc)
print "pwd"
print "aprun -n {} /scratch/scratchdirs/dorislee/ramses/trunk/ramses/bin/ramses2d_mhd_otpatch_hopper ../orszag-tang_mod.nml".format(nproc)
