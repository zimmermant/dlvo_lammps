from mpi4py import MPI
from lammps import lammps
lmp = lammps()
lmp.file("in.colloid")
me=MPI.COMM_WORLD.Get_rank()
nprocs=MPI.COMM_WORLD.Get_size()
print "Proc %d out of %d procs has " % (me,nprocs),lmp
MPI.Finalize()
