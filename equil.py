import os
from simsopt.util.mpi import MpiPartition
from simsopt.mhd import Vmec

mpi = MpiPartition(ngroups=3)

path='/home/matteo/stellarator/pyQSC/config_gen/vmec_files'
vmecfiles=[]

for file in os.listdir(path):
	if file.startswith('input.'):
		vmecfiles.append(file)
 
if not os.path.exists(os.getcwd()+'/outputs'):
	os.mkdir(os.getcwd()+'/outputs')

os.chdir(os.getcwd()+'/outputs')
           
for vmecfile in vmecfiles:
	equil = Vmec(path+'/'+vmecfile, mpi)
	surf = equil.boundary
	surf.fix_all()

