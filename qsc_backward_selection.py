""" Author:	Matteo Notazio

"""

import os
import f90nml
import numpy
import math
from simsopt.util.mpi import MpiPartition
from simsopt.mhd import Vmec
from simsopt.geo.surface import Surface

mpi = MpiPartition(ngroups=3)

# PATH OF VMEC FILES GENERATED USING pyQSC
path='/home/matteo/stellarator/pyQSC/config_gen/vmec_files'

vmecfiles=[]

# SELECT ALL THE VMEC FILES STORED INTO THE FOLDER
for file in os.listdir(path):
	if file.startswith('input.stellqsc'):
		vmecfiles.append(file)

# GENERATE A FOLDER TO STORE THE PATCHED VMEC FILES
if not os.path.exists(os.getcwd()+'/vmec_files_patched'):
	os.mkdir(os.getcwd()+'/vmec_files_patched')
                                                                                                                                                              
os.chdir(os.getcwd()+'/vmec_files_patched')
                                                                               
for vmecfile in vmecfiles:
	
	vmec = Vmec(path+'/'+vmecfile, mpi)
		
	if vmec is not None:
		print('\n'+vmecfile+' opened successfully!\n')
		
		# NUMBER OF TOROIDAL AND POLOIDAL FIELD PERIODS
		ntor=vmec.indata.ntor
		mpol=vmec.indata.mpol

		# GAMMA CURVE DEFINING THE BOUNDARY POINTS (x,y,z)
		gamma=vmec.boundary.gamma()
		
		# BOUNDARY ASPECT RATIO
		ar=vmec.boundary.aspect_ratio()
		
		mean_div=0
		aa=0
		index_out_of_range=False
		
		while mean_div<1 and index_out_of_range==False:
		
			if mean_div<1 and index_out_of_range==False:
			
				if aa<mpol:
								
					# SET TO 0 THE FOURIER COEFFICIENTS STARTING FROM THE LAST ROW
					for phi in range(0,ntor+1):
						for bb in range(aa):
							vmec.boundary.set_rc(mpol-bb,phi,0)
							vmec.boundary.set_zs(mpol-bb,phi,0)
						
					for phi in range(-ntor,1):
						for bb in range(aa):
							vmec.boundary.set_rc(mpol-bb,phi,0)
							vmec.boundary.set_zs(mpol-bb,phi,0)
					
					""" # PRINT THE FOURIER COEFFICIENTS	
					for theta in range(mpol+1):
						for phi in range(0,ntor+1):
							print('rbc['+str(phi)+']['+str(theta)+']='+str(vmec.boundary.get_rc(theta,phi))+'\n')
						for phi in range(-ntor,1):
							print('rbc['+str(phi)+']['+str(theta)+']='+str(vmec.boundary.get_rc(theta,phi))+'\n')
					"""
					
					# GAMMA CURVE DEFINING THE PATCHED BOUNDARY POINTS (x,y,z)
					gamma_patched=vmec.boundary.gamma()
					
					# PATCHED BOUNDARY ASPECT RATIO
					ar_patched=vmec.boundary.aspect_ratio()
					
					# CALCULATE THE POINT-TO-POINT DIVERGENCE
					if len(gamma)==len(gamma_patched) and len(gamma[0])==len(gamma_patched[0]):
						div=[]
						for phi in range(len(gamma)):
							for theta in range(len(gamma[0])):
								x=gamma[phi][theta][0]
								y=gamma[phi][theta][1]
								z=gamma[phi][theta][2]
								point=math.sqrt(x**2+y**2+z**2)
								    
								x_patched=gamma_patched[phi][theta][0]
								y_patched=gamma_patched[phi][theta][1]
								z_patched=gamma_patched[phi][theta][2]
								point_patched=math.sqrt(x_patched**2+y_patched**2+z_patched**2)
									
								div.append((abs(point-point_patched)/point)*100)
						
						# MEAN POINT-TO-POINT DIVERGENCE
						mean_div=numpy.mean(div)

						print('\n'+vmecfile+' ASPECT RATIO: '+str(ar)+'\n')
						print(vmecfile+'_patched ASPECT RATIO: '+str(ar_patched)+'\n')

						print('POINTS MEAN DIVERGENCE: '+str(mean_div)+' %\n')
						
						aa+=1
				else:
					index_out_of_range=True
					print('ALLERT! INDEX m CANNOT BE REDUCED ANYMORE!')
		
		# WRITE A NEW VMEC FILE CONTAINING THE PATCHED BOUNDARY
		vmec.write_input(vmecfile+'_patched')
		
		vmec_patched=Vmec(vmecfile+'_patched')
		
		""" # PLOT THE BOUNDARIES	
		vmec.boundary.plot()
		vmec_patched.boundary.plot()
		"""
