""" 
****************************************************
University of Tuscia & Renaissance Fusion

Author:	Matteo Notazio, University of Tuscia 
****************************************************
"""
from datetime import date,datetime
import os
import numpy
import math
import matplotlib.pyplot as plt
from simsopt.util.mpi import MpiPartition
from simsopt.mhd import Vmec

def readme():
	f=open("README.me","w+")
	f.write('****************************************************\n\
	University of Tuscia & Renaissance Fusion\n\n\
	Author:\tMatteo Notazio, University of Tuscia\n\
	****************************************************\n\n\
	Reference VMEC files path:\t'+vmec_files_path)
	f.close()


####################################################################################################
mpi = MpiPartition(ngroups=3)

path=os.getcwd()
date=date.today()
time=datetime.now()

vmec_files_path=path+'/vmec_files/30052022_104747'
vmecfiles=[]

# SELECT ALL THE VMEC FILES STORED INTO THE FOLDER
for file in os.listdir(vmec_files_path):
	if file.startswith('input.stellqsc'):
		vmecfiles.append(file)
	
# GENERATE A FOLDER TO STORE THREED1 FILES
if not os.path.exists(path+'/threed1_files'):
	os.mkdir(path+'/threed1_files')
	
threed1_files_path=path+'/threed1_files/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
os.mkdir(threed1_files_path)
os.chdir(threed1_files_path)
readme()

# GENERATE A FOLDER TO STORE PATCHED VMEC FILES
if not os.path.exists(path+'/vmec_files_patched'):
	os.mkdir(path+'/vmec_files_patched')

vmec_files_patched_path=path+'/vmec_files_patched/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
os.mkdir(vmec_files_patched_path)
os.chdir(vmec_files_patched_path)
readme()

# GENERATE A FOLDER TO STORE VMEC FILES GRAPHS	 
if not os.path.exists(path+'/vmec_files_graphs'):
	os.mkdir(path+'/vmec_files_graphs')

vmec_files_graphs_path=path+'/vmec_files_graphs/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
os.mkdir(vmec_files_graphs_path)
os.chdir(vmec_files_graphs_path)
readme()
####################################################################################################

norm_min=1e-1 
norm_trend=[]
                                                                                                                                                                                                                                           
for vmecfile in vmecfiles:
	
	os.chdir(threed1_files_path)

	# ORIGINAL BOUNDARY
	vmec = Vmec(vmec_files_path+'/'+vmecfile, mpi)
		
	if vmec is not None:
	
		print('\n'+vmecfile+' opened successfully!\n')

		# NUMBER OF TOROIDAL AND POLOIDAL FIELD PERIODS
		ntor=vmec.indata.ntor
		mpol=vmec.indata.mpol
		
		# TOTAL NUMBER OF FOURIER COEFFICIENTS
		n_fourier_tot=((2*ntor)+1)*mpol

		# GAMMA CURVE DEFINING THE BOUNDARY POINTS (x,y,z)
		gamma=vmec.boundary.gamma()
		
		# BOUNDARY ASPECT RATIO
		ar=vmec.boundary.aspect_ratio()
		
		norm=1e-10
		
		aa=0
		index_out_of_range=False
		
		while norm<norm_min and index_out_of_range==False:
		
			if norm<norm_min and index_out_of_range==False:
			
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
					
					# NUMBER OF FOURIER COEFFICIENTS
					n_fourier=((2*ntor)+1)*(mpol-aa)
					
					os.chdir(vmec_files_patched_path)
					
					# WRITE A NEW VMEC FILE CONTAINING THE PATCHED BOUNDARY		
					vmec.write_input(vmecfile+'_patched')
					
					os.chdir(threed1_files_path)
					
					# PATCHED BOUNDARY
					vmec_patched=Vmec(vmec_files_patched_path+'/'+vmecfile+'_patched')
				
					# GAMMA CURVE DEFINING THE PATCHED BOUNDARY POINTS (x,y,z)
					gamma_patched=vmec_patched.boundary.gamma()
					
					# PATCHED BOUNDARY ASPECT RATIO
					ar_patched=vmec_patched.boundary.aspect_ratio()
					
					# CALCULATE THE POINT-TO-POINT DIVERGENCE (NORM)
					if len(gamma)==len(gamma_patched) and len(gamma[0])==len(gamma_patched[0]):
						
						tmp=[]
						
						gamma_diff=numpy.subtract(gamma,gamma_patched)

						norm=numpy.linalg.norm(gamma_diff)

						tmp.append(n_fourier)
						tmp.append(norm)
						norm_trend.append(tmp)

						print('\n#############################################\n')
						print(vmecfile+' ASPECT RATIO: '+str(ar)+'\n')
						print(vmecfile+'_patched ASPECT RATIO: '+str(ar_patched)+'\n')
						print('NORM: '+str(norm)+'\n')
						print('#############################################\n')
						
						aa+=1
				else:
					index_out_of_range=True
					print('ALLERT! INDEX m CANNOT BE REDUCED ANYMORE!')
		
		
		
		# PLOT BOUNDARIES
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.scatter(gamma[:, :, 0], gamma[:, :, 1], gamma[:, :, 2],c='red',s=10)
		ax.scatter(gamma_patched[:, :, 0], gamma_patched[:, :, 1], gamma_patched[:, :, 2],c='blue',s=10)
			
		ax.set_title(vmecfile+' - '+vmecfile+'_patched')
		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')
			
		ax.set_xlim(-1,1)
		ax.set_ylim(-1,1)
		ax.set_zlim(-1,1)
		
		ax.legend(['Original boundary', 'Patched boundary'])
		
		os.chdir(vmec_files_graphs_path)
		
		plt.savefig(vmecfile+'_bs.jpg')
		plt.close()
		

# PLOT NORM TREND
fig = plt.figure()
ax = plt.axes()
x=[x[0] for x in norm_trend]
y=[y[1] for y in norm_trend]
ax.scatter(x,y,c='blue')
ax.plot(x,[numpy.mean(y)]*len(x),c='red')

ax.set_title('Norm trend')
ax.set_xlabel('Number of Fourier coefficients')
ax.set_ylabel('Norm')

ax.set_xlim(1,max(x))
ax.set_ylim(min(y),max(y))

plt.savefig('n_fourier_norm.jpg')
plt.close()		

