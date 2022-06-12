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
from matplotlib import cm
from simsopt.util.mpi import MpiPartition
from simsopt.mhd import Vmec

#********************FUNCTIONS#********************

def readme(vmec_files_path):
	f=open("README.me","w+")
	
	f.write('****************************************************\n\
	University of Tuscia & Renaissance Fusion\n\n\
	Author:\tMatteo Notazio, University of Tuscia\n\
	****************************************************\n\n\
	Reference VMEC files path:\t'+vmec_files_path)
	
	f.close()
	
def create_folders(vmec_files_dir,date,time):
	path=os.getcwd()

	vmec_files_path=path+vmec_files_dir
		
	# FOLDER TO STORE THREED1 FILES
	if not os.path.exists(path+'/threed1_files'):
		os.mkdir(path+'/threed1_files')
		
	threed1_files_path=path+'/threed1_files/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
	
	os.mkdir(threed1_files_path)
	os.chdir(threed1_files_path)
	
	readme(vmec_files_path)

	# FOLDER TO STORE PATCHED VMEC FILES
	if not os.path.exists(path+'/vmec_files_patched'):
		os.mkdir(path+'/vmec_files_patched')

	vmec_files_patched_path=path+'/vmec_files_patched/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
	
	os.mkdir(vmec_files_patched_path)
	os.chdir(vmec_files_patched_path)
	
	readme(vmec_files_path)

	# FOLDER TO STORE VMEC FILES GRAPHS	 
	if not os.path.exists(path+'/vmec_files_graphs'):
		os.mkdir(path+'/vmec_files_graphs')

	vmec_files_graphs_path=path+'/vmec_files_graphs/'+date.strftime('%d%m%Y')+'_'+time.strftime('%H%M%S')
	
	os.mkdir(vmec_files_graphs_path)
	os.chdir(vmec_files_graphs_path)
	
	readme(vmec_files_path)
	
	return vmec_files_path,threed1_files_path,vmec_files_patched_path,vmec_files_graphs_path


def get_vmec_files(vmec_files_dir):

	path=os.getcwd()

	vmec_files_path=path+vmec_files_dir
	vmec_files=[]

	for file in os.listdir(vmec_files_path):
		if file.startswith('input.stellqsc'):
			vmec_files.append(file)
	
	return vmec_files
	
def hist3d(x,y,z_size,title,zlabel,filename):
	fig = plt.figure()
	ax = plt.axes(projection="3d")

	z = [0]*len(z_size)
	
	x_size = numpy.ones(len(z_size))
	y_size = numpy.ones(len(z_size))
	
	z_max=0.001
	
	z_size=numpy.where(z_size>z_max,z_max,z_size)

	ax.bar3d(x, y, z, x_size, y_size, z_size,color='aqua')
	ax.set_title(title)
	ax.set_xlabel('m')
	ax.set_ylabel('n')
	ax.set_zlabel(zlabel)

	ax.set_zlim(0,z_max)
	
	ax.set_xticks(numpy.arange(0, mpol+1, 1))
	ax.set_yticks(numpy.arange(-ntor, ntor+2, 2))
	ax.set_zticks(numpy.arange(0, z_max, z_max/10))
	
	os.chdir(vmec_files_graphs_path)
	
	#plt.show()	
	plt.savefig(filename,bbox_inches='tight')
	plt.close()

def plot_boundaries(x,y,z,x_patched,y_patched,z_patched,x_lim,y_lim,z_lim,filename):

	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.scatter(x,y,z,c='red',s=10,alpha=0.5)
	ax.scatter(x_patched,y_patched,z_patched,c='blue',s=10,alpha=0.2)
			
	ax.set_title(vmec_file+' - '+vmec_file+'_patched\nBoundary')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
			
	ax.set_xlim(x_lim)
	ax.set_ylim(y_lim)
	ax.set_zlim(z_lim)
		
	ax.legend(['Original boundary', 'Patched boundary'])
		
	os.chdir(vmec_files_graphs_path)

	plt.savefig(filename,bbox_inches='tight')
	plt.close()

	

#********************MAIN********************

mpi = MpiPartition(ngroups=3)

vmec_files_dir='/vmec_files/11062022_115600'

vmec_files=get_vmec_files(vmec_files_dir)

[vmec_files_path,threed1_files_path,vmec_files_patched_path,vmec_files_graphs_path]=create_folders(vmec_files_dir,date.today(),datetime.now())

min_mean_error=1e-1 
error_trend_tot=[]
                                                                                                                                                                                                                                           
for vmec_file in vmec_files:
	
	os.chdir(threed1_files_path)

	# ORIGINAL BOUNDARY
	vmec = Vmec(vmec_files_path+'/'+vmec_file, mpi)
		
	if vmec is not None:
	
		print('\n'+vmec_file+' opened successfully!\n')
		
		# NUMBER OF FIELD PERIODS
		nfp=vmec.indata.nfp

		# NUMBER OF TOROIDAL AND POLOIDAL FIELD DIRECTIONS
		ntor=vmec.indata.ntor
		mpol=vmec.indata.mpol
		
		# TOTAL NUMBER OF FOURIER COEFFICIENTS
		n_fourier_tot=((2*ntor)+1)*mpol
				
		# GAMMA CURVE DEFINING THE BOUNDARY POINTS (x,y,z)
		gamma=vmec.boundary.gamma()
		
		# BOUNDARY ASPECT RATIO
		ar=vmec.boundary.aspect_ratio()
			
		R=numpy.zeros((mpol,2*ntor+1,3)) 
		Z=numpy.zeros((mpol,2*ntor+1,3)) 
		
		for m in range(mpol):
			for n in range(-ntor,ntor+1):
				R[m,n,0]=m
				R[m,n,1]=n
				R[m,n,2]=numpy.abs(vmec.boundary.get_rc(m,n))
				
				Z[m,n,0]=m
				Z[m,n,1]=n
				Z[m,n,2]=numpy.abs(vmec.boundary.get_zs(m,n))

		hist3d(R[:, :, 0].reshape((mpol)*(2*ntor+1)),R[:, :, 1].reshape((mpol)*(2*ntor+1)),R[:, :, 2].reshape((mpol)*(2*ntor+1)),
		vmec_file+'\n$|R_{m,n}(\\theta,\phi)|$','$|R_{m,n}|$',vmec_file+'_R.png')
		
		hist3d(Z[:, :, 0].reshape((mpol)*(2*ntor+1)),Z[:, :, 1].reshape((mpol)*(2*ntor+1)),Z[:, :, 2].reshape((mpol)*(2*ntor+1)),
		vmec_file+'\n$|Z_{m,n}(\\theta,\phi)|$','$|Z_{m,n}|$',vmec_file+'_Z.png')
		
		mean_error=1e-10
		
		trim_m=0
		trim_n=1
		index_out_of_range=False
		
		n_fourier=n_fourier_tot
		error_trend=[]
		
		while mean_error<min_mean_error and index_out_of_range==False:
		
			if mean_error<min_mean_error and index_out_of_range==False:

				if trim_m<mpol:
				
					# SET TO 0 ONE COEFFICIENT AT THE TIME
					if trim_n<2*ntor+1:
					
						if trim_m==0:
					
							if trim_n<=ntor:
						
								for n in range(trim_n+1):
										vmec.boundary.set_rc(mpol-1,-n,0)
										vmec.boundary.set_zs(mpol-1,-n,0)
							else:
							
								for n in range(-ntor,0):
										vmec.boundary.set_rc(mpol-1,-n,0)
										vmec.boundary.set_zs(mpol-1,-n,0)
										
								for n in range((trim_n-ntor)+1):
										vmec.boundary.set_rc(mpol-1,ntor-n,0)
										vmec.boundary.set_zs(mpol-1,ntor-n,0)
							
						else:
							deleted_rows=int(trim_n/(2*ntor+1))
						
							for m in range(trim_m):
								for n in range(-ntor,ntor+1):
									vmec.boundary.set_rc((mpol-m)-1,n,0)
									vmec.boundary.set_zs((mpol-m)-1,n,0)
								
							cols_to_delete=trim_n-(deleted_rows*(2*ntor+1))
								
							if cols_to_delete<=ntor:
						
								for n in range(trim_n+1):
										vmec.boundary.set_rc(mpol-(trim_m+1),-n,0)
										vmec.boundary.set_zs(mpol-(trim_m+1),-n,0)
							else:
							
								for n in range(-ntor,0):
										vmec.boundary.set_rc(mpol-(trim_m+1),-n,0)
										vmec.boundary.set_zs(mpol-(trim_m+1),-n,0)
											
								for n in range((cols_to_delete-ntor)+1):
										vmec.boundary.set_rc(mpol-(trim_m+1),ntor-n,0)
										vmec.boundary.set_zs(mpol-(trim_m+1),ntor-n,0)
	
					else:
						trim_m+=1
						trim_n=0
	
					# NUMBER OF FOURIER COEFFICIENTS
					n_fourier-=1
					
					last_m=mpol-trim_m
					if trim_n<ntor:
						last_n=-(trim_n+1)
					elif trim_n>=ntor and trim_n<=2*ntor+1:
						last_n=2*ntor-trim_n
					os.chdir(vmec_files_patched_path)
					
					# WRITE A NEW VMEC FILE CONTAINING THE PATCHED BOUNDARY		
					vmec.write_input(vmec_file+'_patched')
					
					os.chdir(threed1_files_path)
					
					# PATCHED BOUNDARY
					vmec_patched=Vmec(vmec_files_patched_path+'/'+vmec_file+'_patched')
					
					# NUMBER OF FIELD PERIODS
					nfp_patched=vmec.indata.nfp

					# NUMBER OF TOROIDAL AND POLOIDAL FIELD PERIODS
					ntor_patched=vmec.indata.ntor
					mpol_patched=vmec.indata.mpol
				
					# GAMMA CURVE DEFINING THE PATCHED BOUNDARY POINTS (x,y,z)
					gamma_patched=vmec_patched.boundary.gamma()
					
					# PATCHED BOUNDARY ASPECT RATIO
					ar_patched=vmec_patched.boundary.aspect_ratio()

					
					R_patched=numpy.zeros((mpol,2*ntor+1,3)) 
					Z_patched=numpy.zeros((mpol,2*ntor+1,3)) 
					
					for m in range(mpol):
						for n in range(-ntor,ntor+1):
							R_patched[m,n,0]=m
							R_patched[m,n,1]=n
							R_patched[m,n,2]=numpy.abs(vmec.boundary.get_rc(m,n))
							
							Z_patched[m,n,0]=m
							Z_patched[m,n,1]=n
							Z_patched[m,n,2]=numpy.abs(vmec.boundary.get_zs(m,n))

					
					hist3d(R_patched[:, :, 0].reshape((mpol)*(2*ntor+1)),R_patched[:, :, 1].reshape((mpol)*(2*ntor+1)),R_patched[:, :, 2].reshape((mpol)*(2*ntor+1)),
					vmec_file+'_patched\n$|R_{m,n}(\\theta,\phi)|$','$|R_{m,n}|$',vmec_file+'_patched_R.png')
					
					hist3d(Z_patched[:, :, 0].reshape((mpol)*(2*ntor+1)),Z_patched[:, :, 1].reshape((mpol)*(2*ntor+1)),Z_patched[:, :, 2].reshape((mpol)*(2*ntor+1)),
					vmec_file+'_patched\n$|Z_{m,n}(\\theta,\phi)|$','$|Z_{m,n}|$',vmec_file+'_patched_Z.png')
			
					# CALCULATE THE POINT-TO-POINT ERROR 
					if len(gamma)==len(gamma_patched) and len(gamma[0])==len(gamma_patched[0]):
		
						gamma_diff=numpy.subtract(gamma,gamma_patched)
						
						tmp=[]
						error=[]

						for i in range(gamma_diff.shape[0]):
							norm_diff=numpy.linalg.norm(gamma_diff[i,:,:])
							norm=numpy.linalg.norm(gamma[i,:,:])
							
							error.append(norm_diff/norm*100)

						mean_error=numpy.mean(error)
						min_errror=numpy.min(error)
						max_error=numpy.max(error)
						

						tmp.append(n_fourier)
						tmp.append(mean_error)
						tmp.append(min_errror)
						tmp.append(max_error)
						
						error_trend.append(tmp)
						
						

						print('\n#############################################\n')
						print(vmec_file+' NUMBER OF FOURIER COEFFICIENTS: '+str(n_fourier_tot)+'\n')
						print(vmec_file+' ASPECT RATIO: '+str(ar)+'\n')
						print(vmec_file+'_patched NUMBER OF FOURIER COEFFICIENTS: '+str(n_fourier)+'\n')
						print(vmec_file+'_patched LAST m: '+str(last_m)+' LAST n: '+str(last_n)+'\n')
						print(vmec_file+'_patched ASPECT RATIO: '+str(ar_patched)+'\n')
						print('MEAN ERROR: '+str(mean_error)+' %\n')
						print('MIN ERROR: '+str(min_errror)+' %\n')
						print('MAX ERROR: '+str(max_error)+' %\n')
						print('#############################################\n')

						
					error_trend_tot.append(error_trend)
					trim_n+=1
				else:
					index_out_of_range=True
					print('ALLERT! INDEX m CANNOT BE REDUCED ANYMORE!')
		
								
		fig = plt.figure()
		ax = plt.axes()

		nf=[nf[0] for nf in error_trend_tot[len(error_trend_tot)-1]]
		mean_e=[mean_e[1] for mean_e in error_trend_tot[len(error_trend_tot)-1]]
		min_e=[min_e[2] for min_e in error_trend_tot[len(error_trend_tot)-1]]
		max_e=[max_e[3] for max_e in error_trend_tot[len(error_trend_tot)-1]]

		ax.plot(nf,mean_e)
		ax.plot(nf,min_e)
		ax.plot(nf,max_e)

		ax.set_title(vmec_file+' - '+vmec_file+'_patched\nError trend')
		ax.set_xlabel('Number of Fourier coefficients')
		ax.set_ylabel('Error [%]')

		ax.set_xlim(1,n_fourier_tot)
		ax.set_ylim(0,min_mean_error)
		
		ax.set_xticks(numpy.arange(0, n_fourier_tot, 20))
		ax.set_yticks(numpy.arange(0, min_mean_error+0.005, 0.005))
		
		ax.legend(['Mean error', 'Min error','Max error'])

		plt.savefig(vmec_file+'_error.png',bbox_inches='tight')
		plt.close()


		plot_boundaries(gamma[:, :, 0], gamma[:, :, 1], gamma[:, :, 2],
		gamma_patched[:, :, 0], gamma_patched[:, :, 1], gamma_patched[:, :, 2],
		[-1,1],[-1,1],[-1,1],vmec_file+'_backward_selection.png')
	
		plot_boundaries(gamma[:, :, 0], gamma[:, :, 1], gamma[:, :, 2],
		gamma_patched[:, :, 0], gamma_patched[:, :, 1], gamma_patched[:, :, 2],
		[-1,0],[-1,1],[-1,1],vmec_file+'_backward_selection_zoom1.png')
		
		plot_boundaries(gamma[:, :, 0], gamma[:, :, 1], gamma[:, :, 2],
		gamma_patched[:, :, 0], gamma_patched[:, :, 1], gamma_patched[:, :, 2],
		[0,1],[-1,1],[-1,1],vmec_file+'_backward_selection_zoom2.png')

nf_list=[]
mean_e_list=[]
min_e_list=[]
max_e_list=[]

for i in range(n_fourier_tot,1,-1):
	aa=0
	bb=0
	cc=0
	n_occ=0
	
	for j in range(len(error_trend_tot)):
	
		nf=[nf[0] for nf in error_trend_tot[j]]
		mean_e=[mean_e[1] for mean_e in error_trend_tot[j]]
		min_e=[min_e[2] for min_e in error_trend_tot[j]]
		max_e=[max_e[3] for max_e in error_trend_tot[j]]

		if nf.count(i)>0:
			aa+=mean_e[nf.index(i)]
			bb+=min_e[nf.index(i)]
			cc+=max_e[nf.index(i)]
			n_occ+=1
			
	if aa>0:
		nf_list.append(i)
		mean_e_list.append(aa/n_occ)
		min_e_list.append(bb/n_occ)
		max_e_list.append(cc/n_occ)


fig = plt.figure()
ax = plt.axes()

ax.plot(nf_list,mean_e_list)
ax.plot(nf_list,min_e_list)
ax.plot(nf_list,max_e_list)

ax.set_title('Mean Error trend')
ax.set_xlabel('Number of Fourier coefficients')
ax.set_ylabel('Mean Error [%]')

ax.set_xlim(1,n_fourier_tot)
ax.set_ylim(0,min_mean_error)
		
ax.set_xticks(numpy.arange(0, n_fourier_tot, 20))
ax.set_yticks(numpy.arange(0, min_mean_error+0.005, 0.005))
		
ax.legend(['Mean error', 'Min error','Max error'])

plt.savefig('mean_error.png',bbox_inches='tight')
plt.close()
###################################################################		


