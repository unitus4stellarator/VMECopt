""" 
University of Tuscia & Renaissance Fusion

Author:	Matteo Notazio, University of Tuscia 
"""
from qsc import Qsc
import numpy as np
import random as rnd
import os

#########################################################
#QSC INPUT PARAMETERS (https://landreman.github.io/pyQSC/usage.html#input-parameters)
#########################################################
#rc:     the cosine components of the axis radial component (MINIMAL PARAMETER)
#zs:     the sine components of the axis vertical component (MINIMAL PARAMETER)
#rs:     the sine components of the axis radial component
#zc:     the cosine components of the axis vertical component
#nfp:    the number of field periods (MINIMAL PARAMETER)
#etabar: a scalar that specifies the strength of the first order magnetic field modulus (MINIMAL PARAMETER)
#sigma0: the value of the function sigma at phi=0, which is 0 for stellarator-symmetry
#B0:     the strength of the magnetic field on-axis
#I2:     the second derivative of the current with respect to the radial variable r
#sG:     sign of the Boozer function G
#spsi:   sign of the toroidal flux function Psi
#nphi:   toroidal resolution specifying the number of points in a grid along the axis
#B2s:    a scalar that specifies the strength of the sine component of the second order magnetic field modulus, 0 for stellarator-symmetry
#B2c:    a scalar that specifies the strength of the cosine component of the second order magnetic field modulus
#p2:     the second derivative of the pressure with respect to the radial variable r, usually negative
#order:  a string that specifies the order of the expansion, “r1”, “r2” or “r3”. For “r3” only the X3c1, Y3c1 and Y3s1 components are calculated
#########################################################

numfieldperiods=3
etabar=-0.9

numconfigs=10
rz_length=3

if not os.path.exists(os.getcwd()+'/vmec_files'):
	os.mkdir(os.getcwd()+'/vmec_files')

os.chdir(os.getcwd()+'/vmec_files')


for i in range(numconfigs):
	rcos=np.zeros(rz_length)
	zsin=np.zeros(rz_length)
	
	rcos[0]=1
	zsin[0]=0
	
	filename='input.stellqsc'+str(i+1)
	
	print(f'##################\nGENERATING {filename}...\n')
	
	for j in range(1,rz_length):
		val=rnd.uniform(0,0.1)
		rcos[j]=val
		zsin[j]=-val

	print(rcos)
	print(zsin)
	
	stel=Qsc(rc=rcos,zs=zsin,nfp=numfieldperiods,etabar=etabar)
	stel.to_vmec(filename)
	
	print(f'\n{filename} was succesfully generated!')


