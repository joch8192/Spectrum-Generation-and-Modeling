#!/home/data/joch8192/anaconda3/bin/python
# Generate UV/Vis spectra from QM excitations
import os
##############
# Parameters #
##############
# Create list of files to iterate over
#os.system("ls frame.* > filelist.txt")
os.system("ls frame.ricc2.* > filelist.txt")
# Location of file containing full paths to QM output files (one path per line, no empty lines)
filelist = 'filelist.txt'

# Output settings:
outpath = './'
outfile = 'ret1' # prefix for output files

# range (nm) of generated spectra
lambdamin = 200
lambdamax = 800

##################
# Load Libraries #
##################
import sys                      # for error exit
import os                       # system commands like mkdir
import math as m                # for mathematical operations
import matplotlib.pyplot as plt # for plotting
import numpy as np              # nice arrays
import re                       # for splitting with multiple delimters
import subprocess               # for reading remote file via ssh
from itertools import islice    # to print n lines after match
from operator import itemgetter # extract item from list

######################
# Physical constants #
######################
# according to NIST (https://physics.nist.gov/cuu/Constants/index.html)
c  = 29979245800.0            # speed of light in vacuum (cm/s)
e  = 4.80320425*pow(10,-10)   # elementary charge (esu)
e2 = 1.6021766208*pow(10,-19) # elementary charge (C)
h  = 6.626070040*pow(10,-34)  # Planck constant (Js)
me = 9.10938356*pow(10,-28)   # electron mass (g)
N  = 6.022140857*pow(10,23)   # Avogadro number (1/mol)

#############
# Functions #
#############

# Convert nm to eV and vice-versa
def nm_ev(x):
    return (h*c*pow(10,7)) / (x*e2)

###########################
# Extract data from files #
###########################

# Read output file names
files = []
with open(filelist, 'r') as Eingang:
    for line in Eingang:
        files.append(line.rstrip('\n'))

# Determine file format of output files
# Currently tested: ORCA      --> CIS, sTD-DFT, TD-DFT
#                   Turbomole --> ADC(2), CC2
#                   Gaussian  --> TD
software = open(files[0]).read()
if software.find('R I C C 2') != -1:
    software = 'turbomole'
elif software.find('O   R   C   A') != -1:
    software = 'orca'
elif software.find(' ******************************************\n Gaussian') != -1:
    software = 'gaussian'
elif software.find('Dalton - An Electronic Structure Program') != -1:
    software = 'dalton'
else:
    sys.exit("File format not recognised.")

# Extract wavelength and oscillator strength
data    = [] # holder for wavelengths/oscillator strengths
if software == 'orca':
    # Find out number of excited states: grenze
    with open(files[0]) as Eingang:
        for line in Eingang:
            if 'roots found' in line:
                grenze = 30 # use first 30 states for sTD
                break
            elif 'Number of roots to be determined' in line:
                grenze = int(re.split(' |\n',line)[-2])
                break
    Eingang.close()
    
    # Loop over files
    for i in range(len(files)):
        tmpdata = [] # temporary holder for wavelengths/oscillator strengths
        # Search file for excited states
        with open(files[i]) as Eingang:
            for line in Eingang:
                # If matches string, extract data
                if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                    Eingang.readline()
                    Eingang.readline()
                    Eingang.readline()
                    Eingang.readline()
                    # Extract 'grenze' excitations
                    for k in range(grenze):
                        tmp = list(filter(None,re.split(' +',list(islice(Eingang,1))[0])))
                        tmpdata.append(np.array(itemgetter(2,3)(tmp),dtype=np.float)) #itemgetter(2) is wavelength
                    break
        Eingang.close()
        data.append(tmpdata)
elif software == 'turbomole':
    # Find out number of excited states: grenze
    grenze = open(files[0]).read().count('Transition')
    
    # Loop over files
    for i in range(len(files)):
        tmpdata = [] # temporary holder for wavelengths/oscillator strengths
        wache   =  0 # number of currently found excitations
        # Search file for excited states
        with open(files[i]) as Eingang:
            for line in Eingang:
                # If matches string, extract data
                if 'Transition' in line:
                    tmp = [] # temporary holder for excitations
                    tmp.append(list(filter(None,re.split(' +',list(islice(Eingang,2))[1])))[5])
                    tmp.append(list(filter(None,re.split(' +|\n',list(islice(Eingang,11))[10])))[5])
                    tmp2 = np.array(itemgetter(0,1)(tmp),dtype=np.float)
                    print(tmp)
                    tmp2[0] = nm_ev(tmp2[0]) # convert from eV to nm
                    tmpdata.append(tmp2)
                    wache += 1
                    # stop searching once all excitations extracted
                    if wache == grenze:
                        break
        Eingang.close()
        data.append(tmpdata)
elif software == 'gaussian':
    # Find out number of excited states: grenze
    grenze = open(files[0]).read().count(' Excited State ')
    
    # Loop over files
    for i in range(len(files)):
        tmpdata = [] # temporary holder for wavelengths/oscillator strengths
        wache   =  0 # number of currently found excitations
        # Search file for excited states
        with open(files[i]) as Eingang:
            for line in Eingang:
                # If matches string, extract data
                if ' Excited State ' in line:
                    tmp = list(filter(None,re.split(' |=',line))) # current wavelength/oscillator strength
                    tmpdata.append(np.array(itemgetter(6,9)(tmp),dtype=np.float))
                    wache += 1
                    # stop searching once all excitations extracted
                    if wache == grenze:
                        break
        Eingang.close()
        data.append(tmpdata)
elif software == 'dalton':
    # Find out number of excited states: grenze
    grenze = open(files[0]).read().count(' Excited state ')
    # Loop over files
    for i in range(len(files)):
        tmpdata = [] # temporary holder for wavelengths/oscillator strengths
        wache   =  0 # number of currently found excitations
        # Search file for excited states
        # If matches string, extract data
        os.system('grep -A'+str(grenze)+' "Singlet transition operator label: XDIPLEN" '+str(files[i])+' | grep "STATE NO" > out.log')
        x=np.loadtxt('out.log', dtype=str)
        os.system('grep -A'+str(grenze)+' "Singlet transition operator label: YDIPLEN" '+str(files[i])+' | grep "STATE NO" > out.log')
        y=np.loadtxt('out.log', dtype=str)
        os.system('grep -A'+str(grenze)+' "Singlet transition operator label: ZDIPLEN" '+str(files[i])+' | grep "STATE NO" > out.log')
        z=np.loadtxt('out.log', dtype=str)
        xdip=x[:,6].astype(np.float)
        ydip=y[:,6].astype(np.float)
        zdip=z[:,6].astype(np.float)
        xe=x[:,8].astype(np.float)*0.0367493
        final=((np.power((xdip[:]),2)+np.power((ydip[:]),2)+np.power((zdip[:]),2))*2.0/3.0*xe[:])               
        while (wache != grenze):
            nm=nm_ev(x[wache,8].astype(np.float))
            tmpdata.append(np.array([nm, final[wache]],dtype=np.float))
            wache += 1
        Eingang.close()
        data.append(tmpdata)

####################################
# Generate prefactors for plotting #
####################################
# Based on Gaussian UV/Vis website (http://gaussian.com/uvvisplot)
# Uses wavelength (nm) and oscillator strength

#test
sigma  = (0.15*e2)/(h*c)
sigma2 = pow(10,-7)*sigma

# Prefactor from equation 5, 1.3062974*pow(10,8) in example
eq5 = m.sqrt(np.pi)*m.pow(e,2)*N / (1000*m.log(10)*m.pow(c,2)*me)

# Prefactors from equation 6 and 9
epsi = [] # holder for prefactors
for j in range(len(data)):
    tmp = []
    for k in range(len(data[j])):
        tmp.append(eq5*data[j][k][1]/sigma)
    epsi.append(np.array(tmp,dtype=np.float))

####################
# Generate spectra #
####################

# Wavelength (nm), interval from 'lambdamin' to 'lambdamax'
wavelength = np.linspace(lambdamin, lambdamax, lambdamax-lambdamin+1)

# Absorbance (L * mol^-1 * cm^-1)
absorbance = [] # holder for absorbance
for j in range(len(data)):
    tmpe=[]
    for n in range(len(wavelength)):
        # Equation 9
        tmp = 0.0
        for k in range(len(data[j])):
            tmp = tmp + epsi[j][k]*m.exp(-m.pow(((1/wavelength[n])-(1/data[j][k][0]))/sigma2,2))
        tmpe.append(tmp)
    absorbance.append(np.array(tmpe,dtype=np.float))

# Generate average spectrum
fritz = np.mean(absorbance,axis=0)

#########################
# Save spectra as array #
#########################

np.savetxt('{}/{}.log'.format(outpath,outfile),np.column_stack((wavelength, fritz))) # wavelength array (x axis)
np.save('{}/{}'.format(outpath,outfile),absorbance)        # absorbance arrays for individual snapshots (y axis)

exit
