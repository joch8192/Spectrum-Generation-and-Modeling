#!/home/data/joch8192/anaconda3/bin/python

#exit()

#############
# Libraries #
#############
import os     # system commands
import socket # get hostname
import math as m
import numpy as np
import re #for splitting with multiple delimters
import subprocess #to read remote file via ssh
from itertools import islice    # to print following n lines after match
from scipy.stats import norm
from operator import itemgetter # extract item from list
import scipy.signal as signal

##############
# Atom names #
##############

single = []
double = []

# QM, prot.
c5  = '18-C'
c6  = '19-C'
c7  = '20-C'
c8  = '21-C'
c9  = '22-C'
c10 = '23-C'
c11 = '24-C'
c12 = '25-C'
c13 = '26-C'
c14 = '27-C'
c15 = '28-C'
n16 = '12-N'
#atoms = [c15,c14,c13,c12,c11,c10,c9,c8,c7,c6,c5]
atoms = [c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15]
tmps   = []
tmpd   = []
for i in range(0,len(atoms)-1,2):
    tmps.append('{} , {}'.format(atoms[i],atoms[i+1]))
    tmpd.append('{} , {}'.format(atoms[i+1],atoms[i+2]))
single.append(tmps)
double.append(tmpd)

###########
# Get BOA #
###########
username = 'joch8192'
hostname = 'landau'
system   = ['4kly_pD97']
qmregion = ['std']

# data holders
data = [] # bond orders
x    = [] # x axis for plotting

# check if remote host
if hostname in socket.getfqdn():
    #locally stored
    for i in range(len(system)):
        datatmp_i = []
        print(system[i])
        for j in range(len(qmregion)):
            datatmp_j = []
            root = '/home/data/joch8192/TDDFT_new_TDA_off/4knf/results'  ## change based on your path                                                       
            # open files
            for k in range(1,100,1):
                tmp = []
                with open('{}/frame.{}'.format(root,k)) as Eingang:
                   Eingang.readline() # skip titles line
                   for line in Eingang:
                      if 'B(' in line:
                        tmp.append(list(filter(None,re.split('B\( |\n',line.replace(' ) ', '')))))
                   Eingang.close()
                datatmp_j.append(np.concatenate(tmp))
            data.append(datatmp_j)
boa = []
# print (data, len(data))
for i in range(len(data)): # systems
#    print (len(data[i]))
    tmp_i = []
    for j in range(len(data[i])): # QM regions
#        print (len(data[i][j]))
        tmp_j = []
#           print (len(data[i][j][k]))
#            print ((data[i][j][k]))
        av_single = []
        av_double = []
        N = 0
#            print (single[N])
#            print (double[N])
        for l in range(len(single[N])):
#                print (data[i][j])
#                print (single[N][l])
          av_single.append(float(list(filter(lambda x: re.search(single[N][l], x), data[i][j]))[0][-7:-1]))
          av_double.append(float(list(filter(lambda x: re.search(double[N][l], x), data[i][j]))[0][-7:-1]))
        tmp_j.append(np.average(av_single) - np.average(av_double))
        print(tmp_j[0])
        np.savetxt('/home/data/joch8192/boa_'+str(system[i]),tmp_j) ## change based on your path   
    boa.append(tmp_j)           


    
np.savetxt('/home/data/joch8192/boa_Dat',boa) ## change based on your path   

exit()

