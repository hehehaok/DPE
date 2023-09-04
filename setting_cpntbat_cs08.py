# -*- coding: utf-8 -*-
# oak静态数据
import numpy as np

refname    = 'cleanStatic_gps_cpntbat_cs08'  # Simulated
filename = '08_8bit_100Msps.dat'
datname    = '2023'
descriptor = 'test2'
fs = 100e6
fi = 62.58e6
datatype = np.dtype([('i', np.int8)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 130
proc_time     = 40
max_lead_time = 0

acq_only      = False
# prn_list = [8, 10, 11, 12, 14, 15, 20, 21, 24, 25, 27, 31, 32] # Simulated
prn_list = range(1,33) # Simulated
# prn_list = [15] # Simulated

datpath  = 'D:/DATA/CPNTBAT/'
predir   = './pre-simulator/cpntbat_cs08/'
postdir  = './post-simulator/cpntbat_cs08/'



### No need to change: house keeping ###

import os
dir_req  = [datpath,predir,postdir]#,'netcsv','html']
for d in dir_req:
    if not os.path.exists(d):
        os.makedirs(d)

### More housekeeping...
init_time = start_time
run_time = proc_time

prefix   = datname    + '_skip%d_start%d/'%(init_time,start_time)
prepath  = predir     + prefix
postfix  = descriptor + 'proc%ds/'%run_time
postpath = postdir    + postfix
try:
    os.makedirs(postpath)
except:
    if os.listdir(postpath):
        print('Warning:',postpath,'not an empty directory.')