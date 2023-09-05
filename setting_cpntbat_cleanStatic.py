# -*- coding: utf-8 -*-
# setting_cpntbat_cleanStatic
import numpy as np

refname    = '210A_cleanStatic_test1'  # Simulated
filename = '210A_cleanStatic_test1.bin'
datname    = '2023'
descriptor = 'test1'
fs = 100e6
fi = 62.58e6
datatype = np.dtype([('i', np.int8)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 10
proc_time     = 40
max_lead_time = 0

acq_only      = False
prn_list = [4, 16, 22, 26, 27, 31] # Simulated

datpath  = 'D:/academic/DPEdata/'
predir   = './pre-simulator/CPNTBAT_cleanStatic/'
postdir  = './post-simulator/CPNTBAT_cleanStatic/'



### No need to change: house keeping ###

import os
dir_req  = [datpath,predir,postdir]#,'netcsv','html']
for d in dir_req:
    if not os.path.exists(d):
        os.makedirs(d)

### More housekeeping...
init_time = start_time
run_time = proc_time

prefix   = descriptor    + '_skip%ds_proc%ds/'%(init_time,run_time)
prepath  = predir     + prefix
postfix  = descriptor + '_skip%ds_proc%ds/'%(init_time,run_time)
postpath = postdir    + postfix
try:
    os.makedirs(postpath)
except:
    if os.listdir(postpath):
        print('Warning:',postpath,'not an empty directory.')