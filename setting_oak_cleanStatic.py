# -*- coding: utf-8 -*-
# oak静态数据
import numpy as np

refname    = 'cleanStatic_gps_oak'  # Simulated
filename = 'cleanStatic_gps_oak.bin'
datname    = '2023'
descriptor = 'test3'
fs = 5e6
fi = 0.0e6
datatype = np.dtype([('i', np.int16), ('q', np.int16)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 30
proc_time     = 60
max_lead_time = 0

acq_only      = False
# prn_list = [8, 10, 11, 12, 14, 15, 20, 21, 24, 25, 27, 31, 32] # Simulated
prn_list = [8, 10, 15, 20, 21, 32] # Simulated
# prn_list = [15] # Simulated

datpath  = 'D:/academic/DPEdata/'
predir   = './pre-simulator/oak_cleanStatic/'
postdir  = './post-simulator/oak_cleanStatic/'



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