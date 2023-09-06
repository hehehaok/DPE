# -*- coding: utf-8 -*-
# 芬兰数据
import numpy as np

refname    = 'cleanStatic_gps_finland'  # Simulated
filename = 'cleanStatic_gps_finland.dat'
datname    = '2023'
descriptor = 'test0905'
fs = 26e6
fi = 6.39e6
datatype = np.dtype([('i', np.int8)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 0
proc_time     = 40
max_lead_time = 0

acq_only      = False
# prn_list = [10, 11, 13, 15, 17, 19, 20, 24, 28, 30] # Simulated
prn_list = [13, 15, 17, 24, 28] # Simulated
# prn_list = [15] # Simulated

datpath  = 'D:/academic/DPEdata/'
predir   = './pre-simulator/finland_cleanStatic/'
postdir  = './post-simulator/finland_cleanStatic/'



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