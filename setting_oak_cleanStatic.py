# -*- coding: utf-8 -*-
# oak静态数据
import numpy as np
import os

refname    = 'oak_cleanStatic'  # Simulated
filename = 'cleanStatic_gps_oak.bin'
datname    = '2023'
descriptor = 'test0908_1'
fs = 5e6
fi = 0.0e6
datatype = np.dtype([('i', np.int16), ('q', np.int16)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 50
proc_time     = 40
max_lead_time = 0
DPE_run_time = 10

acq_only      = False
load_acq = True
load_trk = True
# prn_list = [8, 10, 11, 12, 14, 15, 20, 21, 24, 25, 27, 31, 32] # Simulated
# prn_list = [8, 10, 15, 20, 21, 32] # Simulated
# prn_list = [15] # Simulated
prn_list = range(1,33) # Simulated

datpath  = 'D:/academic/DPEdata/'
predir   = './pre-simulator/' + refname + '/'
postdir  = './post-simulator/' + refname + '/'



### No need to change: house keeping ###
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

acq_dir = prepath + 'acq'
if not os.path.exists(acq_dir):
    os.makedirs(acq_dir)