# -*- coding: utf-8 -*-
# CPNTBAT_cleanStatic
import numpy as np
import os

refname    = 'CPNTBAT_cleanStatic'  # Simulated
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
DPE_run_time = 10

acq_only      = False
load_acq = False
load_trk = False
prn_list = [4, 16, 22, 26, 27, 31] # Simulated

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