# -*- coding: utf-8 -*-
# texbat_ds4
import numpy as np
import os
import scipy.io as sio

refname    = 'texbat_ds4'  # Simulated
filename = 'texbat_ds4.bin'
datname    = '2023'
descriptor = 'test2'
fs = 25e6
fi = 0.0e6
datatype = np.dtype([('i', np.int16), ('q', np.int16)])

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6' # Simulated
ip_list       = [6] # Simulated
weekno        = 2258 # If running simulated data
#########未用到的参数

start_time    = 90
proc_time     = 40
DPE_start_time = 100
DPE_lead_time = DPE_start_time - start_time  # must make 3<=DPE_lead_time<=proc_time
DPE_run_time = 50
DPE_interval = 0.02  # 进行DPE的间隔
DPE_corr_save_interval = 1  # 保存DPE流形结果图的间隔

acq_only      = False
load_acq = False
load_trk = False
# prn_list = [8, 10, 11, 12, 14, 15, 20, 21, 24, 25, 27, 31, 32] # Simulated
prn_list = range(1, 33)  # Simulated
# prn_list = [15] # Simulated

datpath  = 'D:/DATA/TEXBAT/'
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

# 保存运行时间相关参数
save_dic = {'start_time': start_time, 'proc_time': proc_time,
            'DPE_start_time': DPE_start_time, 'DPE_run_time': DPE_run_time,
            'DPE_interval': DPE_interval, 'DPE_corr_save_interval': DPE_corr_save_interval}
sio.savemat(os.path.join(postpath, 'dpe_runtime.mat'), save_dic)