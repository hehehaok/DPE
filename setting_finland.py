# -*- coding: utf-8 -*-
# 芬兰数据
import numpy as np
import os
import scipy.io as sio

refname    = 'finland_cleanStatic'  # Simulated
filename = 'cleanStatic_gps_finland.dat'
datname    = '2023'
descriptor = 'test1012'
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
DPE_start_time = 10
DPE_lead_time = DPE_start_time - start_time  # must make 3<=DPE_lead_time<=proc_time
DPE_run_time = 1
DPE_interval = 0.02  # 进行DPE的间隔
DPE_corr_save_interval = 0.2  # 保存DPE流形结果图的间隔

acq_only      = False
load_acq = True
load_trk = True
# prn_list = [10, 11, 13, 15, 17, 19, 20, 24, 28, 30] # Simulated
prn_list = [13, 15, 17, 24, 28] # Simulated
# prn_list = [15] # Simulated

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

# 保存运行时间相关参数
save_dic = {'start_time': start_time, 'proc_time': proc_time,
            'DPE_start_time': DPE_start_time, 'DPE_run_time': DPE_run_time,
            'DPE_interval': DPE_interval, 'DPE_corr_save_interval': DPE_corr_save_interval}
sio.savemat(os.path.join(postpath, 'dpe_runtime.mat'), save_dic)