# -*- coding: utf-8 -*-
#******************公共默认参数，无需修改******************
import os
import scipy.io as sio

predir = './pre-simulator/' + refname + '/'
postdir = './post-simulator/' + refname + '/'

#########未用到的参数
cudarecv_handoff = 'handoff_params_usrp6'  # Simulated
ip_list = [6]  # Simulated
weekno = 2258  # If running simulated data
#########未用到的参数

dir_req = [datpath, predir, postdir]
for d in dir_req:
    if not os.path.exists(d):
        os.makedirs(d)

init_time = start_time
run_time = proc_time

prefix = 'scalar_skip%ds_proc%ds' % (init_time, run_time) + '/'
prepath = predir + prefix
postfix = 'dpe_skip%ds_proc%.2fs_' % (DPE_start_time, DPE_run_time) + descriptor + '/'
postpath = postdir + prefix + postfix
try:
    os.makedirs(postpath)
except:
    if os.listdir(postpath):
        print('Warning:', postpath, 'not an empty directory.')

# 标量跟踪结果保存路径
dir_acq = prepath + 'acq'
dir_trk = prepath + 'trk'
dir_eph = prepath + 'eph'
dir_scalar = [dir_acq, dir_trk, dir_eph]
for dir in dir_scalar:
    if not os.path.exists(dir):
        os.makedirs(dir)

# 保存运行时间相关参数
save_dic = {'start_time': start_time, 'proc_time': proc_time,
            'DPE_start_time': DPE_start_time, 'DPE_run_time': DPE_run_time,
            'DPE_interval': DPE_interval, 'DPE_corr_save_interval': DPE_corr_save_interval}
sio.savemat(os.path.join(postpath, 'dpe_runtime.mat'), save_dic)