# -*- coding: utf-8 -*-
# execfile('setting_finland.py')
execfile('setting_oak_cleanStatic.py')

### Main code starts
from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos,ephemeris
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver
import printer
import threading

import numpy as np
import scipy.io as sio
import time,os

# 实例化接收机
scalar_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                               abspath=datpath + filename, # 数据所在路径
                                               fs=fs, fi=fi, ds=1.0, # 采样率、中频频率、比例系数（不用管）
                                               datatype=datatype, # 单个采样点的数据类型
                                               notes='Data set ' + refname),
                                               mcount_max=run_time * 1000 + 1) # 总的历元数

first_dir = 'end-of-3sec' # 前三秒标量跟踪结果保存路径
second_dir = 'end-of-%dsec' % proc_time # 所有时间，即proc_time标量跟踪结果保存路径

# 捕获
if not load_acq:
    print('Start scalar acquisition @ %ds' % init_time)
    scalar_sdr.add_channels(prn_list) # 初始化各通道
    scalar_sdr.rawfile.seek_rawfile(init_time * 1000 * scalar_sdr.rawfile.S) # 跳过init_time对应的时间
    acq_results = scalar_sdr.scalar_acquisition(prn_list) # 捕获
    del_acq_list = []
    for prn_idx, prn in enumerate(prn_list):
        if not acq_results[prn_idx][0]:
            del_acq_list += [prn]
    scalar_sdr.del_channels(del_acq_list) # 如果acq_results对应卫星捕获结果为False，则删除该通道，即跟踪阶段不再跟踪该卫星
    save_dic = {'acq_results': acq_results, 'prn_list': prn_list, 'init_time': init_time}
    sio.savemat(os.path.join(acq_dir, 'acq.mat'), save_dic) # 保存捕获结果
else:
    scalar_sdr.load_acq_results(acq_dir) # 载入现有捕获结果
print ('Scalar acquisition completed.')

# 跟踪
if not load_trk:
    print('Start scalar tracking @ %ds for %ds' % (init_time, run_time))
    scalar_sdr.scalar_track(mtrack=3000) # 先跟踪3秒然后储存结果作为后续DPE处理的起点
    scalar_sdr.save_measurement_logs(dirname = prepath,subdir= first_dir) # 储存前3秒的跟踪结果
    scalar_sdr.scalar_track(mtrack=run_time * 1000 - 3000) # 完成剩余时间的跟踪
    scalar_sdr.save_measurement_logs(dirname = prepath,subdir= second_dir) # 储存run_time时间的跟踪结果
else:
    scalar_sdr.load_measurement_logs(dirname = prepath,subdir= second_dir) # 载入现有跟踪结果

# 解算星历
dir_eph = prepath + 'eph'
if not os.path.exists(dir_eph):
    os.makedirs(dir_eph)
for prn in scalar_sdr.channels:
    try:
        scalar_sdr.parse_ephemerides(prn_list=[prn], m_start=40)
        scalar_sdr.channels[prn].ephemerides.save_ephemerides(dir_eph + '/channel%d.mat' % prn, dir_eph + '/channel%d.csv' % prn)
    except:
        pass
print ('Scalar tracking completed.  Launching DP.')

# --------------------------------------------------------DPE-----------------------------------------------------------

# 初始化DPE接收机
DPE_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                            abspath=datpath + filename,
                                            fs=fs, fi=fi, ds=1.0,
                                            datatype=datatype,
                                            notes='Data set ' + refname),
                                            mcount_max = DPE_run_time * 50 + 5000)

DPE_sdr.load_measurement_logs(dirname = prepath, subdir= first_dir)

del_clist = []
for prn in DPE_sdr.channels:
    try:
        DPE_sdr.channels[prn].ephemerides = scalar_sdr.channels[prn].ephemerides
    except:
        del_clist += [prn]
DPE_sdr.del_channels(del_clist)
print ('DP Channels')
print (DPE_sdr.channels.keys())

DPE_sdr.rawfile.set_rawsnippet_settings(T=0.020, T_big=0.020)
DPE_sdr.init_dp()
print ('Init at', utils.ECEF_to_LLA(DPE_sdr.ekf.X_ECEF))

counter_max = int(DPE_run_time / DPE_sdr.rawfile.T_big)
X_list = []
rxTime_list = []
csvfile = open(postpath+'usrp.csv','w')

print ('DP Launched')
start = time.time()

printer.header(csvfile)
DPE_sdr.counter = 0
DPE_sdr.initGridInfo(counter_max)
for mc in range(counter_max):
    DPE_sdr.dp_track(1)
    DPE_sdr.counter += 1
    printer.printer(DPE_sdr.counter,weekno,DPE_sdr.rxTime_a,
                    DPE_sdr.ekf.X_ECEF,csvfile)
    X_list += [DPE_sdr.ekf.X_ECEF.copy()]
    rxTime_list += [DPE_sdr.rxTime_a]
    if DPE_sdr.counter % 100 == 0:
        np.save(postpath + 'usrp_X', np.array(X_list))
        np.save(postpath + 'usrp_t', np.array(rxTime_list))
        DPE_sdr.save_measurement_logs(dirname=postpath, subdir='end-of-dp')
        print ('DP running')
        print ('Current time: %d/%d'%(DPE_sdr.counter, counter_max))
        print ('DP File saved, continue running.')

elapse = time.time() - start
print ('DP success!')
print ('%.2f seconds elapsed for data proc.'%elapse)
np.save(postpath + 'usrp_X', np.array(X_list))
np.save(postpath + 'usrp_t', np.array(rxTime_list))
DPE_sdr.save_measurement_logs(dirname=postpath, subdir='end-of-dp')
csvfile.close()


