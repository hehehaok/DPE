# -*- coding: utf-8 -*-
execfile('paramSettings\\setting_texbat_ds4.py')
# setting_finland.py
# setting_texbat_ds4.py

### Main code starts
from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile, utils, satpos, ephemeris
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter, naveng
from pythonreceiver import receiver
import printer
import threading
import h5py
import numpy as np
import scipy.io as sio
import time
import os

# 初始化标量跟踪接收机
scalar_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                               abspath=datpath + filename,  # 数据所在路径
                                               fs=fs, fi=fi, ds=1.0,  # 采样率、中频频率、比例系数（不用管）
                                               datatype=datatype,  # 单个采样点的数据类型
                                               notes='Data set ' + refname),
                                               mcount_max=run_time * 1000 + 1)  # 总的历元数，每个历元0.001秒，即一个伪码周期

# 捕获
if not load_acq:
    print('Start scalar acquisition @ %ds' % init_time)
    scalar_sdr.add_channels(prn_list)  # 初始化各通道
    scalar_sdr.rawfile.seek_rawfile(init_time * 1000 * scalar_sdr.rawfile.S)  # 跳过init_time对应的时间
    acq_results = scalar_sdr.scalar_acquisition(prn_list)  # 捕获
    del_acq_list = []
    for prn_idx, prn in enumerate(prn_list):
        if not acq_results[prn_idx][0]:
            del_acq_list += [prn]
    scalar_sdr.del_channels(del_acq_list)  # 如果acq_results对应卫星捕获结果为False，则删除该通道，即跟踪阶段不再跟踪该卫星
    save_dic = {'acq_results': acq_results, 'prn_list': prn_list, 'init_time': init_time}
    sio.savemat(os.path.join(dir_acq, 'acq.mat'), save_dic)  # 保存捕获结果
else:
    scalar_sdr.load_acq_results(dirname=dir_acq)  # 载入现有捕获结果
print ('Scalar acquisition completed.')

# 跟踪
if not load_trk:
    print('Start scalar tracking @ %ds for %ds' % (init_time, run_time))
    scalar_sdr.scalar_track(mtrack=run_time * 1000)  # 跟踪
    scalar_sdr.save_trk_results(dirname=dir_trk)  # 储存run_time时间的跟踪结果
else:
    scalar_sdr.load_trk_results(dirname=dir_trk)  # 载入现有跟踪结果

# 解算星历
for prn in scalar_sdr.channels:
    try:
        scalar_sdr.parse_ephemerides(prn_list=[prn], m_start=40)  # 解码获得星历参数
        scalar_sdr.channels[prn].ephemerides.save_ephemerides(dir_eph + '/channel%d.mat' % prn, dir_eph + '/channel%d.csv' % prn)
    except:
        pass
print ('Scalar tracking completed.  Launching DPE.')  # 标量跟踪结束

# --------------------------------------------------------DPE-----------------------------------------------------------

# 初始化DPE接收机
DPE_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                            abspath=datpath + filename,
                                            fs=fs, fi=fi, ds=1.0,
                                            datatype=datatype,
                                            notes='Data set ' + refname),
                                            mcount_max=int(DPE_run_time*50)+5000)  # 总的历元数，每个历元0.02秒，即20个伪码周期

DPE_sdr.load_trk_results(dirname=dir_trk, DPE_lead_time=DPE_lead_time)

del_clist = []
for prn in DPE_sdr.channels:
    try:
        DPE_sdr.channels[prn].ephemerides = scalar_sdr.channels[prn].ephemerides  # 载入星历参数用于解算DPE初始位置
    except:
        del_clist += [prn]
DPE_sdr.del_channels(del_clist)  # 仅有包含星历的通道被保留，若跟踪后未能解算出星历则删除对应通道
print ('DP Channels')
print (DPE_sdr.channels.keys())

DPE_sdr.rawfile.set_rawsnippet_settings(T=0.020, T_big=DPE_interval)  # 每次读取数据设置为0.02秒，即一个解算历元为0.02秒
DPE_sdr.init_dp()  # 标量跟踪结果进行初始化
print ('Init at', utils.ECEF_to_LLA(DPE_sdr.ekf.X_ECEF))

counter_max = int(DPE_run_time / DPE_sdr.rawfile.T_big)  # DPE总的历元数
corr_interval = int(DPE_corr_save_interval / DPE_interval)  # DPE结果保存间隔 单位：历元数

X_list = []
rxTime_list = []
csvfile = open(postpath+'usrp.csv','w')

print ('DP Launched')
start = time.time()

printer.header(csvfile)  # 在表格中打印表头(即第一行的各列标题)
DPE_sdr.initPlanInfo(counter_max, corr_interval, dpe_plan, grid_param, ars_param)  # 初始化储存DPE相关结果的矩阵以及网格
for mc in range(counter_max):
    DPE_sdr.dp_track(1)  # DPE解算，每次只处理一个历元
    DPE_sdr.counter += 1
    printer.printer(DPE_sdr.counter,weekno,DPE_sdr.rxTime_a,
                    DPE_sdr.ekf.X_ECEF,csvfile)  # 打印表内容，即打印完表头后，在对应表头下方打印每个历元的数据
    X_list += [DPE_sdr.ekf.X_ECEF.copy()]  # 记录位置信息
    rxTime_list += [DPE_sdr.rxTime_a]  # 记录接收机本地时间
    if DPE_sdr.counter % 10 == 0:  # 每跑完20ms*10=0.2s储存一次结果，这样就算中断也有一部分结果
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
DPE_sdr.save_measurement_logs(dirname=postpath, subdir='end-of-dp')  # 最后再完整的存储一次数据
csvfile.close()


