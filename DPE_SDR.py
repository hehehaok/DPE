# -*- coding: utf-8 -*-
execfile('setting_finland.py')

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

# 初始化接收机
scalar_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                               abspath=datpath + filename,
                                               fs=fs, fi=fi, ds=1.0,
                                               datatype=datatype,
                                               notes='Data set ' + refname),
                                               mcount_max=run_time * 1000 + 1)

first_dir = 'end-of-1sec'
second_dir = 'end-of-%dsec' % proc_time

# 捕获
if not load_acq:
    print('Start scalar acquisition @ %ds' % init_time)
    scalar_sdr.add_channels(prn_list)
    scalar_sdr.rawfile.seek_rawfile(init_time * 1000 * scalar_sdr.rawfile.S)
    acq_results = scalar_sdr.scalar_acquisition(prn_list)
    if save_acq:
        save_dic = {'acq_results': acq_results, 'prn_list': prn_list}
        acq_dir = prepath + 'acq'
        if not os.path.exists(acq_dir):
            os.makedirs(acq_dir)
        sio.savemat(os.path.join(acq_dir, 'acq.mat'), save_dic)
else:
    acq_dir = prepath + 'acq'
    load_dic = sio.loadmat(os.path.join(acq_dir, 'acq.mat'))
    acq_results = load_dic['acq_results']
    prn_list = load_dic['prn_list']
    scalar_sdr.add_channels(prn_list[0])
    for prn_idx, prn in enumerate(prn_list[0]):
        found, rc, ri, fc, fi, cppr, cppm = acq_results[prn_idx, :]
        state = 'True' if found else 'False'
        scalar_sdr.channels[prn].set_scalar_params(rc=rc, ri=ri, fc=fc, fi=fi)
        print('PRN: %d, Found: %s, Code: %.2f chips, Carrier: %.2f cycles, '
              'Doppler: %.2f Hz, Cppr: %.2f, Cppm: %.2f'
              % (prn, state, rc, ri, fi, cppr, cppm))
print ('Scalar acquisition completed.')

# 跟踪
if not load_trk:
    print('Start scalar tracking @ %ds for %ds' % (init_time, run_time))
    scalar_sdr.scalar_track(mtrack=1000)
    scalar_sdr.save_measurement_logs(dirname = prepath,subdir= first_dir)
    scalar_sdr.scalar_track(mtrack=run_time * 1000 - 1000)
    scalar_sdr.save_measurement_logs(dirname = prepath,subdir= second_dir)
else:
    scalar_sdr.load_measurement_logs(dirname = prepath,subdir= second_dir)

# 解算星历
dir = prepath + 'eph'
if not os.path.exists(dir):
    os.makedirs(dir)
for prn in scalar_sdr.channels:
    try:
        scalar_sdr.parse_ephemerides(prn_list=[prn], m_start=40)
        scalar_sdr.channels[prn].ephemerides.save_ephemerides(dir + '/channel%d.mat' % prn, dir + '/channel%d.csv' % prn)
    except:
        pass
print ('Scalar tracking completed.  Launching DP.')

# DPE

# 初始化DPE接收机
DPE_sdr = receiver.Receiver(rawfile.RawFile(metafile=None,
                                               abspath=datpath + filename,
                                               fs=fs, fi=fi, ds=1.0,
                                               datatype=datatype,
                                               notes='Data set ' + refname),
                                               mcount_max = run_time * 50 + 5000)

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

rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI = naveng.calculate_nav_soln(DPE_sdr)
DPE_sdr.rawfile.set_rawsnippet_settings(T=0.020, T_big=0.020)
DPE_sdr.init_dp()
print ('Init at', utils.ECEF_to_LLA(DPE_sdr.ekf.X_ECEF))

counter = 0
counter_max = int(run_time / DPE_sdr.rawfile.T_big)
X_list = []
rxTime_list = []
csvfile = open(postpath+'usrp.csv','w')

print ('DP Launched')
start = time.time()

printer.header(csvfile)
DPE_sdr.initGridInfo(counter_max)
for mc in range(counter_max):
    counter += 1
    DPE_sdr.dp_track(1)
    printer.printer(counter,weekno,DPE_sdr.rxTime_a,
                    DPE_sdr.ekf.X_ECEF,csvfile)
    X_list += [DPE_sdr.ekf.X_ECEF.copy()]
    rxTime_list += [DPE_sdr.rxTime_a]
    if counter % 100 == 0:
        np.save(postpath + 'usrp_X', np.array(X_list))
        np.save(postpath + 'usrp_t', np.array(rxTime_list))
        DPE_sdr.save_measurement_logs(dirname=postpath, subdir='end-of-dp')
        print ('DP running')
        print ('Current time: %d/%d'%(counter, counter_max))
        print ('DP File saved, continue running.')


elapse = time.time() - start
print ('DP success!')
print ('%.2f seconds elapsed for data proc.'%elapse)
np.save(postpath + 'usrp_X', np.array(X_list))
np.save(postpath + 'usrp_t', np.array(rxTime_list))
DPE_sdr.save_measurement_logs(dirname=postpath, subdir='end-of-dp')
csvfile.close()


