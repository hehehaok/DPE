# -*- coding: utf-8 -*-
execfile('setting_cpntbat_cleanStatic.py')

from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

#import printer,os
#from pythonreceiver.libgnss import pygmaps
from pythonreceiver.vector.ekf import ExtendedKalmanFilter

usrp = []

for ip in ip_list:
    rfile = rawfile.RawFile(
        metafile= None,
        abspath=datpath + filename,
        fs=fs, fi=fi, ds=1.0,
        datatype=datatype,
        notes='Data set ' + refname
    )

    usrp += [receiver.Receiver(rfile)]

    rx = usrp[-1]
    rx.load_measurement_logs(dirname = prepath, subdir= 'end-of-%d_usrp'%proc_time+ str(ip))

    dir = prepath + 'eph%d'%ip
    try:
        os.makedirs(dir)
    except:
        pass

    for prn in rx.channels:
        rx.parse_ephemerides(prn_list = [prn], m_start=40)
        if rx.channels[prn].ephemerides is not None:
           rx.channels[prn].ephemerides.save_ephemerides(dir + '/channel%d.mat' % prn,
                                                         dir + '/channel%d.csv' % prn)



