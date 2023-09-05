# -*- coding: utf-8 -*-
execfile('setting_finland.py')
# execfile('setting_210Atest1.py')
### Main code starts

from pythonreceiver.libgnss.constants import *
from pythonreceiver.libgnss import rawfile,utils,satpos
from pythonreceiver.scalar import channel, correlator, discriminator, loopfilter,naveng
from pythonreceiver import receiver
#import printer
import threading

import numpy as np
import scipy.io as sio
import time,os
import csv

scalar_usrp = [
    receiver.Receiver(\
        rawfile.RawFile(
            metafile = None,
            abspath  = datpath + filename,
            fs = fs, fi = fi, ds = 1.0,
            datatype = datatype,
            notes = 'Data set '+ refname
        )
        , mcount_max = run_time * 1000 + 1
    ) for ip in ip_list
]

print('Start scalar tracking @',init_time,'for',run_time)

class scalar_thread (threading.Thread):
    def __init__(self,rx,ip,lock):
        threading.Thread.__init__(self)
        self.rx = rx
        self.ip = ip
        self.running = True
        self.counter = 0
        self.total = 0

    def run(self):
        print('Thread Launched')
        first_dir  =  'end-of-1_usrp'+ str(self.ip)
        second_dir =  'end-of-%d_usrp'%proc_time+ str(self.ip)

        self.rx.add_channels(prn_list)
        self.rx.rawfile.seek_rawfile(init_time * 1000 * self.rx.rawfile.S)
        self.rx.scalar_acquisition(prn_list)

        if acq_only:
            self.running = False
            return

        try:
            self.rx.scalar_track(mtrack=run_time * 1000)
        finally:
            lock.acquire()
            self.rx.save_measurement_logs(dirname = prepath,subdir= second_dir)
            lock.release()
        self.rx.save_measurement_logs(dirname = prepath,subdir= second_dir)

        self.rx.store_ref_mcount()

        self.running = False

        del_clist = []
        for prn in self.rx.channels:
            try:
                self.rx.parse_ephemerides(prn_list = [prn],m_start = 40)
            except:
                del_clist += [prn]
        self.rx.del_channels(del_clist)

        self.rx.save_scalar_handoff(prn_list, dirname=prepath, subdir=first_dir)

        # 修改
        # os.makedirs(prepath + 'eph%d'%self.ip)
        dir = prepath + 'eph%d' % self.ip
        if not os.path.exists(dir):
            os.makedirs(dir)
        for prn in self.rx.channels:
            if self.rx.channels[prn].ephemerides is not None:
                self.rx.channels[prn].ephemerides.save_ephemerides(prepath + 'eph%d/channel%d.mat'%(self.ip,prn), prepath + 'eph%d/channel%d.csv'%(self.ip,prn))

        print ('Scalar Tracking concluded.')
        return

lock= threading.Lock()
s_threads = [scalar_thread(rx,ip_list[i],lock) for i,rx in enumerate(scalar_usrp)]
for st in s_threads:
    st.start()

### Block for scalar
while any([t.running for t in s_threads]):
    print ('Scalar running; total time',run_time)
    print ('Current time',[rx._mcount/1000.0 for rx in scalar_usrp])
    time.sleep(30)

print ('Scalar tracking completed.  Ephemeris saved.')

