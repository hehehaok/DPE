#YYYYMMDD_HHMMSS_usrp%d_xxxxkHz.dat

refname    = 'static_opensky_'  # Simulated
datname    = '20230421_190000'
descriptor = 'Talbot_rooftop_4hr_T_big=1'
fs = 26e6
fi = 6.39e6
# fi = 0.0e6
cudarecv_handoff = 'handoff_params_usrp6' # Simulated

ip_list       = [6] # Simulated
start_time    = 0
#  this will run from start_time to (start_time+proc_time)
proc_time     = 40
max_lead_time = 0
weekno        = 2258 # If running simulated data
acq_only      = False
# prn_list = [10, 11, 13, 15, 17, 19, 20, 24, 28, 30] # Simulated
prn_list = [13, 15, 17, 24, 28] # Simulated
# prn_list = [15] # Simulated

datpath  = 'D:/study/DPEdata/'
predir   = './pre-simulator/'
postdir  = './post-simulator/'



### No need to change: house keeping ###

descriptor = descriptor + '_%dUSRP' %len(ip_list)

import os
dir_req  = [datpath,predir,postdir]#,'netcsv','html']
for d in dir_req:
    if not os.path.exists(d):
        os.makedirs(d)

### More housekeeping...

if start_time < 4:
    init_time = start_time
elif start_time < (max_lead_time + 4):
    init_time = 4
else:
    init_time = start_time - max_lead_time

run_time = proc_time  + (start_time - init_time)
prefix   = datname    + '_skip%d_start%d/'%(init_time,start_time)+'6605_7200/'
prepath  = predir     + prefix
postfix  = descriptor + 'proc%ds_6605_7200/'%run_time
postpath = postdir    + postfix
try:
    os.makedirs(postpath)
except:
    if os.listdir(postpath):
        print 'Warning:',postpath,'not an empty directory.'