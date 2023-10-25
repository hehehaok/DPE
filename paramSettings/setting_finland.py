# -*- coding: utf-8 -*-
# 芬兰数据
import numpy as np

refname    = 'finland_cleanStatic'  # Simulated
datpath  = 'D:/academic/DPEdata/'
filename = 'cleanStatic_gps_finland.dat'
descriptor = 'test1025GRID'
fs = 26e6
fi = 6.39e6
datatype = np.dtype([('i', np.int8)])

start_time    = 0
proc_time     = 40
DPE_start_time = 30
DPE_lead_time = DPE_start_time - start_time  # must make 3<=DPE_lead_time<=proc_time
DPE_run_time = 0.2
DPE_interval = 0.02  # 进行DPE的间隔
DPE_corr_save_interval = 0.02  # 保存DPE流形结果图的间隔

dpe_plan = 'GRID'  # 'GRID'-网格法  'ARS' - accelerated random search
# ************* GRID方法参数 *****************
grid_param = {'N': 25**4}
# ************* GRID方法参数 *****************
# ************* ARS方法参数 *****************
ars_param = {'dmax': 10, 'dmin': 1,
             'dmaxt': 5, 'dmint': 0.5,
             'cf': 2, 'N_Iter': 1000}
# ************* ARS方法参数 *****************

acq_only      = False
load_acq = True
load_trk = True
prn_list = [10, 11, 13, 15, 17, 19, 20, 24, 28, 30] # Simulated
# prn_list = [13, 15, 17, 24, 28] # Simulated
# prn_list = [15] # Simulated

execfile('paramSettings\\default_settings.py')