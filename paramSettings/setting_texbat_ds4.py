# -*- coding: utf-8 -*-
# texbat_ds4
import numpy as np

refname = 'texbat_ds4'  # 数据名
datpath = 'D:\\academic\\DPEdata\\'
filename = 'texbat_ds4.bin'
descriptor = 'test1026GRID_PRN'
fs = 25e6
fi = 0.0e6
datatype = np.dtype([('i', np.int16), ('q', np.int16)])

start_time = 160
proc_time = 40
DPE_start_time = 190
DPE_lead_time = DPE_start_time - start_time  # must make 3<=DPE_lead_time<=proc_time
DPE_run_time = 1
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

acq_only = False
load_acq = True
load_trk = True
# prn_list = [8, 10, 11, 12, 14, 15, 20, 21, 24, 25, 27, 31, 32] # Simulated
prn_list = range(1, 33)  # Simulated
# prn_list = [15] # Simulated

execfile('paramSettings/default_settings.py')