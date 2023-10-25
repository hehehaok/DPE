# -*- coding: utf-8 -*-


from .constants import *

import numpy as np
from math import ceil
import os.path

class RawFile():
    """A container for a binary-formatted (.bin or .dat) GNSS raw signal file.

    RawFile records the data collection settings of the GNSS raw signal file as attributes.
    RawFile returns snippets of the GNSS raw signal as numpy ndarrays.
    """

    def __init__(self, metafile = None, rawfile = None,
                 abspath = None, fs = None, fi = None, ds = None, datatype = None, notes = None, verbose = True):
        """Constructs a RawFile object with a binary-formatted (.bin or .dat) GNSS raw signal file.

        Sets data collection settings as attributes, open file stream, sets rawsnippet settings as attributes.
        """

        # Reads metafile
        if metafile is not None:
            abspath, fs, fi, ds, datatype, notes = self.get_settings_from_metafile(
                metafile_name = metafile, rawfile_name = rawfile, verbose = verbose)

        # Sets rawfile settings
        self.set_settings(abspath, fs, fi, ds, datatype, notes)  # 读取数据相关的设置

        # Open rawfile
        self.open_rawfile()  # 打开指定路径下的数据文件

        # Sets rawsnippet settings
        self.set_rawsnippet_settings(T_CA,T_CA)

        return

    def get_settings_from_metafile(self, metafile_name, rawfile_name, verbose = True):

        metafile = open(metafile_name, 'r')

        while not metafile.closed:

            newline = metafile.readline()

            #if newline[11:-1] == rawfile_name:
            meta_fname = newline[newline.find('=')+1:].strip()
            if meta_fname == rawfile_name:

                newlines = metafile.readlines(5)

                abspath = os.path.join(os.path.dirname(metafile.name), meta_fname)
                assert os.path.exists(abspath), 'Generated path \'%s\' does not exist.'%(abspath)
                fs = float(newlines[0][5:-1])
                fi = float(newlines[1][5:-1])
                ds = float(newlines[2][5:-1])
                datatype = eval(newlines[3][11:-1])
                notes = newlines[4][8:-1]
                #fs = float(newlines[0].split(' ')[2][:-1])
                #fi = float(newlines[1].split(' ')[2][:-1])
                #ds = float(newlines[2].split(' ')[2][:-1])

                metafile.close()

                if verbose:
                    print('Rawfile settings:')
                    print('  Absolute path = %s'%(abspath))
                    print('  Sampling frequency = %.4fe6'%(fs/1.0e6))
                    print('  Intermediate frequency = %.4fe6'%(fi/1.0e6))
                    print('  Doppler sign = %.1f'%(ds))
                    print('  Data type = %s'%(str(datatype)))
                    print('  Notes = %s'%(notes))
                    print('Closing metafile.\n')

                return abspath, fs, fi, ds, datatype, notes

            if newline == '':

                metafile.close()

                print('File descriptor not found.')
                print('  Please check text formatting in metafile.')
                print('Closing metafile.\n')

        print 'Default values applied:'
        return None, None, None, None, None, None

    def set_settings(self, abspath, fs, fi, ds, datatype, notes):

        if fi != 0:
            print('Warning, not zero-IF complex sampling.')  # 如果零中频则进行警告(原代码是默认信号是零中频的)

        self.abspath = abspath  # 数据文件的绝对路径
        self.fs = fs  # 采样频率
        self.fi = fi  # 中频频率
        self.ds = ds
        self.fcaid = ds*F_CA/F_L1  # 码频率与载波频率的比值，后面用来通过载波多普勒频移fi修正码频率fc
        self.datatype = datatype  # 采样点数据的数据类型
        self.notes = notes

        if self.datatype.fields is not None:  # 根据不同的数据类型，使用不同的读取数据的方法
            if (self.datatype.fields.keys() == ['arg_pi4']):
                self.format_rawsnippet = self.format_rawsnippet_datatype_arg_pi4
            if (self.datatype.fields.keys() == ['i','q']):
                self.format_rawsnippet = self.format_rawsnippet_datatype_complex
            # 修改
            if (self.datatype.fields.keys() == ['i']):
                self.format_rawsnippet = self.format_rawsnippet_datatype_real
            else:
                print('Unknown datatype: %s'%(str(self.datatype)))

        return

    def open_rawfile(self):

        if (not hasattr(self, 'rawfile')) or (hasattr(self, 'rawfile') and self.rawfile.closed):
            self.rawfile = open(self.abspath, 'rb')
            print('Opened rawfile = \'%s\'.\n'%self.abspath)
        else:
            print('Rawfile = \'%s\' is already opened.\n'%self.rawfile.name)
        return

    def close_rawfile(self):

        if (hasattr(self, 'rawfile')) and (not self.rawfile.closed):
            self.rawfile.close()
            print('Closed rawfile = \'%s\'.\n'%self.abspath)
        else:
            print('Rawfile = \'%s\' is already closed.\n'%self.rawfile.name)
        return

    def seek_rawfile(self, S_skip, whence = 1):
        self.rawfile.seek(int(S_skip*self.datatype.itemsize),whence)
        self.rawfile_samp = self.rawfile.tell()/self.datatype.itemsize
        self.rawfile_time = self.rawfile_samp/self.fs
        #print('Rawfile time = %.3fs'%(self.rawfile_time))
        return

    def return_num_bytes_read(self):
        return self.rawfile.tell()

    def update_rawsnippet(self):

        self.rawfile_samp = self.rawfile.tell()/self.datatype.itemsize
        self.rawfile_time = self.rawfile_samp/self.fs

        rawsnippet = np.fromfile(self.rawfile, self.datatype, self.S)
        self.rawsnippet = self.format_rawsnippet(rawsnippet)
        return

    def format_rawsnippet(self):
        pass

    def format_rawsnippet_datatype_arg_pi4(self, rawsnippet):
        # assumes fi = 0, else shift to baseband
        return np.exp(1j*(rawsnippet['arg_pi4']*(np.pi/4.0)))

    def format_rawsnippet_datatype_complex(self, rawsnippet):
        # assumes fi = 0, else shift to baseband
        return rawsnippet['i']+1j*rawsnippet['q']

    def format_rawsnippet_datatype_real(self, rawsnippet):
        # assumes fi = 0, else shift to baseband
        return rawsnippet['i']

    def set_rawsnippet_settings(self, T, T_big, verbose=True):

        self.T = T                                 # float (s) 每个历元读取数据的时间
        self.N = int(ceil(self.T/T_CA))           # number of 1 millisecond blocks T占多少个伪码周期
        self.S = int(ceil(T*self.fs))             # int   (samples) T时间内的采样点个数
        if self.N == 10:
            self.S = int(ceil(self.S/10.0)*10)        # 由于210A的采样率为16367667，为了让S最后为10的倍数，不得已这样做
        self.samp_idc = np.arange(0,self.S)        # int   array of indices (samples)
        # 构造一个长度为S的数组(0 - S-1)，表示采样点的序列
        self.time_idc = self.samp_idc / self.fs    # float array of indices (s)
        # 构造一个长度为S的数组(0 - S-1)，表示采样点时间的序列，T时间内有S个采样点
        # 这个数组就是存储这S个采样点的采样时间，该序列中的每个值表示相应采样点对应的采样时间
        self.code_idc = self.time_idc * F_CA       # float array of indices (chips)
        # 构造一个长度为S的数组(0 - S-1)，表示采样点伪码的序列，T时间内有S个采样点，这个数组就是存储这S个采样点在一个伪码序列的第几个码片
        # 注意这里没取模，所以有可能大于伪码序列长度1023，该序列中的每个值表示相应采样点对应的码片值

        code_idc = np.arange(0,int(ceil(T_CA*self.fs))) / self.fs * F_CA # (chips)
        # 与self.code_idc类似，但是是T_CA(0.001秒)，即一个伪码周期时间内采样点的伪码码片序列
        code_idc = np.fft.fftshift(np.where(code_idc>=L_CA/2.0,code_idc-L_CA,code_idc))
        self.code_fftidc = code_idc
        # 例：[0~499] -> [0~249 -250~-1] -> [-250~249] 伪码FFT的频率轴，让频域的零频率移到FFT序列的中心

        self.carr_fftpts = 8 * (1 << (self.S).bit_length())
        # (self.S).bit_length()表示采样点个数的位数，即采样点数量占二进制多少位
        # 1 << (self.S).bit_length()表示将数字1左移(self.S).bit_length()位，即求得对采样点个数进行FFT所需要的点数(2的幂次方)
        # 最后乘以8猜测是对应博士论文的第40页3.3.4部分，是经验上的取值
        self.carr_fftidc = np.fft.fftshift(np.fft.fftfreq(n = self.carr_fftpts, d = 1.0/self.fs))
        # 载波FFT的频率轴，且让频域的零频率移到FFT序列的中心，该序列中的每个值表示相应采样点对应的频率值

        assert T_big >= T, 'Duty-cycle interval T_big=%.2fs has to be larger than or equal to T=%.2fs.'%(T,T_big)

        self.T_big  = T_big
        self.T_skip = self.T_big - self.T

        self.S_big  = int(ceil(self.T_big*self.fs))
        if self.N == 10:
            self.S_big = int(ceil(self.S_big / 10.0) * 10)
        self.S_skip = self.S_big - self.S

        if verbose:
            print('Rawsnippet settings:')
            print('T_big = %.3fs, T = %.3fs, T_skip = %.3fs.'%(self.T_big, self.T, self.T_skip))
            print('S_big = %d, S = %d, S_skip = %dsamples.\n'%(self.S_big, self.S, self.S_skip))

        return

    # Need to write save and load functions
