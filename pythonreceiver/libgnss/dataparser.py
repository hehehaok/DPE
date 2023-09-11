# -*- coding: utf-8 -*-


from . import ephemeris as eph
import numpy as np

preamble = np.array([-1,1,1,1,-1,1,-1,-1])
preamble_cp = np.kron(preamble, np.ones(20))
        
def parse_ephemerides(channel, m_start, m_end):
    
    cp_start = int(channel.cp[m_start])  # 星历解算的起始码周期
    cp_end = int(channel.cp[m_end])  # 星历解算的终止码周期
    
    if len(np.where(np.diff(channel.cp[m_start:m_end])!=1)[0]) != 0:  # 确保中间码周期cp是连续的，中间没有出现跳过的情况
        print('WARNING: possible code period slip occured.')

    cp = np.arange(cp_start,cp_end)
    iP = channel.cp_sign[cp]  # 获得各码周期的导航电文符号

    preamble_correlations = np.correlate(iP,preamble_cp,'valid')  # 获得各码周期的导航电文符号与导航电文的同步码做相关来找帧头
    preamble_locations = np.where(np.abs(preamble_correlations)>153)[0][:]  # 相关结果大于阈值153即表示可能为帧头，记录码周期所在位置
    preamble_locations_set = set(preamble_locations)  # 注意一下preamble_locations变为集合类型之后会存在乱序，即不一定再严格按照升序排列
    subframe_found = False

    for test in preamble_locations:
        if test < 40:
            test = test + 6000  # 由于当前子帧的奇偶校验需要上一个子帧的最后两个导航电文比特，因此需要确保起始帧头位置大于2个导航电文比特的位置，即40个伪码周期
        test_set = set([test,test+6000,test+6000*2,test+6000*3,test+6000*4])
        # 一个子帧6秒，即6000个码周期，一个完整的帧包含5个子帧，如果能够找到连续的5个子帧即说明已具备解算所需的所有星历参数
        # 因此假设找到5个连续的子帧帧头
        print(test_set,test_set.issubset(preamble_locations_set))
        # 如果假设的连续5个子帧帧头所在集合为已找到的所有帧头的子集，则说明假设成立，找到了连续的5个子帧
        if test_set.issubset(preamble_locations_set):
            subframe_found = True
            subframe_locations = sorted(list(test_set))  # 对乱序的帧头位置集合重新按升序进行排序
            break

    if subframe_found == False:
        print('Ephemerides decoding unsuccessful: Unable to locate 5 preambles within given range.')
        return

    subframe_polarities = np.sign(preamble_correlations[subframe_locations])  # 得到各码周期的导航电文符号与导航电文的同步码做相关的结果的符号
    print('polarities: '+str(subframe_polarities))

    if not (np.all(subframe_polarities == [-1,-1,-1,-1,-1]) or np.all(subframe_polarities == [1,1,1,1,1])):
        print('Ephemerides decoding warning: Bit flip occured in between subframes')  # 相关结果的符号应该是一样的

    navbits_all = np.reshape(iP[subframe_locations[0]:(subframe_locations[0]+6000*5)],(1500,20))
    navbits_all = np.sign(np.sum(navbits_all,1))
    # 每20个伪码周期构成一个导航电文，这里从帧头开始每20个伪码周期进行一次累加从而得到真正的导航电文比特

    subframes = []
    first_d29d30 = np.reshape(iP[(subframe_locations[0]-40):(subframe_locations[0])],(2,20))
    first_d29d30 = np.sign(np.sum(first_d29d30,1))
    d29 = first_d29d30[0]
    d30 = first_d29d30[1]
    # 将第一个帧头之前的40个伪码周期进行累加，得到第一个帧头的上一个子帧的最后两个导航电文比特，用于第一个子帧的奇偶校验
    # 具体奇偶校验可以参考谢钢 GPS原理与接收机设计P34

    reload(eph)
    for subframeNr in range(5):
        words = []
        wordNr = 0    
        polarity = subframe_polarities[subframeNr]
        bits = navbits_all[(subframeNr*300+(wordNr*30)):(subframeNr*300+(wordNr*30+30))]     
        words.append(eph.Word(polarity,d29,d30,bits))
        for wordNr in range(1,10):
            bits = navbits_all[(subframeNr*300+(wordNr*30)):(subframeNr*300+(wordNr*30+30))]
            words.append(eph.Word(words[-1].d30,words[-1].d29,words[-1].d30,bits))
        d29 = words[-1].d29
        d30 = words[-1].d30
        subframes.append(eph.Subframe(subframe_locations[subframeNr]+cp[0],words))

    parsed_ephemerides  = eph.Ephemerides(subframes)
    channel.ephemerides = parsed_ephemerides
    
    return
