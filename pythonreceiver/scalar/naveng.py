# -*- coding: utf-8 -*-



from ..libgnss import utils, satpos
from ..libgnss.constants import *

import numpy as np

def calculate_nav_soln(receiver, prn_list=None, mc=None, rxTime0=None, rxPos0=None, pOut=False):

    if prn_list is None:
        prn_list = sorted(receiver.channels.keys())

    if mc is None:
        mc = receiver._mcount

    numChannels = len(prn_list)

    codeIntDiff  = np.zeros(numChannels)
    codeFracDiff = np.zeros(numChannels)
    transmitTime = np.zeros(numChannels)

    sats_ECEF    = np.matrix(np.zeros((8,numChannels)))

    doppler      = np.zeros(numChannels)

    for prn_idx, prn in enumerate(prn_list):

        codeIntDiff[prn_idx]  = (receiver.channels[prn].cp[mc] \
                                 - receiver.channels[prn].ephemerides.timestamp['cp'])*T_CA
        # 伪码周期个数对应的时间=(当前所在伪码周期序号-子帧起始位置处的伪码周期序号)*伪码周期
        codeFracDiff[prn_idx] = (receiver.channels[prn].rc[mc])/F_CA
        # 伪码码相位对应的时间=码相位(chip)/码频率=码相位(chip)/伪码长度(1023)*伪码周期(0.001s)
        transmitTime[prn_idx] = receiver.channels[prn].ephemerides.timestamp['TOW'] \
                                 + codeIntDiff[prn_idx] + codeFracDiff[prn_idx]
        # 卫星信号的发送时间=子帧起始位置处的时间(周内秒)+从子帧起始位置开始经过的伪码周期个数对应时间+当前伪码周期内的码相位对应时间
        # 具体可以参考谢钢 GPS原理与接收机设计 P72

        clkb, clkd = satpos.satellite_clock_correction(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx])
        # 根据星历参数和信号发送时间计算钟差钟漂，注意这里的钟差钟漂是针对卫星的，即用于修正卫星信号发送时间的，而不是接收机本地时钟的钟差钟漂

        sats_ECEF[:,prn_idx] = satpos.locate_satellite(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx]-clkb, clkb, clkd)
        # 计算卫星状态向量[位置状态XYZ 钟差 速度状态VxVyVz 钟漂]（在地心地固坐标系(ECEF)下）

        doppler[prn_idx]      = (receiver.channels[prn].fi[mc])*receiver.rawfile.ds
        # 多普勒频移

    # Channels with a larger codeIntDiff + codeFracDiff are from satellites
    # which are closer in distance. Thus at a fixed received time instance,
    # the signals in those Channels correspond to a later transmitted time.

    # We need to remove satellite clock errors from the pseudorange
    # To reduce the number of for-loops, we determine that together with the
    # satellite positions in the Earth-Centered-Earth-Fixed frame (ECEF)

    if rxTime0 is None:
        # It takes about 68ms for the signal to travel from the closest
        # satellite to the receiver.
        rxTime = max(transmitTime)+0.068
        # 接收机接收信号的时间=发射时间+68ms(估算)
        # GPS GEO卫星到地球的时间为68ms左右
        #print('Initializing rxTime as: %.3f'%(rxTime))
    else:
        rxTime = rxTime0

    # organize observations into pseudoranges and pseudorates
    pseudoranges = (C)*(rxTime - transmitTime) + (C)*(np.asarray(sats_ECEF)[3,:])
    # 计算伪距=光速*(接收时间-发射时间)+光速*钟差
    pseudorates  = (-C/F_L1)*(doppler) + (C)*(np.asarray(sats_ECEF)[7,:])
    # 计算伪距速率=-光速*(多普勒频移/载波频率)+光速*钟漂

    transmitTime = transmitTime - np.asarray(sats_ECEF)[3,:]
    # 用钟差修正信号发送时间

    # rotate satellite positions to Earth-Centered-Inertial (ECI)
    t_c = rxTime
    # 信号接收时间

    sats_ECI = np.matrix(np.zeros((8,numChannels)))

    for prn_idx, prn in enumerate(prn_list):
        sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx],t_gps=transmitTime[prn_idx], t_c=t_c)
    # 地心地固坐标系->地心惯性坐标系(在转换过程中要考虑地球的自转角=地球自转角速度*信号传输时间)
    # 这里感觉不像是完整的地心地固->地心惯性的坐标系转换，只考虑了地球自转角，应该只是为了后续在最小二乘的时候不用再考虑地球自转，参考鲁郁书P308

    # perform navigation solution estimation through iterative least squares in ECI,
    # get actual GPS received time, rotate back to ECEF
    posvel_ECI  = perform_least_sqrs(sats_ECI, pseudoranges, pseudorates=pseudorates, rxPos0=rxPos0)
    # 使用最小二乘计算接收机的8维状态向量，rxPos0为参考位置
    # 注意这里得到的posvel_ECI的第4个元素不再是以秒为单位的钟差，而是钟差乘以光速，且这里的钟差是接收机本地时钟的钟差

    rxTime_a    = rxTime - posvel_ECI[3,0]/C
    # 修正接收机本地时钟的钟差(一开始我们对接收机的本地时钟采用的是发射时间+68ms的估计值，在这里得到修正)
    posvel_ECEF = utils.ECI_to_ECEF(posvel_ECI, t_gps=rxTime_a, t_c=t_c)
    # 将接收机状态向量从地心惯性系转换回地心地固系(此时地心地固系的接收机坐标对应的时间为修正接收机时钟后的时间)

    # rotate all states to receiver's ECI
    t_c = rxTime_a
    posvel_ECI = utils.ECEF_to_ECI(posvel_ECEF, t_gps=rxTime_a, t_c=t_c)
    for prn_idx, prn in enumerate(prn_list):
        sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx], t_gps=transmitTime[prn_idx], t_c=t_c)
    # 将卫星和接收机的状态向量全部转换到地心惯性系

    if pOut == True:
        return rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI, pseudoranges, pseudorates

    return rxTime_a, rxTime, posvel_ECEF, posvel_ECI, sats_ECI

def get_satellite_positions(receiver, prn_list=None, mc=None, t_c=None):

    if prn_list is None:
        prn_list = sorted(receiver.channels.keys())

    if mc is None:
        mc = receiver._mcount

    numChannels = len(prn_list)

    codeIntDiff  = np.zeros(numChannels)
    codeFracDiff = np.zeros(numChannels)
    transmitTime = np.zeros(numChannels)

    sats_ECEF    = np.matrix(np.zeros((8,numChannels)))

    for prn_idx, prn in enumerate(prn_list):

        codeIntDiff[prn_idx]  = (receiver.channels[prn].cp[mc] \
                                 - receiver.channels[prn].ephemerides.timestamp['cp'])*T_CA
        codeFracDiff[prn_idx] = (receiver.channels[prn].rc[mc])/F_CA
        transmitTime[prn_idx] = receiver.channels[prn].ephemerides.timestamp['TOW'] \
                                 + codeIntDiff[prn_idx] + codeFracDiff[prn_idx]

        clkb, clkd = satpos.satellite_clock_correction(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx])
        sats_ECEF[:,prn_idx] = satpos.locate_satellite(receiver.channels[prn].ephemerides,
                                                               transmitTime[prn_idx]-clkb, clkb, clkd)
        transmitTime[prn_idx] = transmitTime[prn_idx] - clkb

    if t_c is not None:

        sats_ECI = np.matrix(np.zeros((8,numChannels)))

        for prn_idx, prn in enumerate(prn_list):
            sats_ECI[:,prn_idx] = utils.ECEF_to_ECI(sats_ECEF[:,prn_idx],
                                                    t_gps=transmitTime[prn_idx], t_c=t_c)

        return sats_ECI, transmitTime

    return sats_ECEF, transmitTime

def perform_least_sqrs(sats, pseudoranges, pseudorates=None, iterations=10, rxPos0=None):

    assert np.shape(sats)[0] == 8,\
         "Error: Expected satPos to be formatted as (8,N) np.matrix,\
         given satPos is formatted as:"+str(np.shape(satPos))
    satPos  = sats[0:3,:]
    satVel  = sats[4:7,:]
    numSats = np.shape(satPos)[1]

    assert np.shape(pseudoranges) == (numSats,), \
        "Error: Expected pseudoranges to be of shape ("+str(numSats)+",),\
         given pseudoranges has shape:"+str(np.shape(pseudoranges))

    if pseudorates is not None:
        assert np.shape(pseudorates) == (numSats,), \
            "Error: Expected pseudorates to be of shape ("+str(numSats)+",),\
             given pseudorates has shape:"+str(np.shape(pseudorates))

    # Initialize the helper variables.
    if rxPos0 is None:
        rxPos = np.matrix(np.zeros((4,1)))
    else:
        assert np.shape(rxPos0) == (4,1), \
             "Error: Expected rxPos0 to be of shape (4,1), \
             given rxPos0 has shape:"+str(np.shape(rxPos0))
        rxPos = rxPos0

    A = np.matrix(np.zeros((numSats,4)))
    A[:,3] = np.matrix(np.ones((numSats,1)))
    b = np.matrix(np.zeros((numSats,1)))

    # Iteratively get the receiver location.
    for i in range(iterations):

        # Gather information from each satellite to initialize the b vector and A matrix.
        for idx in range(numSats):

            # Fill in values to b vector.
            b[idx,0] = pseudoranges[idx] - (np.linalg.norm(satPos[:,idx]-rxPos[0:3,0]) + rxPos[3,0])

            # Fill in values to A matrix (linearized equation).
            A[idx,0:3] = (-(satPos[:,idx] - rxPos[0:3,0])/(np.linalg.norm(satPos[:,idx]-rxPos[0:3]))).T

        # If A is not full rank, we cannot proceed.
        if np.linalg.matrix_rank(A) != 4:
            print('Error: the linear least squares estimate of the receiver',
                  'position requires that the geometry A matrix be of rank 4. ',
                  'A is currently of rank: '+str(np.linalg.matrix_rank(A)))
            print('A matrix:'+str(A))

        # Run linear least squares to get the position update.
        x, residuals, rank, s = np.linalg.lstsq(A, b, rcond= None)

        # Now apply the position update.
        rxPos = rxPos + x
        #print('position_calculation',i,np.linalg.norm(x), residuals, rank)

        if np.linalg.norm(x) < 1.0e-7:
            break

    if np.linalg.norm(x) > 1.0e-5:
        print('Error: update is greater than 1.0e-5m after '+str(iterations)+' iterations.')

    # Initialize the helper variable.
    rxVel = np.matrix(np.zeros((4,1)))

    for idx in range(numSats):
        # Compute the updated los.
        los = ((satPos[:,idx] - rxPos[0:3,0])/(np.linalg.norm(satPos[:,idx]-rxPos[0:3]))).T

        # Fill in values to A matrix (unit Line-Of-Sight vectors)
        A[idx,0:3] = -los

        # Fill in values to the b vector
        b[idx] = pseudorates[idx] - los.dot(satVel[:,idx])

    # If A is not full rank, we cannot proceed and simply return zeros.
    if np.linalg.matrix_rank(A) != 4:
        print('Error: the linear least squares estimate of the receiver',
              'velocity requires that the geometry A matrix be of rank 4. ',
              'A is currently of rank: '+str(np.linalg.matrix_rank(A)))
        print('A matrix:'+str(A))

    # Run linear least squares to get the velocity.
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond = None)
    #print('velocity calculation',np.linalg.norm(x), residuals, rank)

    rxVel = x

    retMat = np.matrix(np.zeros((8,1)))
    retMat[0:4,0] = rxPos
    retMat[4:8,0] = rxVel
    return retMat
