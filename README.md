# DPE
## 环境信息
- python       2.7.13
- numpy        1.15.1
- matplotlib   2.2.5
- scipy        1.1.0
## 使用数据
- 数据来源：FGI-GSRx提供的参考数据 [数据链接](https://tiedostopalvelu.maanmittauslaitos.fi/tp/julkinen/lataus/tuotteet/FGI-GSRx-OS-DATAFILES/FGI-GSRx%20Raw%20IQ%20Data/GPSL1GalileoE1) 【注意下载后将文件名改为：static_opensky_20230421_190000_usrp6_26000kHz.dat】
- 对应配置参数文件：settings2.py 【注意要将其中的datpath改为数据所在路径】
- 中频频率：6.39MHz
- 采样频率：26MHz
- IQ采样：否
- 单个采样点数据大小：8bit
- 捕获卫星列表：[10, 11, 13, 15, 17, 19, 20, 24, 28, 30]
- 真实位置：[trueLat, trueLong, trueHeight] = [60.161086788889, 24.545448080556, 54.1640000026673]  

## 运行步骤

1. 配置好需要的环境信息
2. 在`paramSettings`中创建配置文件，配置好待处理中频数据的数据信息
3. 修改`DPE_SDR.py`顶部需要执行的配置文件
4. 运行`DPE_SDR.py`即可

## 目录结构

- dataProcess DPE结果画图
- paramSettings 用于存放配置文件
- pre-simulator 标量跟踪结果(运行完一次DPE后生成)
- post-simulator DPE跟踪结果(运行完一次DPE后生成)
- pythonreceiver 接收机源码
