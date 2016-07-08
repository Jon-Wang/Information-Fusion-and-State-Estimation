#Information Fusion and State Estimation
##信息融合与状态估计

主要是针对多传感器多时滞（包括状态之后和观测滞后）系统，基于**Kalman滤波**和**现代时间序列分析方法**，利用**集中式融合估计**、**分布式融合估计**（按矩阵加权、按对角阵加权、按标量加权）、
**协方差交叉融合**等方法实现对状态的融合估计。

程序包括：
1. [Covariance Intersection Fusion](https://github.com/Jon-Wang/Information-Fusion-and-State-Estimation/tree/master/Covariance%20Intersection%20Fusion)
Sequential Covariance Intersection Fusion Kalman Filter for Multi-Sensor Systems with Multiple Time Delayed Measurements
利用斜方差交叉(CI)融合方法进行状态估计。
2. [Improved Covariance Intersection](https://github.com/Jon-Wang/Information-Fusion-and-State-Estimation/tree/master/Improved%20Covariance%20Intersection)
Improved covariance intersection fusion Kalman filter for multi-sensor systems with multiple time delayed measurements(带观测滞后多传感器系统的改进协方差交叉融合Kalman滤波器)
利用改进后的协方差交叉融合(ICI)方法实现对状态的估计，相较于原来的CI融合算法，可以提高精度。
3. [Modern Time Series Analysis Method](https://github.com/Jon-Wang/Information-Fusion-and-State-Estimation/tree/master/Modern%20Time%20Series%20Analysis%20Method)
Modern Time Series Analysis Method for Multi-Sensor Systems with Time Delayed Measurements
基于现代时间序列分析方法，对局部传感器构造ARMA信息模型，利用射影定理和白噪声估值器，得到局部状态估计，然后进行融合。
4. [Time Delayed with Correlated Noise](https://github.com/Jon-Wang/Information-Fusion-and-State-Estimation/tree/master/Time%20Delayed%20with%20Correlated%20Noise)
SCI Fusion Estimations for Multi-Sensor Time-Delay Systems with Correlated Noise(带相关噪声多传感器多时滞系统的SCI融合估值器)
为了避免噪声相关带来的推导上的复杂性，先将带相关噪声的系统转化为带不相关白噪声的系统，然后再进行融合。
