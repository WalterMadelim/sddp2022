import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
from matplotlib import pyplot as plt
# from paper to get a pump characteristic
# 20220923
cl_cloR   = 138/256,33 /256,66/256
cl_lotusR = 207/256,131/256,185/256
cl_leaveG = 96 /256,124/256,56/256
cl_fenceG = 66 /256,100/256,93/256
cl_hairB  = 32 /256,25/256,25/256

color1 = (0.1647, 0.7176, 0.6157)  # 海绿
color2 = (0.2667, 0.7804, 0.9569,)  # 海蓝
color3 = (0.9059, 0.6627, 0.1255,)  # 脏橙
color4 = (0.8784, 0.7882, 0.2431,)  # 暗黄
color5 = (0.7686, 0.3529, 0.1804,)  # 血红
color6 = (0.6902, 0.6784, 0.9765)  # 淡紫色
color7 = (0.3922, 0.4078, 0.8549)  # 紫色
color8 = (0.4784, 0.8000, 0.4392,)  # 草绿色
color9 = (0.0863, 0.7843, 0.8353,)  # 电吉他青
color10 = (0.9961,0,0.5882,) # 桃红
color11 = (0.5412,0.1059,0.6431,) # 紫罗兰
color12 = (0.36862,1,0.003921) # 亮绿
color14 = (0,0.40784,0.70588) # intel蓝
color15 = (0.46274,0.7255,0) # Nvidia绿

Pchp = .1*np.array([2.894+4+2+4,3.348+3.281+2.5+4,2.5+2.5+2.659+2.5,4+2.5+3.139+2.5])
Qchp = .1*np.array([1.91543+0+3.464+0,1.130+0+2.60+0,0.746+0+2.323+2.598,0+2.598+1.490+2.598])
hchp = .1*np.array([0+0+0+0,0.522+0.575+1.2+0,1.2+1.2+1.073+1.2,0+1.2+0.688+1.2])
pc = .5*np.array([.989,.945,0,0],dtype=np.float64)  # negative is charging to CAES; CAES do not provide heat power at the combined small case
pd = .5*np.array([0,0,.99,1],dtype=np.float64)  # negative is charging to CAES; CAES do not provide heat power at the combined small case
Ma = np.array([74.8,98.4,74.0,52.0])
Mh = np.array([75.6,98.4,74.8,50.4])
htou = np.array([100,67.1,100,64.3])
elerd = np.array([0.312,0,0,0.320])
heatrd = np.array([0.240,0.010,0.133,0.651])
Qelerd = 0
Ppump = np.array([.547+.16,.198+.242,.256+.256,.194+.239],dtype=np.float64)

vle, fle, ple = .4, 1.0, 1.6
lE = np.array([vle,vle,fle,ple])
lH = 3/5*lE
vP, fP, pP = 308.9, 626.8, 1044 # yuan/MWh
Pe = np.array([vP,pP,fP,pP],dtype=np.float64)

x = np.arange(1,5)  # the label locations
width = 0.2  # the width of the bars

k = 1.2
fig, ax = plt.subplots(figsize=(3,2.4))
ax.bar(x- k*width, Pchp, width,color = cl_cloR) # ‘$P^{CHP}$’:cl_cloR
ax.bar(x- k*width,pd,width,bottom=Pchp, color=color6) # '$P^{CS,d}$':color6
ax.bar(x- k*width,elerd,width,bottom=Pchp+pd,color=color7) # '$P^{grid}$':color7
ax.bar(x- k*width,-lE,width,color=color2) # '$P^{load}$':color2
ax.bar(x- k*width,-pc,width,bottom=-lE,color=color3) # '$P^{CS,c}$':color3
ax.bar(x- k*width,-Ppump,width,bottom=-lE-pc,color=cl_lotusR) # '$P^{PM}$':cl_lotusR

ax.bar(x , Qchp, width,color=color1) # ‘$Q^{CHP}$’:color1
ax.bar(x , -.25*lE, width,color=color12) # '$Q^{load}$':color12
ax.bar(x , -np.tan(np.arccos(.85))*Ppump, width,bottom=-.25*lE,color=color15) # ‘$Q^{PM}$’:color15

ax.bar(x + k*width, hchp, width,color=color10) # '$H^{CHP}$':color10
ax.bar(x + k*width, heatrd, width,bottom=hchp,color=color11) # '$H^{DHN}$':color11
ax.bar(x + k*width, -lH, width,color=color14) # '$H^{load}$':color14

axright = ax.twinx()
axright.step([.5,1.5,2.5,3.5,4.5],[vP,vP,pP,fP,pP],alpha=.4,color = cl_hairB)
axright.set_ylabel('Price (yuan/MWh)')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Power (MW)')
# ax.set_xticks(x, labels)
fig.tight_layout()
plt.show()

# pump results
# pump,t=       0,p0=( 1033.33| 159.167| 546.667),p1=(     100|     138|     160)
# pump,state=(0|1|0|0|__|0|0|1|0)
# pump,t=       1,p0=( 336.119| 129.605| 198.304),p1=( 650.547| 111.031| 241.836)
# pump,state=(0|0|1|0|__|0|0|0|1)
# pump,t=       2,p0=( 806.046| 94.0677| 255.501),p1=( 807.287| 93.9323|  255.61)
# pump,state=(0|0|0|1|__|0|0|0|1)
# pump,t=       3,p0=( 309.881| 130.538| 194.047),p1=( 623.452| 113.987| 239.455)
# pump,state=(0|0|1|0|__|0|0|0|1)

#
# fig,ax = plt.subplots(1,1,figsize=(3,2.4))
# axright = ax.twinx()
# ax.plot(q*3600,h,label='hw,h',color=cl_cloR)
# ax.plot(q1*3600,h1,label='lw,h',color=cl_leaveG)
# axright.plot(q*3600,p,'--',label='hw,P',color=cl_cloR)
# axright.plot(q1*3600,p1,'--',label='lw,P',color=cl_leaveG)
# # ax.plot(q*3600,hloss,color = cl_cloR,linewidth=3,label='EXP')
# # ax.plot(q1*3600,hloss1,color = cl_hairB,linewidth = 1, marker = 'o',label='PWL')
# ax.set_xlabel('Volumetric flow rate ($m^3/h$)')
# ax.set_ylabel('Head gain (m)')
# axright.set_ylabel('Power consumption (kW)')
# # ax.legend(loc = 'center left')
# # ax.legend()
# # axright.legend(loc = 'lower center')
# # axright.legend()
# fig.tight_layout()
# plt.show()