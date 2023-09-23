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
# primal data for a pump
q = np.array([1,7,11])*100/3600 # m3/s
h = np.array([229,210,149]) # h
p = np.array([330,500,556]) # input electric power of pump (kW)
q1 = np.array([1,4.5,8.8])*100/3600
h1 = np.array([138,130,86])
p1 = np.array([160,220,262])
qpipe = np.array([1,200,400,700,1000])/3600 # m3/s
hlosspipe = np.array([1.164,1.54e4,5.569e4,1.597e5,3.15e5])/1000/9.8 # m
fig,ax = plt.subplots(1,1,figsize=(3,2.4))
axright = ax.twinx()
ax.plot(q*3600,h,label='hw,h',color=cl_cloR)
ax.plot(q1*3600,h1,label='lw,h',color=cl_leaveG)
ax.plot(q* 3600,(1547.2653) * q**2+10, color=cl_hairB)
# ax.plot(qpipe*3600,hlosspipe)
axright.plot(q*3600,p,'--',label='hw,P',color=cl_cloR)
axright.plot(q1*3600,p1,'--',label='lw,P',color=cl_leaveG)
# ax.plot(q*3600,hloss,color = cl_cloR,linewidth=3,label='EXP')
# ax.plot(q1*3600,hloss1,color = cl_hairB,linewidth = 1, marker = 'o',label='PWL')
ax.set_xlabel('Volumetric flow rate ($m^3/h$)')
ax.set_ylabel('Head gain (m)')
axright.set_ylabel('Power consumption (kW)')
# ax.legend(loc = 'center left')
# ax.legend()
# axright.legend(loc = 'lower center')
# axright.legend()
fig.tight_layout()
plt.show()