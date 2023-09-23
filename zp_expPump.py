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
q = np.array([1,2,3,4,5,6,7,8,8.5,9,9.216,10,11])*100/3600 # m3/s
h = np.array([229,228,227.3,225,223,216,210,193,186.5,180,177,165,149]) # h
p = np.array([330,360,390,420,450,475,500,515,525,530,532,548,556]) # input electric power of pump (kW)
q1 = np.array([1,2,3.5,4.5,5.5,6.5,7.164,7.5,8.5,8.8])*100/3600
h1 = np.array([138,137,135,130,122,113,106,102,90,86])
p1 = np.array([160,178,204,220,233,245,252,253.5,260,262])

fig,ax = plt.subplots(1,1,figsize=(3,2.4))
axright = ax.twinx()
ax.plot(q*3600,h,label='hw,h',color=cl_cloR)
ax.plot(q1*3600,h1,label='lw,h',color=cl_leaveG)
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