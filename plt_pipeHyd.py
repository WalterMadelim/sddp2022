import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
from matplotlib import pyplot as plt
# from simulink experimental data for hydro character of a pipe
# 20220923

cl_cloR   = 138/256,33 /256,66/256
cl_lotusR = 207/256,131/256,185/256
cl_leaveG = 96 /256,124/256,56/256
cl_fenceG = 66 /256,100/256,93/256
cl_hairB  = 32 /256,25/256,25/256

ppl = np.array([400,1300,3700,1000,4300,2600,800,5000,3500,2300])
ppd = np.array([250,250,300,250,300,250,250,250,250,250])/1000
ppK = np.array([416.5714,1353.85705,1547.2653,1041.4285,1798.1732,2707.7141,833.1428,5207.1425,3644.99975,2395.28555])

ppA = np.pi*(ppd/2)**2
print(ppd)

axLabelFont = 11
# for the pipe0: collect experimental data
# q = np.array([300])/3600 # m3/s
# ploss = np.array([3.258e4,]) # pa
q = np.array([1,10,50,100,200,300,400,500,600,700,800,900,1000])/3600 # m3/s
hloss = np.array([1.164,72.39,1244,4341,1.54e4,3.258e4,5.569e4,8.461e4,1.193e5,1.597e5,2.058e5,2.575e5,3.15e5])/1000/9.8 # m
q1 = np.array([1,200,400,700,1000])/3600 # m3/s
hloss1 = np.array([1.164,1.54e4,5.569e4,1.597e5,3.15e5])/1000/9.8 # m

ppK = np.array([416.5714,1353.85705,1547.2653,1041.4285,1798.1732,2707.7141,833.1428,5207.1425,3644.99975,2395.28555])


fig,ax = plt.subplots(1,1,figsize=(4,2.5))
ax.plot(q*3600,hloss,color = cl_cloR,linewidth=3,label='EXP')
ax.plot(q1*3600,hloss1,color = cl_hairB,linewidth = 1, marker = 'o',label='PWL')
ax.set_xlabel('Volumetric flow rate ($m^3/h$)',fontsize=axLabelFont)
ax.set_ylabel('Head loss (m)',fontsize=axLabelFont)
ax.legend()
fig.tight_layout()
plt.show()
exit(5)
