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
exit(5)

qt_pump = np.array([[1,7,11],[1,5.5,8.8]])*100/3600 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]])
pt_pump = np.array([[330,500,556],[160,233,262]])
ppK = np.array([1547.2653,1041.4285,833.1428,5207.1425,3644.99975,2395.28555])# h=Kq**2, q in (m3/s), h in (m)
qt_pipe = np.array([100,200,500,1000])/3600
ht_pipe = np.zeros((len(ppK),4),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
pipesNum = 6
pumpNum = 2

# load at node J5: range: 300~1300 (m3/h)
T = 2 # dispatch periods
m = gp.Model('pwl_of_pump(single time)')
s_pump = m.addVars(T, pumpNum, 4, vtype=GRB.BINARY) # represent 4 finite states
q_pump,h_pump,p_pump = m.addVars(T,pumpNum),m.addVars(T,pumpNum),m.addVars(T,pumpNum)
l_pump = m.addVars(T, pumpNum, ub=1)
for t in range(T):
    for n in range(pumpNum):
        m.addConstr(gp.quicksum(s_pump[t, n, i] for i in range(4)) == 1)
        m.addConstr((s_pump[t, n, 0] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[0, 0] + (1 - l_pump[t, n]) * qt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 0] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[0, 0] + (1 - l_pump[t, n]) * ht_pump[0, 1]))
        m.addConstr((s_pump[t, n, 0] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[0, 0] + (1 - l_pump[t, n]) * pt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[0, 2] + (1 - l_pump[t, n]) * qt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[0, 2] + (1 - l_pump[t, n]) * ht_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[0, 2] + (1 - l_pump[t, n]) * pt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[1, 0] + (1 - l_pump[t, n]) * qt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[1, 0] + (1 - l_pump[t, n]) * ht_pump[1, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[1, 0] + (1 - l_pump[t, n]) * pt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[1, 2] + (1 - l_pump[t, n]) * qt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[1, 2] + (1 - l_pump[t, n]) * ht_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[1, 2] + (1 - l_pump[t, n]) * pt_pump[1, 1]))
s_pipe = m.addVars(T,pipesNum,3,vtype=GRB.BINARY) # represent 3 finite states
q_pipe,h_pipe = m.addVars(T,pipesNum),m.addVars(T,pipesNum)
l_pipe = m.addVars(T,pipesNum)
for t in range(T):
    for n in range(pipesNum):
        m.addConstr(gp.quicksum(s_pipe[t, n, i] for i in range(3)) == 1)
        m.addConstr((s_pipe[t,n,0]==1)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[0]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n,0]==1)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,0]+(1-l_pipe[t,n])*ht_pipe[n,1]))
        m.addConstr((s_pipe[t,n,2]==1)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[3]+(1-l_pipe[t,n])*qt_pipe[2]))
        m.addConstr((s_pipe[t,n,2]==1)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,3]+(1-l_pipe[t,n])*ht_pipe[n,2]))
        m.addConstr((s_pipe[t,n,1]==1)>>(l_pipe[t,n]<=1))
        m.addConstr((s_pipe[t,n,1]==1)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[2]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n,1]==1)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,2]+(1-l_pipe[t,n])*ht_pipe[n,1]))
hrd = m.addVars(T,3)
h0,h1,h2,h5,h6 = m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T)
htank = m.addVars(T)

t = 0
# flow constrs
m.addConstr(q_pump[t,0] == q_pipe[t,2] + q_pipe[t,0]) # node 0
# m.addConstr(q_pipe[t,0] == q_pipe[t,1]) # node 1
m.addConstr(q_pipe[t,0]+q_pump[t,1] == q_pipe[t,1]) # node 1
m.addConstr(q_pipe[t,1] == q_pipe[t,5]) # node 2
m.addConstr(q_pipe[t,4]+q_pipe[t,5] == 400/3600) # node 5: load level: 140~1000
m.addConstr(q_pipe[t,2] == q_pipe[t,3]) # node tank
m.addConstr(q_pipe[t,4] == q_pipe[t,3]) # node 6

# pressure constrs
m.addConstr(h_pump[t,0] == h0[t]) # node 0_with pump
m.addConstr(h_pump[t,1] == h1[t]) # node 1_with pump
m.addConstr(h0[t]-hrd[t,0]-h_pipe[t,0] == h1[t])# node 1
m.addConstr(h1[t]-h_pipe[t,1] == h2[t])# node 2
m.addConstr(h2[t]-h_pipe[t,5]-hrd[t,1] == h5[t])# node 5
m.addConstr(h0[t]-h_pipe[t,2] == htank[t])# node tank
m.addConstr(htank[t]-h_pipe[t,3] == h6[t])# node 6
m.addConstr(h6[t]-h_pipe[t,4]-hrd[t,2] == h5[t])# node 5

m.setObjective(-p_pump[t,0])
m.optimize()
if m.status != GRB.OPTIMAL:
    print('opt fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
    sys.exit(3)

# print(hrd[t,0].X,hrd[t,1].X,hrd[t,2].X)
# print(s_pump[t,0,0].X,s_pump[t,0,1].X,s_pump[t,0,2].X,s_pump[t,0,3].X)
# print(s_pump[t,1,0].X,s_pump[t,1,1].X,s_pump[t,1,2].X,s_pump[t,1,3].X)
# print(h0[t].X,h1[t].X,h2[t].X,h5[t].X)
# print(q_pipe[t,1].X)
# print(h0[t].X,htank[t].X,h6[t].X,h5[t].X)
# print(q_pipe[t,4].X)
# print(q_pipe[t,0].X)
