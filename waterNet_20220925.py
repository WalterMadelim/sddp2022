import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
from matplotlib import pyplot as plt

# water distribution network with 2 pumps, 3 hrds, 1 load and 1 tank
# doorvanbei
# 20220925

# waterNetwork:
# load at node J5: range: 300~1300 (m3/h)
wLoad = np.array([368,395,395,401,380,358,357,319,314,343,357,335,334,324,316,319,331,341,338,331,398,422,385,352])
vP,fP,pP = 308.9,626.8,1044
eCost = np.array([vP,vP,vP,vP,vP,vP,vP,vP,pP,pP,pP,pP,fP,fP,fP,fP,fP,pP,pP,pP,pP,fP,fP,fP])
tkCoef = 3600/20
# tabuleu data
qt_pump = np.array([[1,7,11],[1,5.5,8.8]])*100/3600 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]])
pt_pump = np.array([[330,500,556],[160,233,262]]) # kW
qt_pipe = np.array([100,500,1000])/3600
ppK = np.array([1547.2653,1041.4285,833.1428,5207.1425,3644.99975,2395.28555])# h=Kq**2, q in (m3/s), h in (m)
ht_pipe = np.zeros((len(ppK),3),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
pipesNum = 6
pumpNum = 2

T = 24 # dispatch periods
# BASIC DAILY PLOT
# hTank,wLoad,eCost
# fig,ax = plt.subplots(1,1)
# ax.step(np.arange(T),hTank/100)
# ax.step(np.arange(T),wLoad/400)
# ax.step(np.arange(T),eCost/600)
# plt.show()
# exit(4)

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
s_pipe = m.addVars(T,pipesNum,vtype=GRB.BINARY) # represent 3 finite states
q_pipe,h_pipe = m.addVars(T,pipesNum),m.addVars(T,pipesNum)
l_pipe = m.addVars(T,pipesNum)
for t in range(T):
    for n in range(pipesNum):
        m.addConstr((s_pipe[t,n]==0)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[0]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n]==0)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,0]+(1-l_pipe[t,n])*ht_pipe[n,1]))
        m.addConstr((s_pipe[t,n]==1)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[2]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n]==1)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,2]+(1-l_pipe[t,n])*ht_pipe[n,1]))

hrd = m.addVars(T,3)
h0,h1,h2,h5,h6 = m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T)
htin,htou = m.addVars(T),m.addVars(T)

m.addConstr(htou[0] == 100 + tkCoef * (q_pipe[0, 2] - q_pipe[0, 3]))
m.addConstrs(htou[t] == htou[t-1] + tkCoef * (q_pipe[t, 2] - q_pipe[t, 3]) for t in range(1, T))
m.addConstr(htou[T-1] >= 80)
for t in range(T):
    # flow constrs
    m.addConstr(q_pump[t,0] == q_pipe[t,2] + q_pipe[t,0]) # node 0
    m.addConstr(q_pipe[t,0]+q_pump[t,1] == q_pipe[t,1]) # node 1
    m.addConstr(q_pipe[t,1] == q_pipe[t,5]) # node 2
    m.addConstr(q_pipe[t,4]+q_pipe[t,5] == wLoad[t]/3600) # node 5: load level: 140~1000
    m.addConstr(q_pipe[t,4] == q_pipe[t,3]) # node 6
    # pressure constrs
    m.addConstr(h_pump[t,0] == h0[t]) # node 0_with pump
    m.addConstr(h_pump[t,1] == h1[t]) # node 1_with pump
    m.addConstr(h0[t]-hrd[t,0]-h_pipe[t,0] == h1[t])# node 1
    m.addConstr(h1[t]-h_pipe[t,1] == h2[t])# node 2
    m.addConstr(h2[t]-h_pipe[t,5]-hrd[t,1] == h5[t])# node 5
    m.addConstr(h0[t]-h_pipe[t,2] == htin[t])# node tank
    m.addConstr(htou[t]-h_pipe[t,3] == h6[t])# node 6
    m.addConstr(h6[t]-h_pipe[t,4]-hrd[t,2] == h5[t])# node 5
m.setObjective( 10*gp.quicksum(hrd[t,i] for i in range(3) for t in range(T)) + gp.quicksum(eCost[t] * p_pump[t,p] for p in range(2) for t in range(T)))

m.optimize()
if m.status != GRB.OPTIMAL:
    print('opt fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
    sys.exit(3)
print('\nThe ObjVal is %g' % m.ObjVal)
print('\nThe tank water level')
for t in range(T):
    print('%8d|%8g' % (t,htou[t].X))
print('\nThe lower|upper flow')
for t in range(T):
    print('%8d|%8g|%8g' % (t,3600*q_pipe[t,4].X,3600*q_pipe[t,5].X))
print('\nThe pump0|pump1 w level')
for t in range(T):
    print('%8d:%3d|%3d|%3d|%3d| |%3d|%3d|%3d|%3d' % (t,s_pump[t,0,0].X,s_pump[t,0,1].X,s_pump[t,0,2].X,s_pump[t,0,3].X,s_pump[t,1,0].X,s_pump[t,1,1].X,s_pump[t,1,2].X,s_pump[t,1,3].X))
    print('%8d:%8g|%8g (kW)' % (t,p_pump[t,0].X,p_pump[t,1].X))
    print('%8d:%8g|%8g (m)' % (t,h_pump[t,0].X,h_pump[t,1].X))
    print('%8d:%8g|%8g (m3/h)' % (t,3600*q_pump[t,0].X,3600*q_pump[t,1].X))
print('\nThe hrd')
for t in range(T):
    print('%8d|%12g|%12g|%12g' % (t,hrd[t,0].X,hrd[t,1].X,hrd[t,2].X))

