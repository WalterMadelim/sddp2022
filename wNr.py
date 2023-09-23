import gurobipy as gp
from gurobipy import GRB
import numpy as np

def gsum(a):
    return gp.quicksum(a)

htoulast = 100
tkCoef = 1/20
qt_pump = np.array([[1,7,11],[1,5.5,8.8]],dtype=np.float64)*100 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]],dtype=np.float64)
pt_pump = np.array([[330,500,556],[160,233,262]],dtype=np.float64) # kW

qt_pipe = np.array([100,500,1000],dtype=np.float64)/3600
ppK = np.array([1547.2653,1041.4285,833.1428,5207.1425,3644.99975,2395.28555],dtype=np.float64)# h=Kq**2, q in (m3/s), h in (m)
ht_pipe = np.zeros((len(ppK),3),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
qt_pipe = np.array([100,500,1000],dtype=np.float64)


pipesNum = 6
pumpNum = 2
hrdNum = 3
s=0

m = gp.Model('water_t')
# -----model part-----
s_pump = m.addVars(pumpNum,4,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # represent 4 finite states
q_pump,h_pump,p_pump = m.addVars(pumpNum,lb=100,ub=1100),m.addVars(pumpNum,lb=86,ub=229),m.addVars(pumpNum,lb=160,ub=556)
l_pump = m.addVars(pumpNum, ub=1)
s_pipe = m.addVars(pipesNum,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
q_pipe,h_pipe = m.addVars(pipesNum,ub=1500),m.addVars(pipesNum,ub=704)
l_pipe = m.addVars(pipesNum)
hrd = m.addVars(hrdNum) # RDV
h0, h1, h2, h4 = m.addVar(),m.addVar(),m.addVar(),m.addVar()

htin = m.addVar()
htou = m.addVar() # discretized state variable

for n in range(pumpNum):
    m.addConstr(gsum(s_pump[n, i] for i in range(4)) == 1)
    m.addConstr(l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000*(1-s_pump[n, 3]))
    m.addConstr(l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000*(1-s_pump[n, 3]))
    m.addConstr(l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143*(1-s_pump[n, 3]))
    m.addConstr(l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143*(1-s_pump[n, 3]))
    m.addConstr(l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396*(1-s_pump[n, 0]))
    m.addConstr(l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396*(1-s_pump[n, 1]))
    m.addConstr(l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396*(1-s_pump[n, 2]))
    m.addConstr(l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396*(1-s_pump[n, 3]))
    m.addConstr(l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396*(1-s_pump[n, 3]))
for n in range(pipesNum):
    m.addConstr(l_pipe[n] * qt_pipe[0] + (1 - l_pipe[n]) * qt_pipe[1] - q_pipe[n] >= -1500*s_pipe[n])
    m.addConstr(l_pipe[n] * qt_pipe[0] + (1 - l_pipe[n]) * qt_pipe[1] - q_pipe[n] <= 1500*s_pipe[n])
    m.addConstr(l_pipe[n]*qt_pipe[2]+(1-l_pipe[n])*qt_pipe[1] - q_pipe[n] >= -1500*(1-s_pipe[n]))
    m.addConstr(l_pipe[n]*qt_pipe[2]+(1-l_pipe[n])*qt_pipe[1] - q_pipe[n] <= 1500*(1-s_pipe[n]))
    m.addConstr(l_pipe[n] * ht_pipe[n, 0] + (1 - l_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] >= -704 * s_pipe[n])
    m.addConstr(l_pipe[n] * ht_pipe[n, 0] + (1 - l_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] <= 704 * s_pipe[n])
    m.addConstr(l_pipe[n]*ht_pipe[n,2]+(1-l_pipe[n])*ht_pipe[n,1] - h_pipe[n] >= -704*(1-s_pipe[n]))
    m.addConstr(l_pipe[n]*ht_pipe[n,2]+(1-l_pipe[n])*ht_pipe[n,1] - h_pipe[n] <= 704*(1-s_pipe[n]))
    # flow constrs

m.addConstr(q_pump[0] == q_pipe[2] + q_pipe[0]) # node 0
m.addConstr(q_pipe[0]+q_pump[1] == q_pipe[1]) # node 1
m.addConstr(q_pipe[1] == q_pipe[5]) # node 2
m.addConstr(q_pipe[4]+q_pipe[5] == lW) # node 3: load level: 140~1000
m.addConstr(q_pipe[4] == q_pipe[3]) # node 4
# pressure constrs
m.addConstr(h_pump[0] == h0) # node 0_with pump
m.addConstr(h_pump[1] == h1) # node 1_with pump
m.addConstr(h0-hrd[0]-h_pipe[0] == h1)# node 1
m.addConstr(h1-h_pipe[1] == h2)# node 2
m.addConstr(h2 - h_pipe[5] - hrd[1] == 0)# node 5
m.addConstr(h0-h_pipe[2] == htin)# node tank
m.addConstr(htou - h_pipe[3] == h4)# node 6
m.addConstr(h4 - h_pipe[4] - hrd[2] == 0)# node 5

m.addConstr(ptkpm == 1000 * q_pipe[2] * 9.8 * htoulast)

activePM = 1e-3*gsum(p_pump[p] for p in range(pumpNum))
theObj = 10*gsum(hrd[i] for i in range(hrdNum)) + Pe[t]*activePM*(1 + np.tan(np.arccos(.85)))
# theObj = 1e-5*Pe[t]*(q_pump[0] + q_pump[1])
if t == nss-1: # penalty for not-satisfying daily balance
    endPen = m.addVar()
    m.addConstr(endPen >= iniHtou-gdot(fbhtou,htou))
    theObj += 10*endPen
