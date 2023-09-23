import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB

# single line supply waterNet
# doorvanbei

def gsum(a):
    return gp.quicksum(a)
lup = .77
ppKf = 1.016 # ppK = ppKf * (km) * (m) ** -5, length, diameter
l = np.array([3.7,lup,.521,.521,.521,lup])
d = np.array([.3,.25,.25,.25,.25,.25])
ppK = ppKf * l * d ** -5
ppK = ppK[[1,5]]
qt_pipe = np.array([0,100,500,1100],dtype=np.float64)/3600
ht_pipe = np.zeros((len(ppK),4),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
qt_pipe = np.array([0,100,500,1100],dtype=np.float64)

qt_pump = np.array([[1,7,11],[1,5.5,8.8]],dtype=np.float64)*100 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]],dtype=np.float64)
pt_pump = np.array([[330,500,556],[160,233,262]],dtype=np.float64) # kW

m = gp.Model('small_water')
pipeNum = 2
pumpNum = 1
s_pipe = m.addVars(pipeNum,3,vtype=GRB.BINARY)
q_pipe = m.addVars(pipeNum,ub=1500)
h_pipe = m.addVars(pipeNum,ub=704)
l0_pipe = m.addVars(pipeNum,ub=1)
l1_pipe = m.addVars(pipeNum)
for n in range(pipeNum):
    m.addConstr(gsum(s_pipe[n,i] for i in range(3)) == 1)
    m.addConstr(l0_pipe[n] * qt_pipe[0] + (1 - l0_pipe[n]) * qt_pipe[1] - q_pipe[n] >= -1500 * (1-s_pipe[n,0]))
    m.addConstr(l0_pipe[n] * qt_pipe[0] + (1 - l0_pipe[n]) * qt_pipe[1] - q_pipe[n] <= 1500 * (1-s_pipe[n,0]))
    m.addConstr(l0_pipe[n] * ht_pipe[n,0] + (1 - l0_pipe[n]) * ht_pipe[n,1] - h_pipe[n] >= -704 * (1-s_pipe[n,0]))
    m.addConstr(l0_pipe[n] * ht_pipe[n,0] + (1 - l0_pipe[n]) * ht_pipe[n,1] - h_pipe[n] <= 704 * (1-s_pipe[n,0]))
    m.addConstr(l0_pipe[n] * qt_pipe[1] + (1 - l0_pipe[n]) * qt_pipe[2] - q_pipe[n] >= -1500 * (1-s_pipe[n,1]))
    m.addConstr(l0_pipe[n] * qt_pipe[1] + (1 - l0_pipe[n]) * qt_pipe[2] - q_pipe[n] <= 1500 * (1-s_pipe[n,1]))
    m.addConstr(l0_pipe[n] * ht_pipe[n,1] + (1 - l0_pipe[n]) * ht_pipe[n,2] - h_pipe[n] >= -704 * (1-s_pipe[n,1]))
    m.addConstr(l0_pipe[n] * ht_pipe[n,1] + (1 - l0_pipe[n]) * ht_pipe[n,2] - h_pipe[n] <= 704 * (1-s_pipe[n,1]))
    m.addConstr(l1_pipe[n] * qt_pipe[3] + (1 - l1_pipe[n]) * qt_pipe[2] - q_pipe[n] >= -1500 * (1-s_pipe[n,2]))
    m.addConstr(l1_pipe[n] * qt_pipe[3] + (1 - l1_pipe[n]) * qt_pipe[2] - q_pipe[n] <= 1500 * (1-s_pipe[n,2]))
    m.addConstr(l1_pipe[n] * ht_pipe[n,3] + (1 - l1_pipe[n]) * ht_pipe[n,2] - h_pipe[n] >= -704 * (1-s_pipe[n,2]))
    m.addConstr(l1_pipe[n] * ht_pipe[n,3] + (1 - l1_pipe[n]) * ht_pipe[n,2] - h_pipe[n] <= 704 * (1-s_pipe[n,2]))
s_pump = m.addVars(pumpNum,4,vtype=GRB.BINARY) # represent 4 finite states, the former 2 represent high speed
q_pump,h_pump,p_pump = m.addVars(pumpNum,lb=100,ub=1100),m.addVars(pumpNum,lb=86,ub=229),m.addVars(pumpNum,lb=160,ub=556)
l_pump = m.addVars(pumpNum, ub=1)
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
m.addConstr(q_pipe[0] == q_pipe[1])
m.addConstr(q_pipe[0] == q_pump[0])
m.addConstr(q_pipe[0] == 100) # demand lb = 100, ub = 1098,
rdv = m.addVar()
# m.addConstr(rdv == 0)
m.addConstr(h_pump[0]-h_pipe[1]-h_pipe[0]-rdv == 0)
m.setObjective(rdv)
m.optimize()
if m.status != GRB.OPTIMAL:
    print('opt Fail >>>>>>>>>>>>>')
    exit(3)
print('head redeced by valve %8g' % (rdv.X))
print('pump head gain: %8g' % (h_pump[0].X))
print('flow rate: %8g' % (q_pipe[0].X))
print('pump state: %8d|%8d|%8d|%8d' % (s_pump[0,0].X,s_pump[0,1].X,s_pump[0,2].X,s_pump[0,3].X))
print('pipe0:')
print('head loss: %8g' % (h_pipe[0].X))
print('pipe1:')
print('head loss: %8g' % (h_pipe[1].X))




