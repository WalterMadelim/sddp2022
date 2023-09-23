import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np
# xb2,xb1 = 1,0 # x = 2
xb21 = np.array([1,0])
# m = gp.Model("subproblem")
# y = m.addVar(vtype=GRB.INTEGER,ub=2)
# m.addConstr(1.5*y >= 2*xb2+xb1)
# m.setObjective(y)
# m.setParam('OutputFlag', 0)
# m.optimize()
#

pai21 = 10*np.array([-.4,.5]) # initial pai
for i in range(18):
    m = gp.Model("Lag-sub")
    y = m.addVar(vtype=GRB.INTEGER,ub=2)
    z2,z1 = m.addVar(ub=1),m.addVar(ub=1)
    m.addConstr(1.5*y >= 2*z2+z1)
    m.setObjective(y-pai21[0]*z2-pai21[1]*z1)
    m.setParam('OutputFlag', 0)
    m.optimize()
    Obj1st = m.ObjVal+pai21[0]*xb21[0]+pai21[1]*xb21[1]
    cpv = np.zeros(2)
    cpv[0],cpv[1] = z2.X,z1.X
    drc = xb21-cpv
    print('f(%8g,%8g)=%8g' % (pai21[0],pai21[1],Obj1st))
    if sum(abs(drc))<1e-4:
        print('subgrad = 0, opt find')
        break
    pai21 += drc

