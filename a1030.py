import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
# extensive form
prodMax = 35

pai = np.array([0.5,-2.0])
s = 0

m = gp.Model('extensive')
vtype = GRB.BINARY if s != 1 else GRB.CONTINUOUS
y0,y10,y11 = m.addVar(vtype=vtype,ub=1),m.addVar(vtype=vtype,ub=1),m.addVar(vtype=vtype,ub=1)
x0,x10,x11 = m.addVar(),m.addVar(),m.addVar() # only x0 is useful
n0,n10,n11 = m.addVar(),m.addVar(),m.addVar()
p0,p10,p11 = m.addVar(),m.addVar(),m.addVar()
cost0 = 5*y0+.5*x0+2*n0+p0
cost10 = 5*y10+.5*x10+4*n10+p10
cost11 = 5*y11+.5*x11+4*n11+p11
if s != 2:
    m.setObjective(cost0+(cost10+cost11)/2)
else:
    m.setObjective(cost0+(cost10+cost11)/2+pai[0]*(-10+p0+x0-n0+n10-p10)+pai[1]*(-30+p0+x0-n0+n11-p11))
if s != 2:
    m.addConstr(p10 - n10 == -10+p0+x0-n0) # link
    m.addConstr(p11 - n11 == -30+p0+x0-n0) # link
# stage0 local
m.addConstr(p0 - n0 == 0)
m.addConstr(x0 <= prodMax*y0)
m.addConstr(p0+x0 <= 30)
# stage1 local: per scene
m.addConstr(x10 <= prodMax*y10)
m.addConstr(p10+x10 <= 30)

m.addConstr(x11 <= prodMax*y11)
m.addConstr(p11+x11 <= 30)

m.setParam('OutputFlag',0)
m.optimize()
if m.status != GRB.OPTIMAL:
    print('opt Fail >>>>>>>>>>>>>')
    exit(3)
if s == 1:
    l = m.getConstrs()
    print(l[0].Pi, l[1].Pi)

print('cost =',m.ObjVal)
# cost0 = 5*y0.X+.5*x0.X+2*n0.X+p0.X
# cost10 = 5*y10.X+.5*x10.X+4*n10.X+p10.X
# cost11 = 5*y11.X+.5*x11.X+4*n11.X+p11.X
# print(cost0,cost10,cost11)

print(y0.X,y10.X,y11.X)
print(p0.X,p10.X,p11.X)
print(n0.X,n10.X,n11.X)
print(x0.X,x10.X,x11.X)