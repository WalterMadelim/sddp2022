import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np

m = gp.Model('mypn')
x = m.addVar(lb=-GRB.INFINITY)
y = m.addVar(lb=-GRB.INFINITY)
m.addConstr(y==0)
m.setObjective(.5*(x**2-y**2)-y)
m.optimize()
print(m.ObjVal)
print(x.X,y.X)
r = m.getConstrs()
r = r[0]
print(r.Pi)
# prodMax = 35

# m = gp.Model('extensive_NA')
# y00 = m.addVar(vtype=GRB.BINARY)
# p00,n00,x00 = m.addVar(),m.addVar(),m.addVar()
# y10 = m.addVar(vtype=GRB.BINARY)
# p10,n10,x10 = m.addVar(),m.addVar(),m.addVar()
#
# y01 = m.addVar(vtype=GRB.BINARY)
# p01,n01,x01 = m.addVar(),m.addVar(),m.addVar()
# y11 = m.addVar(vtype=GRB.BINARY)
# p11,n11,x11 = m.addVar(),m.addVar(),m.addVar()
#
# cost0 = 5*y00+.5*x00+2*n00+p00 + 5*y10+.5*x10+4*n10+p10
# cost1 = 5*y01+.5*x01+2*n01+p01 + 5*y11+.5*x11+4*n11+p11
# m.setObjective(.5*cost0+.5*cost1)
#
# m.addConstr(p00 - n00 == 0)
# m.addConstr(x00 <= prodMax*y00)
# m.addConstr(p00+x00 <= 30)
# m.addConstr(p10 - n10 == -10+p00+x00-n00)
# m.addConstr(x10 <= prodMax*y10)
# m.addConstr(p10+x10 <= 30)
#
# m.addConstr(p01 - n01 == 0)
# m.addConstr(x01 <= prodMax*y01)
# m.addConstr(p01+x01 <= 30)
# m.addConstr(p11 - n11 == -30+p01+x01-n01)
# m.addConstr(x11 <= prodMax*y11)
# m.addConstr(p11+x11 <= 30)
#
# # NA constraint
# m.addConstr(x00 == 0.5*x00+0.5*x01)
# m.addConstr(x01 == 0.5*x00+0.5*x01)
#
# m.optimize()
# if m.status != GRB.OPTIMAL:
#     print('opt Fail >>>>>>>>>>>>>')
#     exit(3)
#
# print('cost =',m.ObjVal)
# # cost0 = 5*y0.X+.5*x0.X+2*n0.X+p0.X
# # cost10 = 5*y10.X+.5*x10.X+4*n10.X+p10.X
# # cost11 = 5*y11.X+.5*x11.X+4*n11.X+p11.X
# # print(cost0,cost10,cost11)
# #
# # print(y0.X,y10.X,y11.X)
# # print(p0.X,p10.X,p11.X)
# # print(n0.X,n10.X,n11.X)
# # print(x0.X,x10.X,x11.X)