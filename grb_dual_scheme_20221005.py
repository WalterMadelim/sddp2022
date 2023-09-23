import numpy as np
import gurobipy as gp
from gurobipy import GRB
# gurobi dual scheme
# 20221005
x1 = 0
x2 = 0
m = gp.Model("MC-SDDiP")
y = m.addVar(ub=4)
m.addConstr(y>=2.6-x1/4-x2/2)
m.setObjective(4*y)
m.optimize()
l = m.getConstrs()[0].Pi
print(m.ObjVal)
print(l)
print('Pi_x1: %g' % (-l/4))
print('Pi_x2: %g' % (-l/2))
print('Cst: %g' % (m.ObjVal - (-l/4) * x1 - (-l/2) * x2))

# constr: tha >= Pi_x1 * x1 + Pi_x2 * x2 + Cst

