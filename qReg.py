import xlrd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import sys
N = 80
sh = xlrd.open_workbook(r'I:\xlsx\temps.xls').sheet_by_name("sheet1")
dset = np.empty((24,N),dtype=np.float64)
for c in range(1,N+1):
    for r in range(1,25):
        dset[r-1,c-1] = sh.cell_value(r,c)
btree = np.empty((24,2),dtype=np.float64) # q_1/3,q_2/3
q13 = 1/3
q23 = 2/3
for j in range(1,24):
    x = dset[j-1,:]
    y = dset[j,:]
    m = gp.Model('qRegression')
    b = m.addVar(lb=-GRB.INFINITY)
    U = m.addVars(N,lb=-GRB.INFINITY)
    ind = m.addVars(N,vtype=GRB.BINARY)
    m.addConstrs((ind[i]==1)>>(U[i]<=0) for i in range(N))
    m.addConstrs((ind[i]==0)>>(U[i]>=0) for i in range(N))
    m.addConstrs(U[i]==y[i]-x[i]*b for i in range(N))
    m.setObjective(gp.quicksum((q13-ind[i])*U[i] for i in range(N)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>> Q(x,xi) eval problem: No solution!')
        sys.exit(3)
    # print('objval is %g, q25_b is %g' % (m.ObjVal,b.X))
    btree[j,0] = b.X
    m.setObjective(gp.quicksum((q23-ind[i])*U[i] for i in range(N)))
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>> Q(x,xi) eval problem: No solution!')
        sys.exit(3)
    # print('objval is %g, q50_b is %g' % (m.ObjVal,b.X))
    btree[j,1] = b.X
print(btree)
cnt1,cnt2 = 0,0
for i in range(N):
    pass