import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np

# Cuts in SDDIP
# 20221008

def dotgp(a,b,zl):
    return gp.quicksum(a[i]*b[i] for i in range(zl))
def getCpv(z,zl):
    a = np.zeros(zl)
    for i in range(zl):
        a[i] = z[i].X
    return a
def getPai(l,zl):
    a = np.zeros(zl)
    for i in range(zl):
        a[i] = l[i].Pi
    return a
def tmd(s,Al,pai): # triplet model, pai is used only when s == 2
    # s = 0: primal-sub(MILP); 1: LP-relaxed-sub(return pai); 2: copy-relaxed-MILP(return val and z)
    zl = len(Al) # length of copy vector = length of lastAction
    m = gp.Model('triple')
    z = m.addVars(zl,ub=1)
    if s!=2:
        for i in range(zl):
            m.addConstr(z[i]==Al[i])
    # -----model part-----
    y = m.addVar(vtype=GRB.CONTINUOUS if s==1 else GRB.INTEGER,ub=2)
    m.addConstr(1.5*y >= 2*z[0]+z[1])
    theObj = y
    # -----model part-----
    if s == 2:
        theObj -= dotgp(pai,z,zl)
    m.setObjective(theObj)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if s == 0:
        return m.ObjVal # the std value used to check tightness
    elif s == 1:
        return getPai(m.getConstrs(),zl) # initial pai for lag_cut, pai for SB and B_cut
    elif s == 2:
        subVal = m.ObjVal
        tailVal = np.dot(pai,Al)
        maxerVal = subVal + tailVal
        cpv = getCpv(z,zl)
        drc = Al - cpv
        return maxerVal,subVal,tailVal,drc,cpv
Al = np.array([1,0]) # last action is 2
print('std Val = %8g' % (tmd(0,Al,0)))
pai = tmd(1,Al,0)
print('ini_pai = ',pai)
while 1:
    maxerVal,subVal,tailVal,drc,cpv = tmd(2,Al,pai)
    print('\nf(%8g,%8g)=%8g=%8g+%8g' % (pai[0],pai[1],maxerVal,subVal,tailVal))
    print('drc is:',drc,'cpv is:',cpv)
    if sum(abs(drc))<1e-4:
        print('subGrad = 0, opt Attained!')
        break
    pai += drc


# pai21 = pai # initial pai
# for i in range(18):
#     m = gp.Model("Lag-sub")
#     y = m.addVar(vtype=GRB.INTEGER,ub=2)
#     # z2,z1 = m.addVar(ub=1),m.addVar(ub=1)
#     z = m.addVars(2,ub=1)
#     m.addConstr(1.5*y >= 2*z[0]+z[1]) # use copyvec rather than lastAct
#     m.setObjective(y-dotgp(pai21,z,2))
#     m.setParam('OutputFlag', 0)
#     m.optimize()
#     Obj1st = m.ObjVal+pai21[0]*xb21[0]+pai21[1]*xb21[1]
#     cpv = getCpv(z,2)
#     drc = xb21-cpv
#     print('f(%8g,%8g)=%8g=%8g+%8g' % (pai21[0],pai21[1],Obj1st,m.ObjVal,np.dot(pai21,xb21)))
#     if sum(abs(drc))<1e-4:
#         print('subgrad = 0, opt find')
#         break
#     pai21 += drc

