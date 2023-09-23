import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np

# 2-stage-MILP with 1st stage-Binary chaining variables
# use Lag-Cuts
# 20221008

def fbgen(lb,ub): # float basis generation
    fb = []
    for i in range(40):
        fb.append(2 ** i)
    tmp = np.array(fb) / 2 ** (40//2+1)
    tmp = tmp[np.where(tmp <= ub)]
    tmp = tmp[np.where(tmp >= lb)]
    return tmp,len(tmp),tmp[0]
def mb2f(b): # b2f in gurobi model
    return gp.quicksum(fb[i]*b[i] for i in range(B))
def b2f(b): # binary vector to a float number
    return sum(fb[i]*b[i] for i in range(B))
def f2b(f): # float to binary List
    m = gp.Model('f2bv')
    x = m.addVars(B, vtype=GRB.BINARY)
    d = m.addVar(lb=-GRB.INFINITY)
    o = m.addVar()
    m.addConstr(d == mb2f(x) - f)
    m.addGenConstrAbs(o, d)
    m.setObjective(o)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        exit(3)
    xthis = np.zeros(B,dtype=np.int8)
    for i in range(B):
        xthis[i] = x[i].X
    return xthis
def xb2f(b): # binary vector to a float number
    return sum(fb[i]*b[i].X for i in range(B))
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
def tmd(dR,ite,t,s,Al,pai): # triplet model, pai is used only when s == 2
    # s = 0: primal-sub(MILP); 1: LP-relaxed-sub(return pai); 2: copy-relaxed-MILP(return val and z)
    zl = len(Al) # length of copy vector = length of lastAction
    m = gp.Model('triple')
    z = m.addVars(zl,ub=1)
    if s!=2:
        for i in range(zl):
            m.addConstr(z[i]==Al[i])
    # -----model part-----
    y = m.addVar(vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1)
    p,n,x = m.addVar(),m.addVar(),m.addVar()
    m.addConstr(x <= 7*y)
    m.addConstr(gp.quicksum(fb[i]*z[i] for i in range(B)) - gp.quicksum(fb[i]*z[B+i] for i in range(B)) - 10 == p - n - x)
    theObj = .5*y+.1*x+.2*p+.4*n
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
def tmd1st(dR,ite,t,s,Al,pai): # 1st_stage: need only to let s == 0
    zl = len(Al) # length of copy vector = length of lastAction
    m = gp.Model('1st-stage-MILP')
    z = m.addVars(zl,ub=1)
    if s!=2:
        for i in range(zl):
            m.addConstr(z[i]==Al[i])
    # -----model part-----
    tha = m.addVar()
    y = m.addVar(vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1)
    p,n,x = m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY),m.addVar()
    m.addConstr(x <= 7*y)
    m.addConstr(gp.quicksum(fb[i]*z[i] for i in range(B)) - gp.quicksum(fb[i]*z[B+i] for i in range(B)) - 5 == mb2f(p) - mb2f(n) - x)
    theObj = .5*y+.1*x+.2*mb2f(p)+.4*mb2f(n)+tha
    # -----model part-----
    if s == 2:
        theObj -= dotgp(pai,z,zl)
    if t<nss:# if not the final stage:
        for its in range(1,ite+dR):
            m.addConstr(tha>=cut[its,-1]+gp.quicksum(cut[its,i]*p[i] for i in range(B))+gp.quicksum(cut[its,B+i]*n[i] for i in range(B)))
    m.setObjective(theObj)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if s == 0:
        ub = .5 * y.X + .1 * x.X + .2 * xb2f(p) + .4 * xb2f(n)
        lb = ub + tha.X
        return ub,np.append(getCpv(p,B),getCpv(n,B)),lb # the std value used to check tightness
    elif s == 1:
        return getPai(m.getConstrs(),zl) # initial pai for lag_cut, pai for SB and B_cut
    elif s == 2:
        subVal = m.ObjVal
        tailVal = np.dot(pai,Al)
        maxerVal = subVal + tailVal
        cpv = getCpv(z,zl)
        drc = Al - cpv
        return maxerVal,subVal,tailVal,drc,cpv




fb,B,SCL = fbgen(2**-2,2**4)
MaxIte = 8
nss = 2
cut = np.zeros((MaxIte,2*B+1))
p0,n0 = 2,0
p0 = f2b(p0)
n0 = f2b(n0)
Al0 = np.append(p0,n0)
for ite in range(1,MaxIte):
    t=1
    val1,Al1,lb = tmd1st(0,ite,t,0,Al0,0)
    t=2
    val2 = tmd(1,ite,t,0,Al1,0)
    ub = val1 + val2
    print('ub = %8g + %8g = %8g|%8g = lb' % (val1,val2,ub,lb))
    pai = tmd(1,ite,t,1,Al1,0) # print('ini_pai = ',pai)
    while 1:
        maxerVal,subVal,tailVal,drc,cpv = tmd(1,ite,t,2,Al1,pai)
        # print('\nf(pai)=%8g=%8g+%8g' % (maxerVal,subVal,tailVal))
        # print('drc is:',drc,'\ncpv is:',cpv)
        if sum(abs(drc))<SCL/2**10:
            # print('subGrad = 0, opt Attained!')
            cut[ite] = np.append(pai, subVal)
            break
        pai += drc
        # print('pai is:',pai)
    # print(cut[ite])
