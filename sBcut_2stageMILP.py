import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np

# 2-stage-MILP with 1st stage-Binary chaining variables
# use sB-cut
# 20221008
def gsum(a):
    return gp.quicksum(a)
def fbgen(lb,ub):
    fb = np.zeros(40)
    for i in range(-20,20):
        fb[i+20] = 2**i
    fb = fb[(fb<ub)&(fb>lb)]
    return fb,len(fb),fb[0]
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
    return gsum(a[i]*b[i] for i in range(zl))
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
def tmd(dR,ite,t,Al,s,pai): # triplet model, pai is used only when s == 2
    # s = 0: primal-sub(MILP); 1: LP-relaxed-sub(return pai); 2: copy-relaxed-MILP(return val and z)
    m = gp.Model('triple')
    if t != 0:
        zl = len(Al) # length of copy vector = length of lastAction
        z = m.addVars(zl,ub=1)
        if s!=2:
            m.addConstrs(z[i]==Al[i] for i in range(zl))
    # -----model part-----
    y = m.addVar(vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1) # intrinsict binary
    Y,p,n,x = m.addVars(B,vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1),m.addVars(B,vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1),m.addVars(B,vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1),m.addVars(B,vtype=GRB.CONTINUOUS if s==1 else GRB.BINARY,ub=1) # Y is intrinsict
    o,D = m.addVar(ub=2.0175),m.addVar() # 100=base(o,p,n,x,D)
    if t == 0:
        m.addConstr(.8 * xi[0,t] + .2 * Yt[t-1] - dotgp(fb,Y,B) == [-SCL / 2, SCL / 2])
        m.addConstr(dotgp(fb,n,B) - dotgp(fb,p,B) - D - nt[t-1] + pt[t-1] + xt[t-1] == [-SCL / 2, SCL / 2])
    else:
        m.addConstr(.8 * xi[0,t] + .2 * gsum(fb[i] * z[i] for i in range(B)) - dotgp(fb,Y,B) == [-SCL / 2, SCL / 2])
        m.addConstr(dotgp(fb,n,B) - dotgp(fb,p,B) - D - gsum(fb[i] * z[2*B+i] for i in range(B)) + gsum(fb[i] * z[B+i] for i in range(B)) + gsum(fb[i] * z[3*B+i] for i in range(B)) == [-SCL / 2, SCL / 2])
    m.addConstr(D - .6 * expDem[t] * dotgp(fb,Y,B) == .4 * xi[1,t])
    m.addConstr(dotgp(fb,x,B) - o + .224 * y <= 8.07)
    m.addConstr(dotgp(fb,x,B) <= 5.38 * y)
    m.addConstr(dotgp(fb,x,B) + dotgp(fb,p,B) <= 9)
    m.addConstr(dotgp(fb,p,B) <= 9)

    c_end = 1 if t != nss-1 else 5
    theObj = o+.15*dotgp(fb,p,B)+c_end*.3*dotgp(fb,n,B)+2*y # base = 100

    if t != nss-1: # have tail functions
        tha = m.addVar()
        for its in range(ite+dR):
            m.addConstr(tha>=cut[its,-1]+gsum(cut[its,i]*Y[i] for i in range(B))+gsum(cut[its,B+i]*p[i] for i in range(B))+gsum(cut[its,2*B+i]*n[i] for i in range(B))+gsum(cut[its,3*B+i]*x[i] for i in range(B)))
        theObj += tha
    # -----model part-----
    if s == 2:
        theObj -= dotgp(pai,z,zl)
    m.setObjective(theObj)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        exit(3)
    if s == 0:
        stageCost = o.X+.15*xb2f(p)+c_end*.3*xb2f(n)+2*y.X
        lb = m.ObjVal if t == 0 else GRB.INFINITY # if t == 1
        return stageCost, np.concatenate((getCpv(Y,B),getCpv(p,B),getCpv(n,B),getCpv(x,B))), lb  # the std value used to check tightness
    elif s == 1:
        return getPai(m.getConstrs(),zl) # initial pai for lag_cut, pai for SB and B_cut
    elif s == 2:
        subVal = m.ObjVal
        tailVal = np.dot(pai,Al)
        maxerVal = subVal + tailVal
        cpv = getCpv(z,zl)
        drc = Al - cpv
        return maxerVal,subVal,tailVal,drc,cpv

def prtAct(Al):
    Y = b2f(Al[:B])
    p = b2f(Al[B:2*B])
    n = b2f(Al[2*B:3*B])
    x = b2f(Al[3*B:])
    print('Y:%8g|p:%8g|n:%8g|x:%8g' % (Y,p,n,x))

fb,B,SCL = fbgen(1e-5,8000)

nss = 2
MaxIte = 20000

expDem = [.67,1.35]
xi = np.array([[.84,.91],[.63,2.78]],dtype=np.float64) # epsi,delt
Yt = np.array([0,0,1.02],dtype=np.float64) # index: 0,1 [-1]
pt = np.array([0,0,.24],dtype=np.float64)
nt,xt= np.zeros(nss+1),np.zeros(nss+1)
cut = np.zeros((MaxIte,4*B+1))
vldCutNum = 0
for ite in range(MaxIte):
    t=0
    val1,Al0,lb = tmd(0,ite,t,0,0,0)
    prtAct(Al0)
    t=1
    val2,Al1,_ = tmd(0,ite,t,Al0,0,0)
    prtAct(Al1)
    ub = val1 + val2
    print('ub = %8g + %8g = %8g|%8g = lb' % (val1,val2,ub,lb))
    # __________ Start adding cuts __________
    t=1
    pai = tmd(1,ite,t,Al0,1,0) # does it need to add tail function?
    maxerVal,subVal,tailVal,drc,cpv = tmd(1,ite,t,Al0,2,pai)
    cut[ite] = np.append(pai, subVal)
