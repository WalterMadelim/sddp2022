import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy

# pureIntCut: that is too slow
# doorvanbei
# 20221006

def fwdtrial(dir, ite, t, sceThis, xi, lastAct):
    m = gp.Model("MC-SDDP")
    Ythis,pthis,nthis,xthis = m.addVars(B,ub=1),m.addVars(B,ub=1),m.addVars(B,ub=1),m.addVars(B,ub=1) # copy variable of this stage
    Y, p, n, x = m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY) # chaning to the next stage
    y, oi, D = m.addVar(vtype=GRB.BINARY), m.addVar(ub=Ot), m.addVar()
    tha = m.addVar()
    cn_coeff = 1 if t < nss else 5
    # pbd = P[t][:,sceThis] if t > 1 else P[1] # transition prob to the next stage
    m.setObjective(cp * mb2f(p) + cn_coeff * cn * mb2f(n) + cy * y + co * oi + tha) # tha represent E[Qt+1]
    m.addConstrs(Ythis[i] == lastAct[i] for i in range(B))
    m.addConstrs(pthis[i] == lastAct[B+i] for i in range(B))
    m.addConstrs(nthis[i] == lastAct[2*B+i] for i in range(B))
    m.addConstrs(xthis[i] == lastAct[3*B+i] for i in range(B))
    m.addConstr((1 - rho) * xi[0] + rho * mb2f(Ythis) - mb2f(Y) == [-SCL/2,SCL/2])
    m.addConstr(mb2f(n) - mb2f(p) - D - (mb2f(nthis) - mb2f(pthis) - mb2f(xthis)) == [-SCL/2,SCL/2])
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y) == (1 - rhoY) * xi[1])
    m.addConstr(TB * mb2f(x) - oi + TS * y <= Ct)
    m.addConstr(mb2f(x) <= Mt * y)
    m.addConstr(mb2f(p) <= It)
    m.addConstr(mb2f(x) + mb2f(p) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):# t = 1, ite = 2
        for its in range(1, ite + dir):
            if o[its, t] == sceThis: # always true for t=1
                m.addConstr(tha >= cutQt[its][t + 1, -1] + gp.quicksum(
                    cutQt[its][t + 1, i] * Y[i] for i in range(B)) + gp.quicksum(
                    cutQt[its][t + 1, B + i] * p[i] for i in range(B)) + gp.quicksum(
                    cutQt[its][t + 1, 2 * B + i] * n[i] for i in range(B)) + gp.quicksum(
                    cutQt[its][t + 1, 3 * B + i] * x[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    actionThis = np.zeros(4*B,dtype=np.int8)
    for i in range(B):
        actionThis[i],actionThis[B+i],actionThis[2*B+i],actionThis[3*B+i] = Y[i].X,p[i].X,n[i].X,x[i].X
    return m.ObjVal-tha.X,actionThis,m.ObjVal
def icutval(dir, ite, t, sceThis, xi, lastAct):
    m = gp.Model("MC-SDDP")
    Ythis, pthis, nthis, xthis = m.addVars(B, ub=1), m.addVars(B, ub=1), m.addVars(B, ub=1), m.addVars(B,ub=1)  # copy variable of this stage
    Y, p, n, x = m.addVars(B, vtype=GRB.BINARY), m.addVars(B, vtype=GRB.BINARY), m.addVars(B,vtype=GRB.BINARY), m.addVars(B, vtype=GRB.BINARY)  # chaning to the next stage
    y, oi, D = m.addVar(vtype=GRB.BINARY), m.addVar(ub=Ot), m.addVar()
    tha = m.addVar()
    cn_coeff = 1 if t < nss else 5
    # pbd = P[t][:, sceThis] if t > 1 else P[1]  # transition prob to the next stage
    m.setObjective(cp * mb2f(p) + cn_coeff * cn * mb2f(n) + cy * y + co * oi + tha)  # tha represent E[Qt+1]
    m.addConstrs(Ythis[i] == lastAct[i] for i in range(B))
    m.addConstrs(pthis[i] == lastAct[B + i] for i in range(B))
    m.addConstrs(nthis[i] == lastAct[2 * B + i] for i in range(B))
    m.addConstrs(xthis[i] == lastAct[3 * B + i] for i in range(B))
    m.addConstr((1 - rho) * xi[0] + rho * mb2f(Ythis) - mb2f(Y) == [-SCL / 2, SCL / 2])
    m.addConstr(mb2f(n) - mb2f(p) - D - (mb2f(nthis) - mb2f(pthis) - mb2f(xthis)) == [-SCL / 2, SCL / 2])
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y) == (1 - rhoY) * xi[1])
    m.addConstr(TB * mb2f(x) - oi + TS * y <= Ct)
    m.addConstr(mb2f(x) <= Mt * y)
    m.addConstr(mb2f(p) <= It)
    m.addConstr(mb2f(x) + mb2f(p) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):
        for its in range(1,ite+dir): # need to be fixed
            if o[its,t] == sceThis:
                m.addConstr(tha >= cutQt[its][t+1,-1] + gp.quicksum(cutQt[its][t+1,i]*Y[i] for i in range(B)) + gp.quicksum(cutQt[its][t+1,B+i]*p[i] for i in range(B))+ gp.quicksum(cutQt[its][t+1,2*B+i]*n[i] for i in range(B)) + gp.quicksum(cutQt[its][t+1,3*B+i]*x[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    return m.ObjVal
def intcut(dir, ite, t, sceLast, lastAct):
    vN = 0
    pbd = P[t-1][:, sceLast] if t > 2 else P[1]  # transition prob to the next stage
    for no in range(sps):
        vN += pbd[no]*icutval(dir, ite, t, no, [epsi[no, t], delt[no, t]], lastAct)
    return np.concatenate((vN*(2*lastAct-1),np.array([(1-sum(lastAct))*vN])))
def r():
    return np.random.random()
def rpd2o(r,pd):
    t = copy.deepcopy(pd)
    for i in range(1,len(t)):
        t[i] += t[i-1]
    s = -1
    for e in t:
        s += 1
        if r < e:
            break
    return s
def myVar(u): # sample variance
    l = np.size(u)
    return l/(l-1)*np.var(u)
def fbgen(n): # float basis generation
    fb = []
    for i in range(n):
        fb.append(2 ** i)
    return np.array(fb) / 2 ** (n//2+1)
def mb2f(b): # b2f in gurobi model
    return gp.quicksum(fb[i]*b[i] for i in range(B))
def b2f(b): # binary vector to a float number
    return sum(fb[i]*b[i] for i in range(B))
def f2bL(f): # float to binary List
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
        sys.exit(3)
    xthis = [0 for i in range(B)]
    for i in range(B):
        xthis[i] = x[i].X
    return xthis
def xb2f(b): # binary vector to a float number
    return sum(fb[i]*b[i].X for i in range(B))
B = 27 # use a binary vector (len=27) to approx a float
fb = fbgen(B)
SCL = fb[0]
np.random.seed(3)
cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
sps = 3
Max_Ite = 200
nss = 5 # 1,2,3,4,5
expDem = [0,67,135,105,80,79] # [0]: null, [1-5]: the nss stages
epsi = np.array([[0.,0.,0.91491305,0.36565183,1.19035292,1.19236862],
                 [0.,0.,1.26479573,0.56714763,1.47861836,0.88311853],
                 [0.,0.,2.3886352 ,0.60174679,1.38746341,0.99567512]])
delt = np.array([[0.,0.,277.785241  , 74.05659523,107.67428731,112.56921403],
                 [0.,0.,133.9353051 ,221.57483785, 48.64062247,170.8259371 ],
                 [0.,0.,156.53066711, 65.09750764, 65.67303495, 48.7458782 ]])
# specify Markov Chain transition matrix series P(list)
P = [0,0,2,3,4,0] # P[2] is P(2->3), P[5] is fictitious
P[1] = np.array([.1,.8,.1]) # for the initial single point
P[2],P[3],P[4] = np.array([[.8,.1,.1],[.1,.8,.1],[.1,.1,.8]]),np.array([[.8,.1,.1],[.1,.8,.1],[.1,.1,.8]]),np.array([[.8,.1,.1],[.1,.8,.1],[.1,.1,.8]])
P[5] = np.array([[.0,.0,.0],[.0,.0,.0],[.0,.0,.0]]) # fictitious
# define 3 storages: cutQt, rsu, o
rsu = [0 for i in range(nss+1)] # rsu[t] is a 1-d np-array
cutQt = [0 for i in range(Max_Ite)]
for i in range(1,Max_Ite):
    cutQt[i] = np.zeros((nss+1,4*B+1))
o = np.zeros((Max_Ite,nss+1),dtype=np.int8)
lbkeep = 20
olb = 0
lbcnt = 0
M = 300
z_0d05 = 1.64
Y0,p0,n0,x0 = f2bL(Y0),f2bL(p0),f2bL(n0),f2bL(x0) # for SDDiP
for ite in range(1, Max_Ite):
    t = 1
    ub, rsu[1], lb = fwdtrial(0, ite, t, 0, [epsi1, delt1], Y0+p0+n0+x0)
    t = 2
    o[ite,t] = rpd2o(r(), P[1])
    tmp, rsu[t], _ = fwdtrial(0, ite, t, o[ite,t], [epsi[o[ite,t],t], delt[o[ite,t],t]], rsu[t-1])
    ub += tmp
    for t in range(3,nss+1):
        o[ite,t] = rpd2o(r(), P[t-1][:,o[ite,t-1]])
        tmp, rsu[t], _ = fwdtrial(0, ite, t, o[ite,t], [epsi[o[ite, t], t], delt[o[ite, t], t]], rsu[t-1])
        ub += tmp
    print('%8d | %8g | %8g' % (ite, ub, lb))
    # if lb + 5e-5 > ub:
    #     break
    # -------------------backward ite : purely gen cuts -------------------
    # o[1,1] = 0, o[1,2] = 0, o[1,3] = 1, o[1,4] = 1, o[1,5] = 2
    for t in range(nss,1,-1): # t = 4, do Q5
        cutQt[ite][t] = intcut(1, ite, t, o[ite, t-1], rsu[t-1])



# convergence of benders' for LP
      #  1 |   105590 |  1092.46
      #  2 |  8971.86 |  3892.47
      #  3 |  18703.2 |  5416.72
      #  4 |  8450.34 |  5806.44
      #  5 |  8719.59 |  6110.14
      #  6 |  7445.44 |  6261.31
      #  7 |  7522.54 |  6630.46
      #  8 |  8934.12 |  6711.84
      #  9 |  15893.1 |  6830.31
      # 10 |  7414.31 |  7537.12
      # 11 |   7390.9 |  7573.23
      # 12 |  7385.23 |  7694.97
      # 13 |  7283.68 |  7712.31
      # 14 |  7449.36 |  7722.21
      # 15 |  7285.62 |  7724.97
      # 16 |  7281.48 |  7739.85
      # 17 |   7641.2 |  7740.38
      # 18 |  10852.4 |  7741.92
      # 19 |  7452.37 |   7800.5
      # 20 |  8382.29 |  7813.18
      # 21 |  7978.87 |  7879.76
      # 22 |  7278.73 |  7943.14
      # 23 |  10546.5 |  7943.14
      # 24 |  7947.18 |  7943.41
      # 25 |  7278.73 |  7943.41
      # 26 |  7947.18 |  7943.41
      # 27 |  7640.43 |  7943.41
      # 28 |  10538.7 |  7943.41
      # 29 |  9376.82 |  7943.41
      # 30 |  7452.37 |  7943.41
      # 31 |  7278.73 |  7943.41
      # 32 |  7765.47 |  7943.41
      # 33 |  8328.83 |  7943.41
      # 34 |  7278.73 |  7945.07
      # 35 |  7453.14 |  7945.38
      # 36 |  7278.73 |  7946.92
      # 37 |  7278.73 |  7946.92
      # 38 |  7278.73 |  7946.92
      # 39 |  7278.73 |  7946.92
      # 40 |  10538.7 |  7946.92
      # 41 |  7852.45 |  7947.15
      # 42 |  7452.37 |  7947.15
      # 43 |  7278.73 |  7947.15
      # 44 |  7452.37 |  7947.15
      # 45 |  7278.73 |  7947.15
      # 46 |  10538.7 |  7947.15
      # 47 |  7640.43 |  7947.15
      # 48 |  7278.73 |  7947.15
      # 49 |  10043.9 |  7947.15
      # 50 |  7652.14 |  7947.15
      # 51 |  7826.31 |  7947.15
      # 52 |  7278.73 |  7947.15
      # 53 |  7278.73 |  7947.15
      # 54 |  7278.73 |  7947.15
      # 55 |  7640.43 |  7947.15
      # 56 |  7278.73 |  7947.15
      # 57 |  7452.37 |  7947.15
      # 58 |  7278.73 |  7947.15
      # 59 |  7278.73 |  7947.15
      # 60 |  10201.1 |  7947.15
# final: statistical ub:  7983.52 |  7947.15 lb # M = 3000