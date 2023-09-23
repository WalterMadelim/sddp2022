import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy
# MC-sddp: only one traj per trial
# one problem exists: m\in C(n)
# 20221005
def trial(dir, ite, t, sceThis, epsi, delt, lastAct):
    m = gp.Model("MC-SDDP")
    Y, p, n, x, y, oi, D = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    tha = m.addVars(sps)
    Ylast,plast,nlast,xlast = lastAct
    cn_coeff = 1 if t < nss else 5
    pbd = P[t][:,sceThis] if t > 1 else P[1] # transition prob to the next stage
    m.setObjective(cp * p + cn_coeff * cn * n + cy * y + co * oi + gp.quicksum(pbd[i]*tha[i] for i in range(sps))) # tha represent E[Qt+1]
    m.addConstr(Y == (1 - rho) * epsi + rho * Ylast)
    m.addConstr(n - p - D == nlast - plast - xlast)
    m.addConstr(D - rhoY * expDem[t] * Y == (1 - rhoY) * delt)
    m.addConstr(TB * x - oi + TS * y <= Ct)
    m.addConstr(x <= Mt * y)
    m.addConstr(p <= It)
    m.addConstr(x + p <= It)
    m.addConstr(oi <= Ot)
    m.addConstr(y <= 1)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):
        for its in range(1,ite+dir):
            m.addConstr(tha[o[its,t+1]]>=cutQt[its][t+1,-1]+cutQt[its][t+1,0]*Y+cutQt[its][t+1,1]*p+cutQt[its][t+1,2]*n+cutQt[its][t+1,3]*x)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    if dir: # backward
        cutCoeff = np.zeros(5, dtype=np.float64)
        l = m.getConstrs() # Y,p,n,x,const
        cutCoeff[0] = rho * l[0].Pi
        cutCoeff[1] = -l[1].Pi
        cutCoeff[2] = l[1].Pi
        cutCoeff[3] = -l[1].Pi
        cutCoeff[4] = m.ObjVal - sum(cutCoeff[i]*e for i,e in enumerate(lastAct))
        return cutCoeff
    else: # forward
        return cp * p.X + cn_coeff * cn * n.X + cy * y.X + co * oi.X,[Y.X,p.X,n.X,x.X],m.ObjVal
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
rsu = [0 for i in range(nss+1)]
cutQt = [0 for i in range(Max_Ite)]
for i in range(1,Max_Ite):
    cutQt[i] = np.zeros((nss+1,5))
o = np.zeros((Max_Ite,nss+1),dtype=np.int8)
lbkeep = 20
olb = 0
lbcnt = 0
M = 3000
z_0d05 = 1.64
for ite in range(1, Max_Ite):
    t = 1
    ub, rsu[1], lb = trial(0, ite, t, 0, epsi1, delt1, [Y0, p0, n0, x0])
    if lb - olb < 5e-4:
        lbcnt += 1
    else:
        lbcnt = 0
    olb = lb
    if lbcnt == lbkeep: # do a test
        # keep current ite not change
        t = 1
        tmp, rsu[1], _ = trial(0, ite, t, 0, epsi1, delt1, [Y0, p0, n0, x0])
        trjVal = tmp * np.ones(M)
        for m in range(M):
            t = 2
            o[ite, t] = rpd2o(r(), P[1])
            tmp, rsu[t], _ = trial(0, ite, t, o[ite, t], epsi[o[ite, t], t], delt[o[ite, t], t], rsu[t - 1])
            trjVal[m]+=tmp
            for t in range(3, nss + 1):
                o[ite, t] = rpd2o(r(), P[t - 1][:, o[ite, t - 1]])
                tmp, rsu[t], _ = trial(0, ite, t, o[ite, t], epsi[o[ite, t], t], delt[o[ite, t], t], rsu[t - 1])
                trjVal[m]+=tmp
        muhat, sigmahat = np.average(trjVal), myVar(trjVal) ** 0.5
        ub = muhat + z_0d05 * sigmahat / M ** 0.5  # in statistical meaning
        print('final: statistical ub: %8g | %8g lb' % (ub, lb))
        break
    t = 2
    o[ite,t] = rpd2o(r(), P[1])
    tmp, rsu[t], _ = trial(0, ite, t, o[ite,t], epsi[o[ite,t],t], delt[o[ite,t],t], rsu[t-1])
    ub += tmp
    for t in range(3,nss+1):
        o[ite,t] = rpd2o(r(), P[t-1][:,o[ite,t-1]])
        tmp, rsu[t], _ = trial(0, ite, t, o[ite,t], epsi[o[ite,t],t], delt[o[ite,t],t], rsu[t-1])
        ub += tmp
    print('%8d | %8g | %8g' % (ite, ub, lb))
    # if lb + 5e-5 > ub:
    #     break
    # -------------------backward ite : purely gen cuts -------------------
    # o[1,1] = 0, o[1,2] = 0, o[1,3] = 1, o[1,4] = 1, o[1,5] = 2
    for t in range(nss,1,-1):
        cutQt[ite][t] = trial(1, ite, t, o[ite, t], epsi[o[ite, t], t], delt[o[ite, t], t], rsu[t-1])


