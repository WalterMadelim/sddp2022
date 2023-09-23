import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy

# Markov-Chain-SDDP program for inventory
# 20220917
# doorvanbei
# outcome of MC-SDDP
# 1 | 105821 | 1092.46 | 103126
# 2 | 9677.42 | 4300.7 | 9424.68
# 3 | 10129.5 | 7366.82 | 9528.67
# 4 | 8205.49 | 7913.74 | 7830.82
# 5 | 8360.01 | 7944.47 | 7915.04
# 6 | 9034.57 | 7947.15 | 8379.34
# 7 | 8077.36 | 7947.15 | 7835.71
# 8 | 8812.33 | 7947.15 | 8232.56
# 9 | 9190.37 | 7947.15 | 8439.57
# 10 | 8080.64 | 7947.15 | 7719.55
# 11 | 7784.65 | 7947.15 | 7592.4

def myVar(u): # sample variance
    l = np.size(u)
    return l/(l-1)*np.var(u)
def trial(dir, ite, t, sceThis, epsi, delt, lastAct): # t = 4, epsi = [t=4,o=0]
    m = gp.Model("MC-SDDP")
    Y, x, p, n, y, o, D = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    tha = m.addVars(sps)
    Ylast,plast,nlast,xlast = lastAct
    cn_coeff = 1 if t < nss else 5
    pbd = P[t][:,sceThis] if t > 1 else P[1]
    m.setObjective(cp * p + cn_coeff * cn * n + cy * y + co * o + gp.quicksum( pbd[i]*tha[i] for i in range(sps) )) # tha represent E[Qt+1]
    m.addConstr(Y == (1 - rho) * epsi + rho * Ylast)
    m.addConstr(n - p - D == nlast - plast - xlast)
    m.addConstr(D - rhoY * expDem[t] * Y == (1 - rhoY) * delt)
    m.addConstr(TB * x - o + TS * y <= Ct)
    m.addConstr(x <= Mt * y)
    m.addConstr(p <= It)
    m.addConstr(x + p <= It)
    m.addConstr(o <= Ot)
    m.addConstr(y <= 1)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1): # no tail function for the last stage
        for i in range(sps): # all possible scenes of t+1
            # cutQt[t+1][i][ite, trj, :]
            cpc = cutQt[t+1][i].shape[1]
            for c in range(1, ite + dir): # when eval t=1, this is cut of Q2
                for cuts in range(cpc):
                    tmp = cutQt[t+1][i][c,cuts,:]
                    m.addConstr(tmp[0] * Y + tmp[1] * p + tmp[2] * n + tmp[3] * x - tha[i] <= tmp[-1])
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
        cutCoeff[4] -= l[0].Pi * (1 - rho) * epsi + l[2].Pi * (1 - rhoY) * delt + l[3].Pi * Ct + (l[5].Pi + l[6].Pi) * It + l[7].Pi * Ot + l[8].Pi
        if t == nss: # last stage eval: return cuts and cost of last stage
            return cutCoeff
        else:
            for i in range(sps):
                for c in range(ite):
                    for cuts in range(cpc):
                        pass
                        cutCoeff[4] -= l[9+cuts+c*cpc+i*ite*cpc].Pi * cutQt[t+1][i][c+1,cuts,-1]
            return cutCoeff # backward but not last stage: return only cuts is enough
    else: # forward
        # tmp = y.X
        # tmp1 = x.X
        return cp * p.X + cn_coeff * cn * n.X + cy * y.X + co * o.X,[Y.X,p.X,n.X,x.X],m.ObjVal
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

np.random.seed(3)
z_0d05 = 1.64
cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
sps = 3
Max_Ite = 20
nss = 5
TRJ = 2
expDem = [0,67,135,105,80,79]
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
cutQt = [0 for t in range(nss+1)] # global storage
for t in range(2,nss+1):
    cutQt[t] = [0 for o in range(sps)]
for t in range(2,nss+1):
    cpc = 1 if t == 2 else TRJ
    for o in range(sps):
        cutQt[t][o] = np.zeros((Max_Ite,cpc,5),dtype=np.float64)
rsu = [0 for i in range(nss+1)]
for t in range(2,nss+1):
    rsu[t] = [0 for i in range(TRJ)]
o = sps*np.ones((TRJ,nss+1),dtype=np.int32)
for ite in range(1, Max_Ite):
    t = 1
    tmp, rsu[1], lb = trial(0, ite, t, 0, epsi1, delt1, [Y0, p0, n0, x0])
    trjVal = tmp * np.ones(TRJ,dtype=np.float64)
    for trj in range(TRJ):
        t = 2
        o[trj,t] = rpd2o(r(), P[1])
        tmp,rsu[t][trj],_ = trial(0, ite, t, o[trj,t],epsi[o[trj,t],t], delt[o[trj,t],t], rsu[t-1])
        trjVal[trj] += tmp
        for t in range(3,nss+1):
            o[trj,t] = rpd2o(r(),P[t-1][:,o[trj,t-1]])
            tmp,rsu[t][trj],_ = trial(0, ite, t, o[trj,t],epsi[o[trj,t],t], delt[o[trj,t],t], rsu[t-1][trj])
            trjVal[trj] += tmp
    muhat, sigmahat = np.average(trjVal), myVar(trjVal) ** 0.5
    ub = muhat + z_0d05 * sigmahat / TRJ ** 0.5  # in statistical meaning
    print('%8d | %8g | %8g | %8g' % (ite, ub, lb, muhat))
    if lb + 5e-5 > ub:
        break
    # -------------------backward ite : purely gen cuts -------------------
    for t in range(5,2,-1):
        for trj in range(TRJ):
            for ot in range(sps):
                cutQt[t][ot][ite,trj,:] = trial(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1][trj])
    t = 2
    for ot in range(sps):  # this-stage-all-possible-scenes
        cutQt[t][ot][ite,0,:] = trial(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1])

