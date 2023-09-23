import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys

# t = 1: deterministic info, t = 2, t = 3: 4 scenes
# 20220915

def trial(dir, bcwIte, t, epsi, delt, lastAct, cutQtpp): # t = 1,2,...,8
    m = gp.Model("stage_t")
    Y, x, p, n, y, o, D, tha = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    Ylast,plast,nlast,xlast = lastAct
    cn_coeff = 1 if t < nss else 5
    m.setObjective(cp * p + cn_coeff * cn * n + cy * y + co * o + tha)

    m.addConstr(Y == (1 - rho) * epsi + rho * Ylast)
    m.addConstr(n - p - D == nlast - plast - xlast)
    m.addConstr(D - rhoY * expDem[t] * Y == (1 - rhoY) * delt)
    m.addConstr(TB * x - o + TS * y <= Ct)
    m.addConstr(x <= Mt * y)
    m.addConstr(p <= It)
    m.addConstr(x + p <= It)
    m.addConstr(o <= Ot)
    m.addConstr(y <= 1)
    if t < nss:
        for c in range(1,bcwIte+dir): # when eval t=1, this is cut of Q2
            for cuts in range(cutQtpp.shape[1]):
                tmp = cutQtpp[c,cuts,:]
                m.addConstr(tmp[0] * Y + tmp[1] * p + tmp[2] * n + tmp[3] * x - tha <= tmp[-1])
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        sys.exit(3)
    if dir:
        cutCoeff = np.zeros(5, dtype=np.float64)
        l = m.getConstrs() # Y,p,n,x,const
        cutCoeff[0] += rho * l[0].Pi
        cutCoeff[1] += -l[1].Pi
        cutCoeff[2] += l[1].Pi
        cutCoeff[3] += -l[1].Pi
        cutCoeff[4] -= l[0].Pi * (1 - rho) * epsi + l[2].Pi * (1 - rhoY) * delt + l[3].Pi * Ct + (l[5].Pi + l[6].Pi) * It + l[7].Pi * Ot + l[8].Pi
        if t < nss:
            for c in range(0, bcwIte):
                for cuts in range(cutQtpp.shape[1]):
                    cutCoeff[4] -= l[9 + c * sps + cuts].Pi * cutQtpp[c + 1, cuts, -1]
    if dir: # backward
        if t == nss: # last stage eval: return cuts and cost of last stage
            return cutCoeff,m.ObjVal
        else:
            return cutCoeff # backward but not last stage: return only cuts is enough
    else: # forward
        if t == 1: # 1stage cost, action, lb
            return m.ObjVal-tha.X, [Y.X,p.X,n.X,x.X],m.ObjVal
        else:# stage cost, action
            return m.ObjVal-tha.X, [Y.X,p.X,n.X,x.X]


sps = 4
nss = 3
Max_Ite = 15
# cutpo = np.zeros((Max_Ite,nss,5),dtype=np.float64) # global storage
cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
expDem = [0,67,135,200]
epsi2 = [1.26479573,0.56714763,1.47861836,0.88311853]
delt2 = [93.53678368, 75.30641826,9.15385607,162.5951657]
epsi3 = [1.19236862, 0.78056259,1.54051919,1.12249152]
delt3 = [112.56921403, 77.0670528, 52.40610459,127.8916252]
cutQt2 = np.zeros((Max_Ite, 1, 5), dtype=np.float64) # global storage
cutQt3 = np.zeros((Max_Ite, sps, 5), dtype=np.float64) # global storage

val2_A = np.zeros(sps,dtype=np.float64)
rsu2_L = [0 for i in range(sps)]

for bcwIte in range(1,Max_Ite):
    t = 1 # first stage need to do only once
    val1,rsu1,lb = trial(0, bcwIte, t, epsi1, delt1, [Y0, p0, n0, x0], cutQt2)
    t = 2 # second stage:
    for k in range(sps): # trajs --> 4 kinds of actions[2], 4 kinds of stagecost[2]
        val2_A[k],rsu2_L[k] = trial(0, bcwIte, t, epsi2[k], delt2[k], rsu1, cutQt3)
    # -------------------backward ite actually from here-------------------
    t = 3 # end stage evaluation: last stage cost (average) && cutQtend gen
    valk = np.zeros(sps,dtype=np.float64) # cost of traj per traj
    for k in range(sps): # done 4 actions[2], so 4 cuts added to cutQt3
        val3 = 0
        cut = np.zeros(5, dtype=np.float64)
        for o in range(sps): # all possible scenes of end-stage
            ctmp,vtmp = trial(1,bcwIte,3,epsi3[o],delt3[o],rsu2_L[k],0)
            cut += ctmp
            val3 += vtmp
        cutQt3[bcwIte, k, :] = cut / sps
        valk[k] = val1 + val2_A[k] + val3/sps # total cost per traj
    ub = np.average(valk)
    print('%8d | %8g | %8g' % (bcwIte, ub, lb))
    if ub < lb + 5e-4:
        print('Y,p,n,x:',rsu1)
        break
    t = 2
    cut = np.zeros(5,dtype=np.float64)
    for o in range(sps): # all possible scenes of 2-stage
        cut += trial(1, bcwIte, t, epsi2[o], delt2[o], rsu1, cutQt3)
    cutQt2[bcwIte, 0, :] = cut / sps

