import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy
# pre-SDDP: deterministic extensive form: 1-4-4-4_4stage problem
# doorvanbei
# 20220915

def trial(dir, ite, t, epsi, delt, lastAct, cutQtpp): # t = 1,2,...,8
    m = gp.Model("0f1b")
    Y, x, p, n, y, o, D, tha = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    Ylast,plast,nlast,xlast = lastAct
    cn_coeff = 1 if t < nss else 5
    m.setObjective(cp * p + cn_coeff * cn * n + cy * y + co * o + tha) # tha represent E[Qt+1]

    m.addConstr(Y == (1 - rho) * epsi + rho * Ylast)
    m.addConstr(n - p - D == nlast - plast - xlast)

    m.addConstr(D - rhoY * expDem[t] * Y == (1 - rhoY) * delt)
    m.addConstr(TB * x - o + TS * y <= Ct)
    m.addConstr(x <= Mt * y)
    m.addConstr(p <= It)
    m.addConstr(x + p <= It)
    m.addConstr(o <= Ot)
    m.addConstr(y <= 1)
    if t < nss: # no tail function for the last stage
        cpc = cutQtpp.shape[1]
        for c in range(1, ite + dir): # when eval t=1, this is cut of Q2
            for cuts in range(cpc):
                tmp = cutQtpp[c,cuts,:]
                m.addConstr(tmp[0] * Y + tmp[1] * p + tmp[2] * n + tmp[3] * x - tha <= tmp[-1])
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        sys.exit(3)
    if dir: # backward
        cutCoeff = np.zeros(5, dtype=np.float64)
        l = m.getConstrs() # Y,p,n,x,const
        cutCoeff[0] += rho * l[0].Pi
        cutCoeff[1] -= l[1].Pi
        cutCoeff[2] += l[1].Pi
        cutCoeff[3] -= l[1].Pi
        cutCoeff[4] = -(m.ObjVal - cutCoeff[0] * Ylast - cutCoeff[1] * plast - cutCoeff[2] * nlast - cutCoeff[3] * xlast)
        # cutCoeff[4] -= l[0].Pi * (1 - rho) * epsi + l[2].Pi * (1 - rhoY) * delt + l[3].Pi * Ct + (l[5].Pi + l[6].Pi) * It + l[7].Pi * Ot + l[8].Pi
        if t == nss: # last stage eval: return cuts and cost of last stage
            return cutCoeff,m.ObjVal
        else:
        #     for c in range(0, ite):
        #         for cuts in range(cpc):
        #             cutCoeff[4] -= l[9 + c * cpc + cuts].Pi * cutQtpp[c + 1, cuts, -1]
            return cutCoeff # backward but not last stage: return only cuts is enough
    else: # forward
        if t == 1: # 1stage cost, action, lb
            return m.ObjVal-tha.X,[Y.X,p.X,n.X,x.X],m.ObjVal
        else:# stage cost, action
            return m.ObjVal-tha.X,[Y.X,p.X,n.X,x.X],m.ObjVal

sps = 4
nss = 4
Max_Ite = 30
cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
expDem = [0,67,135,100,84]
delt2 = [93.53678368, 75.30641826,9.15385607,162.5951657]
delt3 = [112.56921403, 77.0670528, 52.40610459,127.8916252]
delt4 = [156.53066711, 65.09750764, 65.67303495, 48.7458782]
epsi2 = [1.26479573,0.56714763,1.47861836,0.88311853]
epsi3 = [1.19236862, 0.78056259,1.54051919,1.12249152]
epsi4 = [2.3886352 ,0.60174679,1.38746341,0.99567512]

cutQt2 = np.zeros((Max_Ite,1,5), dtype=np.float64) # global storage
cutQt3 = np.zeros((Max_Ite,sps,5), dtype=np.float64) # global storage
cutQt4 = np.zeros((Max_Ite,sps**2,5), dtype=np.float64) # global storage

rsu2ps2 = [0 for i in range(sps)]
rsu3ps2s3 = [copy.deepcopy(rsu2ps2) for i in range(sps)]

val2ps2 = np.zeros(sps,dtype=np.float64) # val of stage2 per s2
val3ps2s3 = np.zeros((sps,sps),dtype=np.float64) # val of stage2 per s2
val4ps2s3o = np.zeros((sps,sps,sps),dtype=np.float64) # val of

lb234sub = np.zeros(sps,dtype = np.float64) # test Q3 -> is not correct
lb34sub = np.zeros((sps,sps),dtype=np.float64) # test Q4 -> proved to be accurate

for ite in range(1, Max_Ite):
    t = 1 # first stage need to do only once
    val1,rsu1,lb = trial(0, ite, t, epsi1, delt1, [Y0, p0, n0, x0], cutQt2)
    # print('c1 %8g' % (val1))
    for s2 in range(sps): # actually this is enum(trajs)
        t = 2
        val2ps2[s2],rsu2ps2[s2],lb234sub[s2] = trial(0, ite, t, epsi2[s2], delt2[s2], rsu1, cutQt3) # f(x1,k2[s2]) + Q3()
        for s3 in range(sps):
            t = 3
            val3ps2s3[s2,s3], rsu3ps2s3[s2][s3],lb34sub[s2,s3] = trial(0, ite, t, epsi3[s3], delt3[s3], rsu2ps2[s2], cutQt4)
    # -------------------backward ite actually from here-------------------
    t = 4 # end stage evaluation: last stage cost (average) && cutQtend gen
    for s2 in range(sps):
        for s3 in range(sps): # dual fors: number of action3
            cut = np.zeros(5, dtype=np.float64)
            for o in range(sps): # all possible scenes of end-stage (4)
                ctmp,val4ps2s3o[s2,s3,o] = trial(1, ite, t, epsi4[o], delt4[o],rsu3ps2s3[s2][s3], 0)
                cut += ctmp
            cutQt4[ite, sps * s2 + s3, :] = cut/sps
            # print('[%d]: %g || %g ' % ( s2*sps+s3,val3ps2s3[s2,s3]+np.average(val4ps2s3o[s2,s3,:]),lb34sub[s2,s3]))
        # print('[%d]: %g || %g ' % ( s2, val2ps2[s2] + np.average(val3ps2s3[s2,:]) + np.average(val4ps2s3o[s2,:,:]) ,lb234sub[s2]))
    ub = val1 + np.average(val2ps2) + np.average(val3ps2s3) + np.average(val4ps2s3o)
    print('%8d | %8g | %8g' % (ite, ub, lb))
    if lb + 5e-5 > ub:
        break
    t = 3
    for s2 in range(sps): # number of action[2]
        cut = np.zeros(5,dtype=np.float64)
        for o in range(sps): # all possible scenes of stage t
            cut += trial(1, ite, t, epsi3[o], delt3[o], rsu2ps2[s2], cutQt4)
        cutQt3[ite, s2, :] = cut / sps
    t = 2
    cut = np.zeros(5,dtype=np.float64)
    for o in range(sps): # all possible scenes of stage t
        cut += trial(1, ite, t, epsi2[o], delt2[o], rsu1, cutQt3)
    cutQt2[ite, 0, :] = cut / sps


