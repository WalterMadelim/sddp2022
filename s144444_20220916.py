import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy
# pre-SDDP: deterministic extensive form: 1-4-4-4-4-4 stage problem
# results: about 10 min
# Academic license - for non-commercial use only - expires 2023-07-26
#        1 |   125089 |  1092.46
#        2 |  11169.4 |  4154.88
#        3 |  14992.4 |  7804.58
#        4 |  10314.7 |  8590.45
#        5 |    11286 |  9319.51
#        6 |    10304 |  9741.34
#        7 |  10434.2 |  10026.3
#        8 |  10226.5 |  10086.4
#        9 |  10215.3 |  10130.1
#       10 |  10186.9 |  10186.9
#
# Process finished with exit code 0
# doorvanbei
# 20220916
def getUb(val1,val):
    v = val1
    for i in range(2,nss+1):
        v += np.average(val[i])
    return v
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
        cutCoeff[4] -= l[0].Pi * (1 - rho) * epsi + l[2].Pi * (1 - rhoY) * delt + l[3].Pi * Ct + (l[5].Pi + l[6].Pi) * It + l[7].Pi * Ot + l[8].Pi
        if t == nss: # last stage eval: return cuts and cost of last stage
            return cutCoeff,m.ObjVal
        else:
            for c in range(0, ite):
                for cuts in range(cpc):
                    cutCoeff[4] -= l[9 + c * cpc + cuts].Pi * cutQtpp[c + 1, cuts, -1]
            return cutCoeff # backward but not last stage: return only cuts is enough
    else: # forward
        if t == 1: # 1stage cost, action, lb
            return m.ObjVal-tha.X,[Y.X,p.X,n.X,x.X],m.ObjVal
        else:# stage cost, action
            return m.ObjVal-tha.X,[Y.X,p.X,n.X,x.X],m.ObjVal

cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
sps = 4

Max_Ite = 30
nss = 6
expDem = [0,67,135,105,80,79,72]
epsi = np.array([[0.,0.,0.91491305,0.36565183,1.19035292,1.19236862,0.78056259],
                 [0.,0.,1.26479573,0.56714763,1.47861836,0.88311853,0.66679045],
                 [0.,0.,1.16396591,1.15079421,0.91888046,0.4409381 ,0.72007593],
                 [0.,0.,2.3886352 ,0.60174679,1.38746341,0.99567512,0.51351849]])
delt = np.array([[0.,0.,277.785241  , 74.05659523,107.67428731,112.56921403, 77.0670528 ],
                 [0.,0.,182.78231243, 93.53678368, 75.30641826,  9.15385607,162.59516571],
                 [0.,0.,133.9353051 ,221.57483785, 48.64062247,170.8259371 , 21.2438368 ],
                 [0.,0.,156.53066711, 65.09750764, 65.67303495, 48.7458782 , 91.16585795]])

cutQt = [0 for i in range(nss+1)]
for t in range(2,nss+1):
    cutQt[t] = np.zeros((Max_Ite, sps ** (t - 2), 5), dtype=np.float64) # global storage
rsu = [0 for i in range(nss)] # (0,1,) 2,3,
rsu[2] = [0 for i in range(sps)] # simplex layer: 4 actions
for t in range(3,nss):
    rsu[t] = [copy.deepcopy(rsu[t-1]) for i in range(sps)]
val = [0 for i in range(nss+1)] # (0,1,) 2,3,4,
for t in range(2,nss+1):
    val[t] = np.zeros([sps for i in range(t-1)],dtype=np.float64)

# for debugging use
# ckQ6 = copy.deepcopy(val[5]) # tight
# ckQ5 = copy.deepcopy(val[4])
for ite in range(1, Max_Ite):
    t = 1 # first stage need to do only once
    val1,rsu1,lb = trial(0, ite, t, epsi1, delt1, [Y0, p0, n0, x0], cutQt[t+1])
    s = np.arange(nss) # per ite
    for s[2] in range(sps): # actually this is enum(trajs)
        t = 2
        val[t][s[2]],rsu[t][s[2]],_ = trial(0, ite, t, epsi[s[2],t], delt[s[2],t], rsu1, cutQt[t+1]) # f(x1,k2[s2]) + Q3()
        for s[3] in range(sps):
            t = 3
            val[t][s[2],s[3]], rsu[t][s[2]][s[3]],_ = trial(0, ite, t, epsi[s[3],t], delt[s[3],t], rsu[t-1][s[2]], cutQt[t+1])
            for s[4] in range(sps):
                t = 4
                if ite == 2:
                    pass
                val[t][s[2],s[3],s[4]], rsu[t][s[2]][s[3]][s[4]], _ = trial(0, ite, t, epsi[s[4], t], delt[s[4], t],rsu[t-1][s[2]][s[3]], cutQt[t + 1])
                for s[5] in range(sps):
                    t = 5
                    val[t][s[2], s[3], s[4],s[5]], rsu[t][s[2]][s[3]][s[4]][s[5]], _ = trial(0, ite, t, epsi[s[5], t],
                                                                                  delt[s[5], t], rsu[t-1][s[2]][s[3]][s[4]],
                                                                                  cutQt[t + 1])

    # -------------------backward ite actually from here-------------------
    t = 6
    for s[2] in range(sps):
        for s[3] in range(sps):  # dual fors: number of action3
            for s[4] in range(sps):
                for s[5] in range(sps):
                    cut = np.zeros(5, dtype=np.float64)
                    for o in range(sps):
                        ctmp, val[t][s[2], s[3],s[4],s[5],o] = trial(1, ite, t, epsi[o, t], delt[o, t], rsu[t - 1][s[2]][s[3]][s[4]][s[5]],0)
                        cut += ctmp
                    # print('[%d]: %g || %g ' % ( s[5]+s[4]*sps+s[3]*sps**2+s[2]*sps**3,val[5][s[2], s[3],s[4],s[5]]+np.average(val[6][s[2], s[3],s[4],s[5],:]),ckQ6[s[2], s[3], s[4],s[5]]))
                    cutQt[t][ite,s[5]+s[4]*sps+s[3]*sps**2+s[2]*sps**3,:] = cut/sps

    ub = getUb(val1,val)
    print('%8d | %8g | %8g' % (ite, ub, lb))
    if lb + 5e-5 > ub:
        break

    t = 5
    for s[2] in range(sps):
        for s[3] in range(sps):  # dual fors: number of action3
            for s[4] in range(sps):
                cut = np.zeros(5, dtype=np.float64)
                for o in range(sps):
                    cut+=trial(1,ite,t,epsi[o,t],delt[o,t],rsu[t-1][s[2]][s[3]][s[4]],cutQt[t+1])
                # print('[%d]: %g || %g ' % ( s[4]+s[3]*sps+s[2]*sps**2,val[4][s[2],s[3],s[4]]+np.average(val[5][s[2],s[3],s[4],:])+np.average(val[6][s[2],s[3],s[4],:,:]),ckQ5[s[2], s[3], s[4]]))
                cutQt[t][ite,s[4]+s[3]*sps+s[2]*sps**2,:] = cut/sps

    t = 4
    for s[2] in range(sps):
        for s[3] in range(sps):  # dual fors: number of action3
            cut = np.zeros(5,dtype=np.float64)
            for o in range(sps):
                cut+=trial(1,ite,t,epsi[o,t],delt[o,t],rsu[t-1][s[2]][s[3]],cutQt[t+1])
            cutQt[t][ite,s[3]+s[2]*sps,:] = cut/sps

    t = 3
    for s[2] in range(sps): # number of action[2]
        cut = np.zeros(5,dtype=np.float64)
        for o in range(sps): # all possible scenes of stage t
            cut += trial(1, ite, t, epsi[o,t], delt[o,t], rsu[t-1][s[2]], cutQt[t+1])
        cutQt[t][ite, s[2], :] = cut / sps

    t = 2
    cut = np.zeros(5,dtype=np.float64)
    for o in range(sps): # all possible scenes of stage t
        cut += trial(1, ite, t, epsi[o,t], delt[o,t], rsu1, cutQt[t+1])
    cutQt[t][ite, 0, :] = cut / sps


