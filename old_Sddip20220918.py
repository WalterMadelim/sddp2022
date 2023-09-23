import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
import copy
# Markov-Chain-SDDiP using sBcut and intCut
# results: it seems that the intCut is of no use
# You can just remove the intCut and use only sBcut
# 20220917
# doorvanbei
def trial(dir, ite, t, sceThis, epsi, delt, lastAct): # t = 4, epsi = [t=4,o=0]
    Y_bil,p_bil,n_bil,x_bil = lastAct
    m = gp.Model("MC-SDDiP")
    # 4 chain variables, all continuous, but all represented by bin vector.
    Y_bi, x_bi, p_bi, n_bi = m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY),m.addVars(B,vtype=GRB.BINARY)
    y, o, D = m.addVar(vtype=GRB.BINARY), m.addVar(ub=Ot), m.addVar()
    tha = m.addVars(sps) # all possible scenes at next stage
    cn_coeff = 1 if t < nss else 5
    pbd = P[t][:,sceThis] if t > 1 else P[1]
    m.setObjective(cp * mb2f(p_bi) + cn_coeff * cn * mb2f(n_bi) + cy * y + co * o + gp.quicksum(pbd[i]*tha[i] for i in range(sps)))
    # linking constraints
    m.addConstr(mb2f(Y_bi) <= (1 - rho) * epsi + rho * b2f(Y_bil) + SCL/2)
    m.addConstr(mb2f(Y_bi) >= (1 - rho) * epsi + rho * b2f(Y_bil) - SCL/2)
    m.addConstr(mb2f(n_bi) - mb2f(p_bi) - D == b2f(n_bil) - b2f(p_bil) - b2f(x_bil))
    # below are local constraints
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y_bi) - (1 - rhoY) * delt == [-SCL/2,SCL/2])
    m.addConstr(TB * mb2f(x_bi) - o + TS * y <= Ct)
    m.addConstr(mb2f(x_bi) <= Mt * y)
    m.addConstr(mb2f(p_bi) <= It)
    m.addConstr(mb2f(x_bi) + mb2f(p_bi) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):  # no tail function for the last stage
        for i in range(sps):
            for c in range(1, ite + dir):  # when eval t=1, this is cut of Q2
                cst = cutQt[t + 1][i][c][0][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][0][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
                cst = cutQt[t + 1][i][c][1][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][1][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    ythis,pthis,nthis,xthis = [0 for i in range(B)],[0 for i in range(B)],[0 for i in range(B)],[0 for i in range(B)]
    for i in range(B):
        ythis[i],pthis[i],nthis[i],xthis[i] = Y_bi[i].X,p_bi[i].X,n_bi[i].X,x_bi[i].X
    return cp * xb2f(p_bi) + cn_coeff * cn * xb2f(n_bi) + cy * y.X + co * o.X,[ythis,pthis,nthis,xthis],m.ObjVal
def bw1stLP(dir, ite, t, sceThis, epsi, delt, lastAct):
    pbd = P[t][:, sceThis] if t > 1 else P[1] # used to indicate the tail Qt+1
    Y_bil,p_bil,n_bil,x_bil = lastAct
    stt = GRB.CONTINUOUS
    m = gp.Model("bwd-1st-LP")
    zY_bil, zx_bil, zp_bil, zn_bil = m.addVars(B),m.addVars(B),m.addVars(B),m.addVars(B)
    Y_bi, x_bi, p_bi, n_bi = m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt)
    y, o, D = m.addVar(ub=1,vtype=stt), m.addVar(ub=Ot), m.addVar()
    tha = m.addVars(sps)  # all possible scenes at next stage
    cn_coeff = 1 if t < nss else 5
    m.setObjective(cp * mb2f(p_bi) + cn_coeff * cn * mb2f(n_bi) + cy * y + co * o + gp.quicksum(pbd[i] * tha[i] for i in range(sps)))
    # linking constraints
    m.addConstrs(zY_bil[i] == Y_bil[i] for i in range(B))
    m.addConstrs(zp_bil[i] == p_bil[i] for i in range(B))
    m.addConstrs(zn_bil[i] == n_bil[i] for i in range(B))
    m.addConstrs(zx_bil[i] == x_bil[i] for i in range(B))
    # below are local constraints
    m.addConstr(mb2f(Y_bi) <= (1 - rho) * epsi + rho * b2f(zY_bil) + SCL / 2)
    m.addConstr(mb2f(Y_bi) >= (1 - rho) * epsi + rho * b2f(zY_bil) - SCL / 2)
    m.addConstr(mb2f(n_bi) - mb2f(p_bi) - D == b2f(zn_bil) - b2f(zp_bil) - b2f(zx_bil))
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y_bi) - (1 - rhoY) * delt == [-SCL / 2, SCL / 2])
    m.addConstr(TB * mb2f(x_bi) - o + TS * y <= Ct)
    m.addConstr(mb2f(x_bi) <= Mt * y)
    m.addConstr(mb2f(p_bi) <= It)
    m.addConstr(mb2f(x_bi) + mb2f(p_bi) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):  # no tail function for the last stage
        for i in range(sps):
            for c in range(1, ite + dir):  # when eval t=1, this is cut of Q2
                cst = cutQt[t + 1][i][c][0][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][0][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
                cst = cutQt[t + 1][i][c][1][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][1][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    l = m.getConstrs()  # Y,p,n,x,const
    paiY, paix, paip, pain = [0 for i in range(B)], [0 for i in range(B)], [0 for i in range(B)], [0 for i in range(B)]
    for i in range(B):
        paiY[i] = l[i].Pi
    for i in range(B):
        paip[i] = l[B + i].Pi
    for i in range(B):
        pain[i] = l[2 * B + i].Pi
    for i in range(B):
        paix[i] = l[3 * B + i].Pi
    return [paiY, paip, pain, paix]
def bwTrial(dir, ite, t, sceThis, epsi, delt, lastAct):
    paiY, paip, pain, paix = bw1stLP(dir,ite,t,sceThis,epsi,delt,lastAct) # 1st calculate
    # 2nd calculate
    pbd = P[t][:, sceThis] if t > 1 else P[1] # used to indicate the tail Qt+1
    stt = GRB.BINARY # an integer programming
    m = gp.Model("bwd-2nd-iP")
    zY_bil, zx_bil, zp_bil, zn_bil = m.addVars(B,ub=1),m.addVars(B,ub=1),m.addVars(B,ub=1),m.addVars(B,ub=1)
    Y_bi, x_bi, p_bi, n_bi = m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt)
    y, o, D = m.addVar(ub=1,vtype=stt), m.addVar(ub=Ot), m.addVar()
    tha = m.addVars(sps)  # all possible scenes at next stage
    cn_coeff = 1 if t < nss else 5
    m.setObjective(cp * mb2f(p_bi) + cn_coeff * cn * mb2f(n_bi) + cy * y + co * o + gp.quicksum(pbd[i] * tha[i] for i in range(sps)) - gp.quicksum(paiY[i]*zY_bil[i] for i in range(B)) - gp.quicksum(pain[i]*zn_bil[i] for i in range(B)) - gp.quicksum(paip[i]*zp_bil[i] for i in range(B)) - gp.quicksum(paix[i]*zx_bil[i] for i in range(B)) )
    m.addConstr(mb2f(Y_bi) <= (1 - rho) * epsi + rho * b2f(zY_bil) + SCL / 2)
    m.addConstr(mb2f(Y_bi) >= (1 - rho) * epsi + rho * b2f(zY_bil) - SCL / 2)
    m.addConstr(mb2f(n_bi) - mb2f(p_bi) - D == b2f(zn_bil) - b2f(zp_bil) - b2f(zx_bil))
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y_bi) - (1 - rhoY) * delt == [-SCL / 2, SCL / 2])
    m.addConstr(TB * mb2f(x_bi) - o + TS * y <= Ct)
    m.addConstr(mb2f(x_bi) <= Mt * y)
    m.addConstr(mb2f(p_bi) <= It)
    m.addConstr(mb2f(x_bi) + mb2f(p_bi) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):  # no tail function for the last stage
        for i in range(sps):
            for c in range(1, ite + dir):  # when eval t=1, this is cut of Q2
                cst = cutQt[t + 1][i][c][0][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][0][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
                cst = cutQt[t + 1][i][c][1][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][1][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    return m.ObjVal,[paiY, paip, pain, paix]
def bwintCut(dir, ite, t, sceThis, epsi, delt, lastAct):
    Y_bil,p_bil,n_bil,x_bil = lastAct
    pbd = P[t][:, sceThis] if t > 1 else P[1] # used to indicate the tail Qt+1
    stt = GRB.BINARY # an integer programming
    m = gp.Model("bwd-intCut")
    Y_bi, x_bi, p_bi, n_bi = m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt),m.addVars(B,ub=1,vtype=stt)
    y, o, D = m.addVar(ub=1,vtype=stt), m.addVar(ub=Ot), m.addVar()
    tha = m.addVars(sps)  # all possible scenes at next stage
    cn_coeff = 1 if t < nss else 5
    m.setObjective(cp * mb2f(p_bi) + cn_coeff * cn * mb2f(n_bi) + cy * y + co * o + gp.quicksum(pbd[i] * tha[i] for i in range(sps)))
    m.addConstr(mb2f(Y_bi) <= (1 - rho) * epsi + rho * b2f(Y_bil) + SCL / 2)
    m.addConstr(mb2f(Y_bi) >= (1 - rho) * epsi + rho * b2f(Y_bil) - SCL / 2)
    m.addConstr(mb2f(n_bi) - mb2f(p_bi) - D == b2f(n_bil) - b2f(p_bil) - b2f(x_bil))
    m.addConstr(D - rhoY * expDem[t] * mb2f(Y_bi) - (1 - rhoY) * delt == [-SCL / 2, SCL / 2])
    m.addConstr(TB * mb2f(x_bi) - o + TS * y <= Ct)
    m.addConstr(mb2f(x_bi) <= Mt * y)
    m.addConstr(mb2f(p_bi) <= It)
    m.addConstr(mb2f(x_bi) + mb2f(p_bi) <= It)
    if (t < nss and ite > 1) or (t < nss and dir and ite == 1):  # no tail function for the last stage
        for i in range(sps):
            for c in range(1, ite + dir):  # when eval t=1, this is cut of Q2
                cst = cutQt[t + 1][i][c][0][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][0][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
                cst = cutQt[t + 1][i][c][1][0]
                c_Y, c_p, c_n, c_x = cutQt[t + 1][i][c][1][1]
                m.addConstr(tha[i] >= cst + gp.quicksum(c_Y[i] * Y_bi[i] for i in range(B)) + gp.quicksum(
                    c_p[i] * p_bi[i] for i in range(B)) + gp.quicksum(c_x[i] * x_bi[i] for i in range(B)) + gp.quicksum(
                    c_n[i] * n_bi[i] for i in range(B)))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('>>>>>>>>>>>>> opt Fail >>>>>>>>>>>>>>>>>>>>>>>>>>')
        sys.exit(3)
    paiY, paix, paip, pain = [0 for i in range(B)], [0 for i in range(B)], [0 for i in range(B)], [0 for i in range(B)]
    for i in range(B):
        paiY[i] = 2 * Y_bil[i] - 1
        paix[i] = 2 * x_bil[i] - 1
        paip[i] = 2 * p_bil[i] - 1
        pain[i] = 2 * n_bil[i] - 1
    cst = 0
    for i in lastAct:
        cst += sum(i)
    return m.ObjVal*(1-cst),[paiY, paip, pain, paix]


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
z_0d05 = 1.64
cp,cn,cy,co,rho,rhoY,TS,TB,Ct,Mt,It,Ot,Y0,p0,n0,x0,epsi1,delt1 = 15,30,6456,100,.2,.6,22.4167,1,807,538,896.7,201.75,1.02,24,0,0,.84,63
Y0,p0,n0,x0 = f2bL(Y0),f2bL(p0),f2bL(n0),f2bL(x0) # for SDDiP
sps = 3
Max_Ite = 300
nss = 5
TRJ = 1
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
    for o in range(sps):
        cutQt[t][o] = [0 for i in range(Max_Ite)]
for t in range(2,nss+1):
    for o in range(sps):
        for c in range(Max_Ite):
            cutQt[t][o][c]=[0,0]
# cutQt[t][o][ite][0][0] = cst
# cutQt[t][o][ite][0][1][0] = paiY
# cutQt[t][o][ite][0][1][1] = paip
# cutQt[t][o][ite][0][1][2] = pain
# cutQt[t][o][ite][0][1][3] = paix
# 4th []-> 0: sbcut, 1: intcut

rsu = [0 for i in range(nss+1)]
for t in range(2,nss+1):
    rsu[t] = [0 for i in range(TRJ)]
o = sps*np.ones((TRJ,nss+1),dtype=np.int32)

for ite in range(1, Max_Ite):
    t = 1
    tmp, rsu[1], lb = trial(0, ite, t, 0, epsi1, delt1, [Y0, p0, n0, x0])
    trjVal = tmp * np.ones(TRJ,dtype=np.float64) # t=1 value
    for trj in range(TRJ): # trj = 0
        t = 2
        o[trj,t] = rpd2o(r(), P[1])
        tmp,rsu[t][trj],_ = trial(0, ite, t, o[trj,t],epsi[o[trj,t],t], delt[o[trj,t],t], rsu[t-1])
        trjVal[trj] += tmp
        for t in range(3,nss+1):
            o[trj,t] = rpd2o(r(),P[t-1][:,o[trj,t-1]])
            tmp,rsu[t][trj],_ = trial(0, ite, t, o[trj,t],epsi[o[trj,t],t], delt[o[trj,t],t], rsu[t-1][trj])
            trjVal[trj] += tmp
    ub = trjVal[0]
    print('%8d | %8g | %8g' % (ite, ub, lb))
    # -------------------backward ite : purely gen cuts -------------------
    for t in range(5,2,-1):
        for trj in range(TRJ):
            for ot in range(sps):
                cutQt[t][ot][ite][0] = bwTrial(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1][trj])
                cutQt[t][ot][ite][1] = bwintCut(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1][trj])
    t = 2
    for ot in range(sps):  # this-stage-all-possible-scenes
        cutQt[t][ot][ite][0] = bwTrial(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1])
        cutQt[t][ot][ite][1] = bwintCut(1,ite,t,ot,epsi[ot,t],delt[ot,t],rsu[t-1])

