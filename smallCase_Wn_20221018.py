import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
# 4-stage, small-case, water_only 20221018
# model setting:
# pump power level:
# t=       0| 546.667|     160|
# t=       1| 196.176| 240.645|
# t=       2| 255.501|  255.61|
# t=       3| 198.373| 239.455|
def gsum(a):
    return gp.quicksum(a)
def gdot(a,b):
    return gp.quicksum(a[i]*b[i] for i in range(len(a)))
def wNet(dR,ite,t,htoulast,s,pai):
    m = gp.Model('water_t')
    if t != 0:
        zl = len(htoulast)  # length of copy vector = length of lastAction
        z = m.addVars(zl, ub=1)
        if s != 2:
            m.addConstrs(z[i] == htoulast[i] for i in range(zl))
    # -----model part-----
    s_pump = m.addVars(pumpNum,4,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # represent 4 finite states
    q_pump,h_pump,p_pump = m.addVars(pumpNum,lb=100,ub=1100),m.addVars(pumpNum,lb=86,ub=229),m.addVars(pumpNum,lb=160,ub=556)
    l_pump = m.addVars(pumpNum, ub=1)

    s_pipe = m.addVars(pipeNum,3,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
    q_pipe,h_pipe = m.addVars(pipeNum,ub=1500),m.addVars(pipeNum,ub=704)
    l0_pipe = m.addVars(pipeNum, ub=1)
    l1_pipe = m.addVars(pipeNum)
    hrd = m.addVars(hrdNum) # RDV
    htin = m.addVar()
    htou = m.addVars(Bhtou,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # discretized state variable
    m.addConstr(gdot(fbhtou,htou) == (gdot(fbhtou,z) if t>0 else 35) + tkCoef * (q_pipe[2] - q_pipe[3]))
    m.addConstr(gdot(fbhtou,htou) <= 70)
    for n in range(pumpNum):
        m.addConstr(gsum(s_pump[n, i] for i in range(4)) == 1)
        m.addConstr(l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000*(1-s_pump[n, 3]))
        m.addConstr(l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000*(1-s_pump[n, 3]))
        m.addConstr(l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143*(1-s_pump[n, 3]))
        m.addConstr(l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143*(1-s_pump[n, 3]))
        m.addConstr(l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396*(1-s_pump[n, 0]))
        m.addConstr(l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396*(1-s_pump[n, 1]))
        m.addConstr(l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396*(1-s_pump[n, 2]))
        m.addConstr(l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396*(1-s_pump[n, 3]))
        m.addConstr(l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396*(1-s_pump[n, 3]))
    for n in range(pipeNum):
        m.addConstr(gsum(s_pipe[n, i] for i in range(3)) == 1)
        m.addConstr(l0_pipe[n] * qt_pipe[0] + (1 - l0_pipe[n]) * qt_pipe[1] - q_pipe[n] >= -1500 * (1 - s_pipe[n, 0]))
        m.addConstr(l0_pipe[n] * qt_pipe[0] + (1 - l0_pipe[n]) * qt_pipe[1] - q_pipe[n] <= 1500 * (1 - s_pipe[n, 0]))
        m.addConstr(
            l0_pipe[n] * ht_pipe[n, 0] + (1 - l0_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] >= -704 * (1 - s_pipe[n, 0]))
        m.addConstr(
            l0_pipe[n] * ht_pipe[n, 0] + (1 - l0_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] <= 704 * (1 - s_pipe[n, 0]))
        m.addConstr(l0_pipe[n] * qt_pipe[1] + (1 - l0_pipe[n]) * qt_pipe[2] - q_pipe[n] >= -1500 * (1 - s_pipe[n, 1]))
        m.addConstr(l0_pipe[n] * qt_pipe[1] + (1 - l0_pipe[n]) * qt_pipe[2] - q_pipe[n] <= 1500 * (1 - s_pipe[n, 1]))
        m.addConstr(
            l0_pipe[n] * ht_pipe[n, 1] + (1 - l0_pipe[n]) * ht_pipe[n, 2] - h_pipe[n] >= -704 * (1 - s_pipe[n, 1]))
        m.addConstr(
            l0_pipe[n] * ht_pipe[n, 1] + (1 - l0_pipe[n]) * ht_pipe[n, 2] - h_pipe[n] <= 704 * (1 - s_pipe[n, 1]))
        m.addConstr(l1_pipe[n] * qt_pipe[3] + (1 - l1_pipe[n]) * qt_pipe[2] - q_pipe[n] >= -1500 * (1 - s_pipe[n, 2]))
        m.addConstr(l1_pipe[n] * qt_pipe[3] + (1 - l1_pipe[n]) * qt_pipe[2] - q_pipe[n] <= 1500 * (1 - s_pipe[n, 2]))
        m.addConstr(
            l1_pipe[n] * ht_pipe[n, 3] + (1 - l1_pipe[n]) * ht_pipe[n, 2] - h_pipe[n] >= -704 * (1 - s_pipe[n, 2]))
        m.addConstr(
            l1_pipe[n] * ht_pipe[n, 3] + (1 - l1_pipe[n]) * ht_pipe[n, 2] - h_pipe[n] <= 704 * (1 - s_pipe[n, 2]))

    # flow constrs
    m.addConstr(q_pump[0] == q_pipe[2] + q_pipe[0]) # node 0
    m.addConstr(q_pipe[0]+q_pump[1] == q_pipe[1]) # node 1
    m.addConstr(q_pipe[1] == q_pipe[5]) # node 2
    m.addConstr(q_pipe[4]+q_pipe[5] == lW[t])
    m.addConstr(q_pipe[4] == q_pipe[3]) # node 4

    # pressure constrs
    m.addConstr(h_pump[0]-h_pipe[2] == htin)# node tank
    m.addConstr(gdot(fbhtou,htou) - h_pipe[3] - h_pipe[4] - hrd[2] == 0)# node 6
    m.addConstr(h_pump[0] - h_pipe[0] - hrd[0] == h_pump[1])
    m.addConstr(h_pump[1] - h_pipe[1] - h_pipe[5] - hrd[1] == 0)

    activePM = 1e-3*gsum(p_pump[p] for p in range(pumpNum))
    theObj = Pe[t]*activePM*(1 + np.tan(np.arccos(.85)))
    if t == nss-1: # penalty for not-satisfying daily balance
        endPen = m.addVar()
        m.addConstr(endPen >= 35-gdot(fbhtou,htou))
        theObj += 500*endPen

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            m.addConstr(tha>=cut[its,t,-1]+gdot(cut[its,t,:-1],htou))
    theObj += tha
    # -----model part-----
    if s == 2:
        theObj -= gsum(pai[i]*z[i] for i in range(zl))
    m.setObjective(theObj)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        exit(3)
    if dR == 0:
        if ite == 88:
            if t == nss-1:
                print('endPen',endPen.X)
        a = np.zeros(Bhtou)
        for i in range(Bhtou):
            a[i] = htou[i].X
        return m.ObjVal-tha.X,a,m.ObjVal
    else:
        if s == 0:
            return m.ObjVal
        elif s == 1:
            pai = np.zeros(zl)
            l = m.getConstrs()
            for i in range(zl):
                pai[i] = l[i].Pi
            return pai
        elif s == 2:
            subVal = m.ObjVal
            tailVal = np.dot(pai,htoulast)
            maxerVal = subVal + tailVal
            cpv = np.zeros(zl)
            for i in range(zl):
                cpv[i] = z[i].X
            drc = htoulast - cpv
            return maxerVal,subVal,tailVal,drc,cpv
def fbgen(lb,ub):
    fb = np.zeros(40)
    for i in range(-20,20):
        fb[i+20] = 2**i
    fb = fb[(fb<ub)&(fb>lb)]
    return fb,len(fb),fb[0]
def make0(c): # avoid the gurobi small constr warning
    c[abs(c) < 5e-12] = 0
    return c
def cutGen(dR,ite,t,al):
    global cut
    valt = wNet(dR, ite, t, al, 0, 0)
    pai = wNet(dR, ite, t, al, 1, 0)
    if ite <= 1:
        maxerVal, subVal, tailVal, drc, cpv = wNet(1, ite, t, al, 2, pai)
        cut[ite, t - 1] = np.append(pai, subVal)
    else:
        while 1:
            maxerVal, subVal, _, drc, _ = wNet(1, ite, t, al, 2, pai)
            print('\r=%8g|%8g=' % (maxerVal, valt), end='')
            if maxerVal + 1e-4 >= valt or sum(abs(drc)) < 1e-7 or valt < maxerVal:
                if abs(maxerVal - valt) >= 1e-4:
                    print('Gap = %8g' % (valt - maxerVal))
                cut[ite, t - 1] = make0(np.append(pai, subVal))
                break
            pai += (valt - maxerVal) / np.linalg.norm(drc) ** 2 * drc

ppKf = 1.016 # ppK = ppKf * (km) * (m) ** -5, length, diameter
l = np.array([3.7,.77,.521,.521,.521,.77])
d = np.array([.3,.25,.25,.25,.25,.25])
ppK = ppKf * l * d ** -5

qt_pipe = np.array([0,100,500,1100],dtype=np.float64)/3600
ht_pipe = np.zeros((len(ppK),4),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
qt_pipe = np.array([0,100,500,1100],dtype=np.float64)

qt_pump = np.array([[1,7,11],[1,5.5,8.8]],dtype=np.float64)*100 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]],dtype=np.float64)
pt_pump = np.array([[330,500,556],[160,233,262]],dtype=np.float64) # kW

fbhtou,Bhtou,SCLhtou = fbgen(.5,128)


vP, fP, pP = 308.9, 626.8, 1044 # valley,flat,peak price of ele (yuan/MWh)
Pe = np.array([vP,pP,fP,pP])
vW, fW, pW = 200,1000,1600 # 200, 1000, 2190
lW = np.array([vW,pW,fW,pW])

tkCoef = 6/160 # delta t = 6 hours

pipeNum = 6
pumpNum = 2
hrdNum = 3

nss = 4
MaxIte = 50000
cut = np.zeros((MaxIte,nss-1,Bhtou+1)) # last stage do not need a cut
act = np.zeros((nss,Bhtou))
val,vlb = np.zeros(nss),np.zeros(nss)
bestub = 40000
for ite in range(MaxIte):
    dR = 0
    for t in range(nss):
        val[t],act[t],vlb[t] = wNet(dR,ite,t,act[t-1],0,0)
        # print(np.dot(act[t],fbhtou))
    # ----------- statistic -----------
    lb = vlb[0]
    ub = sum(val)
    if bestub > ub:
        bestIte = ite
        bestub = ub
        bestact = copy.deepcopy(act)
    print('ite=%8d\lb=%12g<%12g<%12g=' % (ite, lb, bestub, ub), end='')
    if lb + 1 >= bestub:  # final convergent
        print('converge at ite = %8d' % ite)
        print('best at bestIte = %8d' % bestIte)
        break
    for i in range(nss-1):
        print('%8g + ' % (val[i]),end='')
    print('%8g' % (val[nss-1]))
    # ----------- statistic -----------
    dR = 1
    for t in range(nss-1,0,-1):
        cutGen(dR, ite, t, act[t-1])
