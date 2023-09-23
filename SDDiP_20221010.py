import gurobipy as gp
from gurobipy import GRB
import numpy as np
# waterNet 3-stage schedule: SDDiP
# use sB and Lag cut, convergent after 20 iterations (see the bottom)
# speed is satisfied
# 20221010
# @doorvanbei
def gsum(a):
    return gp.quicksum(a)
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
    s_pipe = m.addVars(pipesNum,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # represent 3 finite states
    q_pipe,h_pipe = m.addVars(pipesNum,ub=1500),m.addVars(pipesNum,ub=704)
    l_pipe = m.addVars(pipesNum)
    hrd = m.addVars(hrdNum) # RDV
    h0, h1, h2, h3, h4 = m.addVar(),m.addVar(),m.addVar(),m.addVar(),m.addVar()
    htin = m.addVar()
    htou = m.addVars(Bhtou,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
    if t == 0:
        m.addConstr(gsum(fbhtou[i]*htou[i] for i in range(Bhtou)) == 100 + tkCoef * (q_pipe[2] - q_pipe[3]))
    else:
        m.addConstr(gsum(fbhtou[i]*htou[i] for i in range(Bhtou)) == gsum(fbhtou[i]*z[i] for i in range(Bhtou)) + tkCoef * (q_pipe[2] - q_pipe[3]))

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
    for n in range(pipesNum):
        m.addConstr(l_pipe[n] * qt_pipe[0] + (1 - l_pipe[n]) * qt_pipe[1] - q_pipe[n] >= -1500*s_pipe[n])
        m.addConstr(l_pipe[n] * qt_pipe[0] + (1 - l_pipe[n]) * qt_pipe[1] - q_pipe[n] <= 1500*s_pipe[n])
        m.addConstr(l_pipe[n]*qt_pipe[2]+(1-l_pipe[n])*qt_pipe[1] - q_pipe[n] >= -1500*(1-s_pipe[n]))
        m.addConstr(l_pipe[n]*qt_pipe[2]+(1-l_pipe[n])*qt_pipe[1] - q_pipe[n] <= 1500*(1-s_pipe[n]))
        m.addConstr(l_pipe[n] * ht_pipe[n, 0] + (1 - l_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] >= -704 * s_pipe[n])
        m.addConstr(l_pipe[n] * ht_pipe[n, 0] + (1 - l_pipe[n]) * ht_pipe[n, 1] - h_pipe[n] <= 704 * s_pipe[n])
        m.addConstr(l_pipe[n]*ht_pipe[n,2]+(1-l_pipe[n])*ht_pipe[n,1] - h_pipe[n] >= -704*(1-s_pipe[n]))
        m.addConstr(l_pipe[n]*ht_pipe[n,2]+(1-l_pipe[n])*ht_pipe[n,1] - h_pipe[n] <= 704*(1-s_pipe[n]))
    # flow constrs
    m.addConstr(q_pump[0] == q_pipe[2] + q_pipe[0]) # node 0
    m.addConstr(q_pipe[0]+q_pump[1] == q_pipe[1]) # node 1
    m.addConstr(q_pipe[1] == q_pipe[5]) # node 2
    m.addConstr(q_pipe[4]+q_pipe[5] == lW[t]) # node 3: load level: 140~1000
    m.addConstr(q_pipe[4] == q_pipe[3]) # node 4
    # pressure constrs
    m.addConstr(h_pump[0] == h0) # node 0_with pump
    m.addConstr(h_pump[1] == h1) # node 1_with pump
    m.addConstr(h0-hrd[0]-h_pipe[0] == h1)# node 1
    m.addConstr(h1-h_pipe[1] == h2)# node 2
    m.addConstr(h2 - h_pipe[5] - hrd[1] == h3)# node 5
    m.addConstr(h0-h_pipe[2] == htin)# node tank
    m.addConstr(gsum(fbhtou[i]*htou[i] for i in range(Bhtou)) - h_pipe[3] == h4)# node 6
    m.addConstr(h4 - h_pipe[4] - hrd[2] == h3)# node 5

    activePM = 1e-3*gsum(p_pump[p] for p in range(pumpNum))
    theObj = 10 * gsum(hrd[i] for i in range(hrdNum)) + Pe[t]*activePM*(1 + np.tan(np.arccos(.85)))
    theObj = 1e-5*Pe[t]*(q_pump[0] + q_pump[1])
    if t == nss-1: # penalty for not-satisfying daily balance
        endPen = m.addVar()
        m.addConstr(endPen >= 100-gsum(fbhtou[i]*htou[i] for i in range(Bhtou)))
        theObj += .3*endPen

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            m.addConstr(tha>=cut[its,t,0,-1]+gsum(cut[its,t,0,i]*htou[i] for i in range(Bhtou)))
            m.addConstr(tha>=cut[its,t,1,-1]+gsum(cut[its,t,1,i]*htou[i] for i in range(Bhtou)))
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
    # print('\nThe ObjVal is %g' % m.ObjVal)
    # print('\nThe tank water level')
    # print('%8d|%8g' % (t, htou.X))
    # print('\nThe lower|upper flow')
    # print('%8d|%8g|%8g' % (t,q_pipe[4].X,q_pipe[5].X))
    # print('\nThe pump0|pump1 w level')
    # print('%8d:%3d|%3d|%3d|%3d| |%3d|%3d|%3d|%3d' % (
    # t, s_pump[ 0, 0].X, s_pump[ 0, 1].X, s_pump[ 0, 2].X, s_pump[ 0, 3].X, s_pump[ 1, 0].X,
    # s_pump[ 1, 1].X, s_pump[ 1, 2].X, s_pump[ 1, 3].X))
    # print('%8d:%8g|%8g (kW)' % (t, p_pump[ 0].X, p_pump[ 1].X))
    # print('%8d:%8g|%8g (m)' % (t, h_pump[ 0].X, h_pump[ 1].X))
    # print('%8d:%8g|%8g (m3/h)' % (t, q_pump[ 0].X, q_pump[ 1].X))
    # print('\nThe hrd')
    # print('%8d|%12g|%12g|%12g' % (t, hrd[ 0].X, hrd[ 1].X, hrd[ 2].X))
    if dR == 0:
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
fbhtou,Bhtou,SCLhtou = fbgen(.5,1000)
vP, fP, pP = 308.9, 626.8, 1044 # valley,flat,peak price of ele (yuan/MWh)
lW = np.array([368,395,395,401,380,358,357,319,314,343,357,335,334,324,316,319,331,341,338,331,398,422,385,352],dtype=np.float64)
Pe = np.array([vP,vP,vP,vP,vP,vP,vP,vP,pP,pP,pP,pP,fP,fP,fP,fP,fP,pP,pP,pP,pP,fP,fP,fP])

lW = np.array([324,395,432])
Pe = np.array([vP,fP,pP])

tkCoef = 1/20
qt_pump = np.array([[1,7,11],[1,5.5,8.8]],dtype=np.float64)*100 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]],dtype=np.float64)
pt_pump = np.array([[330,500,556],[160,233,262]],dtype=np.float64) # kW

qt_pipe = np.array([100,500,1000],dtype=np.float64)/3600
ppK = np.array([1547.2653,1041.4285,833.1428,5207.1425,3644.99975,2395.28555],dtype=np.float64)# h=Kq**2, q in (m3/s), h in (m)
ht_pipe = np.zeros((len(ppK),3),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
qt_pipe = np.array([100,500,1000],dtype=np.float64)

pipesNum = 6
pumpNum = 2
hrdNum = 3
nss = 3
MaxIte = 30000
cut = np.zeros((MaxIte,nss-1,2,Bhtou+1)) # last stage do not need a cut
act = np.zeros((nss,Bhtou))
val = np.zeros(nss)
vlb = np.zeros(nss)
vldCut = 0
for ite in range(MaxIte):
    dR = 0
    print('htouLevel:',end='')
    for t in range(nss):
        val[t],act[t],vlb[t] = wNet(dR,ite,t,act[t-1],0,0)
        print('%3d:%8g'%(t,np.dot(fbhtou,act[t])),end='|')
    lb = vlb[0]
    ub = sum(val)
    print('ite=%8d\lb=%8g<%8g='%(ite,lb,ub),end='')
    for i in range(nss-1):
        print('%8g + ' % (val[i]),end='')
    print('%8g' % (val[nss-1]))
    dR = 1
    for t in range(nss-1,0,-1):
        valt = wNet(dR,ite,t,act[t-1],0,0)
        # print(act[t-1])
        pai = wNet(dR,ite,t,act[t-1],1,0)
        maxerVal, subVal, tailVal, drc, cpv = wNet(dR, ite, t, act[t - 1], 2, pai)
        cut[ite, t - 1, 0] = np.append(pai, subVal)  # sB cut
        while 1:
            maxerVal, subVal, tailVal, drc, cpv = wNet(dR,ite,t,act[t-1],2,pai)
            # print('\nf(pai)=%8g=%8g+%8g' % (maxerVal,subVal,tailVal))
            # print('drc is:',drc,'\ncpv is:',cpv)
            # print('t=%8d: maxer = %8g|stdVal = %8g' % (t,maxerVal,valt))
            if sum(abs(drc)) < SCLhtou / 2 ** 10:# or maxerVal + 1e-4 >= valt:
                # print('subGrad = 0, opt Attained!')
                cut[ite,t-1,1] = np.append(pai, subVal)
                # check tightness
                # print('t=%8d: maxer = %8g|stdVal = %8g' % (t,maxerVal,valt))
                vldCut += 1
                break
            pai += 10*drc

        if vldCut != 2*ite + nss - t: # for lag cuts
            print('one invalid cut may be added ______>______>______>______>')

#
# # convergent curve:
# htouLevel:  0:      98|  1:      92|  2:     100|ite=       0\lb=0.877276< 8.78146=0.877276 +   1.7237 +  6.18048
# htouLevel:  0:      99|  1:      96|  2:     100|ite=       1\lb= 2.61052< 8.38412=0.939056 +  2.09978 +  5.34528
# htouLevel:  0:     100|  1:     122|  2:     114|ite=       2\lb=  2.6723<  9.0743= 1.00084 +  5.23378 +  2.83968
# htouLevel:  0:      98|  1:     108|  2:     100|ite=       3\lb= 5.33274< 7.44642=0.877276 +  3.72946 +  2.83968
# htouLevel:  0:     101|  1:     108|  2:     100|ite=       4\lb= 5.51808< 7.25568= 1.06262 +  3.35338 +  2.83968
# htouLevel:  0:     102|  1:     108|  2:     100|ite=       5\lb= 5.57986<  7.1921=  1.1244 +  3.22802 +  2.83968
# htouLevel:  0:     103|  1:     108|  2:     100|ite=       6\lb= 5.64164< 7.12852= 1.18618 +  3.10266 +  2.83968
# htouLevel:  0:     104|  1:     108|  2:     100|ite=       7\lb= 5.70342< 7.06494= 1.24796 +   2.9773 +  2.83968
# htouLevel:  0:     105|  1:     108|  2:     100|ite=       8\lb=  5.7652< 7.00136= 1.30974 +  2.85194 +  2.83968
# htouLevel:  0:     106|  1:     108|  2:     100|ite=       9\lb= 5.82698< 6.93778= 1.37152 +  2.72658 +  2.83968
# htouLevel:  0:     107|  1:     108|  2:     100|ite=      10\lb= 5.88876<  6.8742=  1.4333 +  2.60122 +  2.83968
# htouLevel:  0:     109|  1:     108|  2:     100|ite=      11\lb= 6.01232< 6.74704= 1.55686 +   2.3505 +  2.83968
# htouLevel:  0:     108|  1:     108|  2:     100|ite=      12\lb= 6.04266< 6.81062= 1.49508 +  2.47586 +  2.83968
# htouLevel:  0:     110|  1:     108|  2:     100|ite=      13\lb=  6.0741< 6.68346= 1.61864 +  2.22514 +  2.83968
# htouLevel:  0:     111|  1:     108|  2:     100|ite=      14\lb= 6.13588< 6.61988= 1.68042 +  2.09978 +  2.83968
# htouLevel:  0:     112|  1:     108|  2:     100|ite=      15\lb= 6.19766<  6.5563=  1.7422 +  1.97442 +  2.83968
# htouLevel:  0:     113|  1:     108|  2:     100|ite=      16\lb= 6.25944< 6.49272= 1.80398 +  1.84906 +  2.83968
# htouLevel:  0:      99|  1:     108|  2:     100|ite=      17\lb= 6.29812< 7.38284=0.939056 +   3.6041 +  2.83968
# htouLevel:  0:     114|  1:     108|  2:     100|ite=      18\lb= 6.32122< 6.42914= 1.86576 +   1.7237 +  2.83968
# htouLevel:  0:     115|  1:     109|  2:     101|ite=      19\lb=   6.383< 6.49092= 1.92754 +   1.7237 +  2.83968
# htouLevel:  0:     114|  1:     108|  2:     100|ite=      20\lb=  6.4001< 6.42914= 1.86576 +   1.7237 +  2.83968