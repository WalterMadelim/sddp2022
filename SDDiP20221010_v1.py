import gurobipy as gp
from gurobipy import GRB
import numpy as np
# waterNet 4-stage schedule: SDDiP
# use sB and Lag cut, convergent after 20 iterations (see the bottom)
# speed is satisfied
# gap is closed
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

lW = np.array([324,395,432,373])
Pe = np.array([vP,vP,fP,pP])

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
nss = 4
MaxIte = 50000
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
                if maxerVal + 1e-4 <= valt:
                    print('gap = %8g' % (valt-maxerVal))
                break
            pai += 19000*drc
# 190 is perfect

        # if vldCut != (nss-1)*ite + nss - t: # for lag cuts
            # print('one invalid cut may be added ______>______>______>______>')

# htouLevel:Set parameter Username
# Academic license - for non-commercial use only - expires 2023-07-26
#   0:      98|  1:      92|  2:      84|  3:     100|ite=       0\lb=0.877276< 10.6666=0.877276 + 0.849475 +   1.7049 +  7.23492
# htouLevel:  0:      99|  1:      93|  2:      99|  3:     100|ite=       1\lb= 3.43426< 9.35139=0.939056 + 0.849475 +  3.45994 +  4.10292
# htouLevel:  0:      98|  1:      94|  2:     112|  3:     107|ite=       2\lb=    3.46< 9.66469=0.877276 + 0.973035 +  4.96426 +  2.85012
# gap = 0.061457
# htouLevel:  0:      98|  1:      92|  2:     100|  3:     100|ite=       3\lb= 6.18222< 9.33153=0.877276 + 0.849475 +  3.71066 +  3.89412
# htouLevel:  0:      99|  1:      95|  2:     104|  3:     100|ite=       4\lb= 6.21826< 8.80703=0.939056 + 0.973035 +  3.83602 +  3.05892
# htouLevel:  0:     100|  1:      96|  2:     105|  3:     100|ite=       5\lb= 6.28004< 8.66001= 1.00084 + 0.973035 +  3.83602 +  2.85012
# htouLevel:  0:     101|  1:      97|  2:     105|  3:     100|ite=       6\lb= 6.34182< 8.59643= 1.06262 + 0.973035 +  3.71066 +  2.85012
# htouLevel:  0:      98|  1:      98|  2:     105|  3:     100|ite=       7\lb= 6.36756< 8.53285=0.877276 +  1.22016 +   3.5853 +  2.85012
# htouLevel:  0:     102|  1:      99|  2:     105|  3:     100|ite=       8\lb=  6.4036< 8.46927=  1.1244 +  1.03482 +  3.45994 +  2.85012
# htouLevel:  0:      99|  1:     100|  2:     105|  3:     100|ite=       9\lb= 6.42934< 8.40569=0.939056 +  1.28193 +  3.33458 +  2.85012
# htouLevel:  0:     103|  1:     101|  2:     105|  3:     100|ite=      10\lb= 6.46538< 8.34211= 1.18618 +  1.09659 +  3.20922 +  2.85012
# htouLevel:  0:     100|  1:     102|  2:     105|  3:     100|ite=      11\lb= 6.49112< 8.27853= 1.00084 +  1.34372 +  3.08386 +  2.85012
# htouLevel:  0:     104|  1:     103|  2:     105|  3:     100|ite=      12\lb= 6.52716< 8.21495= 1.24796 +  1.15838 +   2.9585 +  2.85012
# htouLevel:  0:     101|  1:     104|  2:     105|  3:     100|ite=      13\lb=  6.5529< 8.15137= 1.06262 +  1.40549 +  2.83314 +  2.85012
# htouLevel:  0:     105|  1:     105|  2:     105|  3:     100|ite=      14\lb= 6.58894< 8.08779= 1.30974 +  1.22016 +  2.70778 +  2.85012
# htouLevel:  0:      98|  1:     106|  2:     105|  3:     100|ite=      15\lb= 6.61468< 8.02421=0.877276 +   1.7144 +  2.58242 +  2.85012
# htouLevel:  0:     106|  1:     107|  2:     105|  3:     100|ite=      16\lb= 6.65072< 7.96063= 1.37152 +  1.28193 +  2.45706 +  2.85012
# htouLevel:  0:     102|  1:     108|  2:     105|  3:     100|ite=      17\lb= 6.67646< 7.89705=  1.1244 +  1.59084 +   2.3317 +  2.85012
# htouLevel:  0:     107|  1:     109|  2:     105|  3:     100|ite=      18\lb=  6.7125< 7.83347=  1.4333 +  1.34372 +  2.20634 +  2.85012
# htouLevel:  0:      99|  1:     110|  2:     105|  3:     100|ite=      19\lb= 6.73824< 7.76989=0.939056 +  1.89974 +  2.08098 +  2.85012
# htouLevel:  0:     108|  1:     111|  2:     105|  3:     100|ite=      20\lb= 6.77428< 7.70631= 1.49508 +   1.4055 +  1.95562 +  2.85012
# htouLevel:  0:     103|  1:     112|  2:     105|  3:     100|ite=      21\lb= 6.80002< 7.64273= 1.18618 +  1.77618 +  1.83026 +  2.85012
# htouLevel:  0:     109|  1:     113|  2:     105|  3:     100|ite=      22\lb= 6.83606< 7.57915= 1.55686 +  1.46728 +   1.7049 +  2.85012
# htouLevel:  0:     100|  1:     114|  2:     106|  3:     101|ite=      23\lb=  6.8618< 7.64093= 1.00084 +  2.08507 +   1.7049 +  2.85012
# htouLevel:  0:     110|  1:     113|  2:     105|  3:     100|ite=      24\lb= 6.89784< 7.57915= 1.61864 +   1.4055 +   1.7049 +  2.85012
# htouLevel:  0:     104|  1:     113|  2:     105|  3:     100|ite=      25\lb= 6.92358< 7.57915= 1.24796 +  1.77618 +   1.7049 +  2.85012
# htouLevel:  0:     111|  1:     113|  2:     105|  3:     100|ite=      26\lb= 6.95962< 7.57915= 1.68042 +  1.34371 +   1.7049 +  2.85012
# htouLevel:  0:     101|  1:     113|  2:     105|  3:     100|ite=      27\lb= 6.98536< 7.57915= 1.06262 +  1.96152 +   1.7049 +  2.85012
# htouLevel:  0:     112|  1:     113|  2:     105|  3:     100|ite=      28\lb=  7.0214< 7.57915=  1.7422 +  1.28193 +   1.7049 +  2.85012
# htouLevel:  0:     105|  1:     113|  2:     105|  3:     100|ite=      29\lb= 7.04714< 7.57915= 1.30974 +  1.71439 +   1.7049 +  2.85012
# htouLevel:  0:     113|  1:     113|  2:     105|  3:     100|ite=      30\lb= 7.08318< 7.57915= 1.80398 +  1.22016 +   1.7049 +  2.85012
# htouLevel:  0:      98|  1:     113|  2:     105|  3:     100|ite=      31\lb= 7.10892< 7.57915=0.877276 +  2.14685 +   1.7049 +  2.85012
# htouLevel:  0:     114|  1:     113|  2:     105|  3:     100|ite=      32\lb= 7.14496< 7.57915= 1.86576 +  1.15837 +   1.7049 +  2.85012
# htouLevel:  0:     106|  1:     113|  2:     105|  3:     100|ite=      33\lb=  7.1707< 7.57915= 1.37152 +  1.65261 +   1.7049 +  2.85012
# htouLevel:  0:     115|  1:     113|  2:     105|  3:     100|ite=      34\lb= 7.20674< 7.57915= 1.92754 +  1.09659 +   1.7049 +  2.85012
# htouLevel:  0:     102|  1:     113|  2:     105|  3:     100|ite=      35\lb= 7.23248< 7.57915=  1.1244 +  1.89973 +   1.7049 +  2.85012
# htouLevel:  0:     116|  1:     113|  2:     105|  3:     100|ite=      36\lb= 7.26852< 7.57915= 1.98932 +  1.03482 +   1.7049 +  2.85012
# htouLevel:  0:     107|  1:     113|  2:     105|  3:     100|ite=      37\lb= 7.29426< 7.57915=  1.4333 +  1.59084 +   1.7049 +  2.85012
# htouLevel:  0:     117|  1:     113|  2:     105|  3:     100|ite=      38\lb=  7.3303< 7.57915=  2.0511 + 0.973035 +   1.7049 +  2.85012
# htouLevel:  0:      99|  1:     113|  2:     105|  3:     100|ite=      39\lb= 7.35604< 7.57915=0.939056 +  2.08507 +   1.7049 +  2.85012
# htouLevel:  0:     118|  1:     113|  2:     105|  3:     100|ite=      40\lb= 7.39208< 7.57915= 2.11288 + 0.911255 +   1.7049 +  2.85012
# htouLevel:  0:     108|  1:     113|  2:     105|  3:     100|ite=      41\lb= 7.41782< 7.57915= 1.49508 +  1.52905 +   1.7049 +  2.85012
# htouLevel:  0:     119|  1:     113|  2:     105|  3:     100|ite=      42\lb= 7.45386< 7.57915= 2.17466 + 0.849475 +   1.7049 +  2.85012
# htouLevel:  0:     103|  1:     113|  2:     105|  3:     100|ite=      43\lb=  7.4796< 7.57915= 1.18618 +  1.83796 +   1.7049 +  2.85012
# htouLevel:  0:     120|  1:     115|  2:     107|  3:     102|ite=      44\lb= 7.51564< 7.70271= 2.23644 + 0.911255 +   1.7049 +  2.85012
# htouLevel:  0:     121|  1:     116|  2:     108|  3:     103|ite=      45\lb= 7.57742< 7.76449= 2.29822 + 0.911255 +   1.7049 +  2.85012
# htouLevel:  0:     100|  1:     113|  2:     105|  3:     100|ite=      46\lb= 7.57915< 7.57915= 1.00084 +   2.0233 +   1.7049 +  2.85012
# htouLevel:  0:     100|  1:     113|  2:     105|  3:     100|ite=      47\lb= 7.57915< 7.57915= 1.00084 +   2.0233 +   1.7049 +  2.85012
# htouLevel:  0:     100|  1:     113|  2:     105|  3:     100|ite=      48\lb= 7.57915< 7.57915= 1.00084 +   2.0233 +   1.7049 +  2.85012