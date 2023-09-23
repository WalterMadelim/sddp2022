import gurobipy as gp
from gurobipy import GRB
import numpy as np
# waterNet_schedule: can close the gap but ub do not decrease.
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

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            m.addConstr(tha>=cut[its,t,-1]+gsum(cut[its,t,i]*htou[i] for i in range(Bhtou)))
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
lW = np.array([368,395,395,401,380,358,357,319,314,343,357,335,334,324,316,319,331,341,338,331,398,422,385,352],dtype=np.float64)

tkCoef = 3600/20
vP, fP, pP = 308.9, 626.8, 1044 # valley,flat,peak price of ele (yuan/MWh)
Pe = np.array([vP,vP,vP,vP,vP,vP,vP,vP,pP,pP,pP,pP,fP,fP,fP,fP,fP,pP,pP,pP,pP,fP,fP,fP])
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
vldCutNum = 0
nss = 24
MaxIte = 20000
cut = np.zeros((MaxIte,nss-1,Bhtou+1))
act = np.zeros((nss,Bhtou))
val = np.zeros(nss)
for ite in range(MaxIte):
    dR = 0
    t = 0
    val[t],act[t],lb = wNet(dR,ite,t,act[t-1],0,0)
    for t in range(1,nss):
        val[t],act[t],_ = wNet(dR,ite,t,act[t-1],0,0)
    ub = sum(val)
    print('lb=%8g<%8g='%(lb,ub),end='')
    for i in range(nss-1):
        print('%8g + ' % (val[i]),end='')
    print('%8g' % (val[nss-1]))

    dR = 1
    for t in range(nss-1,0,-1):
        valt = wNet(dR,ite,t,act[t-1],0,0)
        pai = wNet(dR,ite,t,act[t-1],1,0)
        while 1:
            maxerVal, subVal, tailVal, drc, cpv = wNet(dR,ite,t,act[t-1],2,pai)
            # print('\nf(pai)=%8g=%8g+%8g' % (maxerVal,subVal,tailVal))
            # print('drc is:',drc,'\ncpv is:',cpv)
            if sum(abs(drc)) < SCLhtou / 2 ** 10:
                # print('subGrad = 0, opt Attained!')
                cut[ite,t-1] = np.append(pai, subVal)
                # check tightness
                # print('t=%8d: maxer = %8g|stdVal = %8g' % (t,maxerVal,valt))
                break
            pai += drc



