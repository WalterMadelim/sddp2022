import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
# small-case: ele and heat only
# 20221014

def gsum(a):
    return gp.quicksum(a)
def gdot(a,b):
    return gp.quicksum(a[i]*b[i] for i in range(len(a)))
def fmac(moc,Tam):
    return 0.322 + 0.5621 * moc -2.209e-5 * Tam ** 2 + 0.005004 * Tam
def fpsc(mac,Tam):
    return 0.1672 + 0.8257 * mac - 1.006e-5 * Tam ** 2 + 0.0005409 * Tam
def fmadlb(mod): # do 1-dim pwl
    return 1.9094 * mod ** 4 -8.5877 * mod ** 3 + 14.5357 * mod ** 2 -10.9975 * mod + 4
def fpsd(mad,mod): # do 2-dim pwl
    return -8.05+21.2095*mad+1.7761*mod -20.6123*mad**2 +3.3126*mad*mod -2.754*mod**2 +6.8186*mad**3 -1.0925*mod*mad**2 -0.4445*mad*mod**2 +0.8743*mod**3
def fhsh(moh):
    return 0.31265*moh+0.7
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
    valt = cacase(dR, ite, t, al, 0, 0)
    pai = cacase(dR, ite, t, al, 1, 0)
    if ite <= 27:
        maxerVal, subVal, tailVal, drc, cpv = cacase(1, ite, t, al, 2, pai)
        cut[ite, t - 1] = np.append(pai, subVal)
    else:
        while 1:
            maxerVal, subVal, _, drc, _ = cacase(1, ite, t, al, 2, pai)
            print('\r=%8g|%8g=' % (maxerVal, valt), end='')
            if maxerVal + 1e-4 >= valt or sum(abs(drc)) < 1e-7 or valt < maxerVal:
                if abs(maxerVal - valt) >= 1e-4:
                    print('Gap = %8g' % (valt - maxerVal))
                cut[ite, t - 1] = make0(np.append(pai, subVal))
                break
            pai += (valt - maxerVal) / np.linalg.norm(drc) ** 2 * drc

def cacase(dR,ite,t,al,s,pai):
    m = gp.Model('water_t')
    if t != 0:
        zl = len(al)  # length of copy vector = length of lastAction
        z = m.addVars(zl, ub=1)
        if s != 2:
            m.addConstrs(z[i] == al[i] for i in range(zl))
    # -----model part-----
    mu_chp = m.addVars(U_chp,vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)  # 'on' state indicator of chp
    P_chp, Q_chp, h_chp, cst_chp, csh_chp, ng_chp = m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp)
    for u in range(U_chp):
        # cost
        m.addConstr(beta_ng2e * ng_chp[u] == 6*(.1*P_chp[u] / eta_pchp + .1*h_chp[u] / eta_hchp))  # how much ng consumed by chp
        m.addConstr(csh_chp[u] >= ((shPchp*(1 - mu_chp[u])) if t==0 else (shPchp*(z[2*Bca+u] - mu_chp[u]))))  # only at t = 0, shut down cost of chp
        m.addConstr(cst_chp[u] >= (0 if t==0 else (stPchp*(mu_chp[u] - z[2*Bca+u]))))
        # operation
        m.addConstr(h_chp[u] <= hB_chp * mu_chp[u])
        m.addConstr(P_chp[u] <= pA_chp * mu_chp[u])
        m.addConstr(Q_chp[u] <= .5 * 3 ** .5 * pA_chp * mu_chp[u])  # m.addConstrs(np.cos((2*k-1)/6*np.pi) * P_chp[t,u] + np.sin((2*k-1)/6*np.pi) * Q_chp[t,u] - np.cos(np.pi/6) * S_chp[t,u] <= 0 for k in range(6))
        m.addConstr(Q_chp[u] + 3 ** 0.5 * (P_chp[u] - pA_chp) <= 0)  # with ub of Q
        m.addConstr(P_chp[u] - pA_chp - (pA_chp - pB_chp) / (hA_chp - hB_chp) * (h_chp[u] - hA_chp) <= 0)
        m.addConstr(P_chp[u] - pD_chp - (pD_chp - pC_chp) / (hD_chp - hC_chp) * (h_chp[u] - hD_chp) >= -5e3 * (1 - mu_chp[u]))
        m.addConstr(P_chp[u] - pB_chp - (pC_chp - pB_chp) / (hC_chp - hB_chp) * (h_chp[u] - hB_chp) >= 0)  # -5e3 * (1 - mu_chp[t,u]))
    # ------------------------------
    Ma, Mh = m.addVars(Bca, vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1), m.addVars(Bca, vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # Mamin + MaSCL * (gdot(fbca, Ma))
    mac, moc, mad, mod = m.addVar(lb=.9, ub=1), m.addVar(lb=.6, ub=1.3), m.addVar(lb=.7, ub=1), m.addVar(lb=.7, ub=1.3)
    moh, psc, psd, hsh = m.addVar(lb=.8707, ub=1), m.addVar(lb=.8, ub=1), m.addVar(lb=.4, ub=1), m.addVar(lb=.6, ub=1)
    mac_bar, moc_bar, mad_bar, mod_bar, moh_bar, psc_bar, psd_bar, hsh_bar = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    musc, musd, mush = m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1),m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1),m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
    if t == 0:
        m.addConstr(gdot(fbca, Ma) - 64 - 6*3600/MaSCL*(macB * mac_bar - madB * mad_bar) == [-1/2,1/2])
        m.addConstr(gdot(fbca,Mh) - 64 - 6 * 3600/MhSCL * (2 * mocB * moc_bar - 2 * modB * mod_bar - mohB * moh_bar) == [-1/2,1/2])
    else:
        m.addConstr(gdot(fbca, Ma) - gsum(fbca[i]*z[i] for i in range(Bca)) - 6*3600/MaSCL*(macB * mac_bar - madB * mad_bar) == [-1/2,1/2])
        m.addConstr(gdot(fbca,Mh) - gsum(fbca[i]*z[Bca+i] for i in range(Bca)) - 6 * 3600/MhSCL * (2 * mocB * moc_bar - 2 * modB * mod_bar - mohB * moh_bar) == [-1/2,1/2])
    m.addConstr(psc == fpsc(mac, Ta))
    m.addConstr(mac == fmac(moc, Ta))
    m.addConstr(hsh == fhsh(moh))
    xt_madlb = np.linspace(.7, 1.3, 5)
    yt_madlb = fmadlb(xt_madlb)
    l_madlb, f_madlb = m.addVar(ub=1), m.addVars(4, vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
    m.addConstr(gsum(f_madlb[r] for r in range(4)) == 1)
    for r in range(4):
        m.addConstr(l_madlb * xt_madlb[r] + (1 - l_madlb) * xt_madlb[r + 1] - mod <= .6 * (1 - f_madlb[r]))
        m.addConstr(l_madlb * xt_madlb[r] + (1 - l_madlb) * xt_madlb[r + 1] - mod >= -.6 * (1 - f_madlb[r]))
        m.addConstr(l_madlb * yt_madlb[r] + (1 - l_madlb) * yt_madlb[r + 1] - mad <= .237 * (1 - f_madlb[r]))
    xlb_psd, xub_psd, ylb_psd, yub_psd, segs_psd = .7, 1, .7, 1.3, 2  # 2d-pwl of fpsd(mad,mod)
    xt_psd, yt_psd = np.linspace(xlb_psd, xub_psd, segs_psd + 1), np.linspace(ylb_psd, yub_psd, segs_psd + 1)
    tmp1, tmp2 = np.meshgrid(yt_psd, xt_psd)  # this is an intermediate assign (temporary)
    zt_psd = fpsd(tmp1, tmp2)
    l_psd = m.addVars(4, ub=1)
    fx_psd, fy_psd = m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1), m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)
    m.addConstr(gsum(l_psd[i] for i in range(4)) == 1)
    m.addConstr((l_psd[0] + l_psd[2]) * xt_psd[0] + (l_psd[1] + l_psd[3]) * xt_psd[1] - mad <= .3 * fx_psd)
    m.addConstr((l_psd[0] + l_psd[2]) * xt_psd[0] + (l_psd[1] + l_psd[3]) * xt_psd[1] - mad >= -.3 * fx_psd)
    m.addConstr((l_psd[0] + l_psd[2]) * xt_psd[1] + (l_psd[1] + l_psd[3]) * xt_psd[2] - mad <= .3 * (1 - fx_psd))
    m.addConstr((l_psd[0] + l_psd[2]) * xt_psd[1] + (l_psd[1] + l_psd[3]) * xt_psd[2] - mad >= -.3 * (1 - fx_psd))
    m.addConstr((l_psd[0] + l_psd[1]) * yt_psd[0] + (l_psd[2] + l_psd[3]) * yt_psd[1] - mod <= .6 * fy_psd)
    m.addConstr((l_psd[0] + l_psd[1]) * yt_psd[0] + (l_psd[2] + l_psd[3]) * yt_psd[1] - mod >= -.6 * fy_psd)
    m.addConstr((l_psd[0] + l_psd[1]) * yt_psd[1] + (l_psd[2] + l_psd[3]) * yt_psd[2] - mod <= .6 * (1 - fy_psd))
    m.addConstr((l_psd[0] + l_psd[1]) * yt_psd[1] + (l_psd[2] + l_psd[3]) * yt_psd[2] - mod >= -.6 * (1 - fy_psd))
    m.addConstr(l_psd[0] * zt_psd[0, 0] + l_psd[1] * zt_psd[1, 0] + l_psd[2] * zt_psd[0, 1] + l_psd[3] * zt_psd[
        1, 1] - psd <= 1.125 * fx_psd + 1.125 * fy_psd)
    m.addConstr(l_psd[0] * zt_psd[0, 0] + l_psd[1] * zt_psd[1, 0] + l_psd[2] * zt_psd[0, 1] + l_psd[3] * zt_psd[
        1, 1] - psd >= -1.125 * fx_psd - 1.125 * fy_psd)
    m.addConstr(l_psd[0] * zt_psd[0, 1] + l_psd[1] * zt_psd[1, 1] + l_psd[2] * zt_psd[0, 2] + l_psd[3] * zt_psd[
        1, 2] - psd <= 1.125 * fx_psd + 1.125 * (1 - fy_psd))
    m.addConstr(l_psd[0] * zt_psd[0, 1] + l_psd[1] * zt_psd[1, 1] + l_psd[2] * zt_psd[0, 2] + l_psd[3] * zt_psd[
        1, 2] - psd >= -1.125 * fx_psd - 1.125 * (1 - fy_psd))
    m.addConstr(l_psd[0] * zt_psd[1, 0] + l_psd[1] * zt_psd[2, 0] + l_psd[2] * zt_psd[1, 1] + l_psd[3] * zt_psd[
        2, 1] - psd <= 1.125 * (1 - fx_psd) + 1.125 * fy_psd)
    m.addConstr(l_psd[0] * zt_psd[1, 0] + l_psd[1] * zt_psd[2, 0] + l_psd[2] * zt_psd[1, 1] + l_psd[3] * zt_psd[
        2, 1] - psd >= -1.125 * (1 - fx_psd) - 1.125 * fy_psd)
    m.addConstr(l_psd[0] * zt_psd[1, 1] + l_psd[1] * zt_psd[2, 1] + l_psd[2] * zt_psd[1, 2] + l_psd[3] * zt_psd[
        2, 2] - psd <= 1.125 * (1 - fx_psd) + 1.125 * (1 - fy_psd))
    m.addConstr(l_psd[0] * zt_psd[1, 1] + l_psd[1] * zt_psd[2, 1] + l_psd[2] * zt_psd[1, 2] + l_psd[3] * zt_psd[
        2, 2] - psd >= -1.125 * (1 - fx_psd) - 1.125 * (1 - fy_psd))

    m.addConstr(musc + musd <= 1)
    m.addConstr(mac_bar - mac <= 1 - musc)
    m.addConstr(mac_bar - mac >= -(1 - musc))
    m.addConstr(mac_bar <= musc)
    m.addConstr(moc_bar - moc <= 1.3 * (1 - musc))
    m.addConstr(moc_bar - moc >= -1.3 * (1 - musc))
    m.addConstr(moc_bar <= musc)
    m.addConstr(psc_bar - psc <= 1 - musc)
    m.addConstr(psc_bar - psc >= -(1 - musc))
    m.addConstr(psc_bar <= musc)
    m.addConstr(mad_bar - mad <= 1 - musd)
    m.addConstr(mad_bar - mad >= -(1 - musd))
    m.addConstr(mad_bar <= musd)
    m.addConstr(mod_bar - mod <= 1.3 * (1 - musd))
    m.addConstr(mod_bar - mod >= -1.3 * (1 - musd))
    m.addConstr(mod_bar <= musd)
    m.addConstr(psd_bar - psd <= 1 - musd)
    m.addConstr(psd_bar - psd >= -(1 - musd))
    m.addConstr(psd_bar <= musd)
    m.addConstr(moh_bar - moh <= 1 - mush)
    m.addConstr(moh_bar - moh >= -(1 - mush))
    m.addConstr(moh_bar <= mush)
    m.addConstr(hsh_bar - hsh <= 1 - mush)
    m.addConstr(hsh_bar - hsh >= -(1 - mush))
    m.addConstr(hsh_bar <= mush)

    heatrd = m.addVar()
    elerd = m.addVar()
    Qelerd = m.addVar()

    lcaes_e = .5
    lcaes_h = lcaes_e/2

    m.addConstr(lcaes_h*hsh_bar + .1*gsum(h_chp[i] for i in range(U_chp)) == lH[t] - heatrd)
    m.addConstr(lcaes_e*(psd_bar - psc_bar) + .1*gsum(P_chp[i] for i in range(U_chp)) == lE[t] - elerd)
    m.addConstr(.1*gsum(Q_chp[i] for i in range(U_chp)) == 1/4*lE[t] - Qelerd)

    theObj = Png*gsum(ng_chp[u] for u in range(U_chp)) + 6*Pe[t]*(heatrd/2.4+elerd+Qelerd) + gsum(cst_chp[u] + csh_chp[u] for u in range(U_chp))
    if t == nss-1: # penalty for not-satisfying daily balance
        endPen1,endPen2 = m.addVar(),m.addVar()
        m.addConstr(endPen1 >= (Mamin + Mamax) / 2 - (Mamin + MaSCL*(gdot(fbca,Ma))))
        m.addConstr(endPen2 >= (Mhmin + Mhmax) / 2 - (Mhmin + MhSCL*(gdot(fbca,Mh))))
        theObj += endPen1+endPen2

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            m.addConstr(tha >= cut[its,t,-1]+gsum(cut[its,t,i]*Ma[i] for i in range(Bca))+gsum(cut[its,t,Bca+i]*Mh[i] for i in range(Bca))+gsum(cut[its,t,2*Bca+i]*mu_chp[i] for i in range(U_chp)))

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
        if ite == 13:
            print('t = ',t)
            print('purchase P_h_Q = |%8g|%8g|%8g|' % (elerd.X,heatrd.X,Qelerd.X))
            for u in range(U_chp):
                print('chp%1d: P_h_Q = |%8g|%8g|%8g|' % (u,.1*P_chp[u].X,.1*h_chp[u].X,.1*Q_chp[u].X))
            print('caes c/d/h = |%8g|%8g|%8g|' % (lcaes_e*psc_bar.X,lcaes_e * psd_bar.X,lcaes_h * hsh_bar.X))
            print('load: P_h_Q = |%8g|%8g|%8g|' % (lE[t],lH[t],.25*lE[t]))
            if t == nss-1:
                print('penalty = |%8g|%8g|' % (endPen1.X,endPen2.X))
        Mathis,Mhthis = np.zeros(Bca,dtype=np.int8),np.zeros(Bca,dtype=np.int8)
        muchpthis = np.zeros(U_chp,dtype=np.int8)
        for i in range(Bca):
            Mathis[i], Mhthis[i] = Ma[i].X, Mh[i].X
        for i in range(U_chp):
            muchpthis[i] = mu_chp[i].X
        return m.ObjVal-tha.X,np.concatenate((Mathis,Mhthis,muchpthis)),m.ObjVal
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
            tailVal = np.dot(pai,al)
            maxerVal = subVal + tailVal
            cpv = np.zeros(zl)
            for i in range(zl):
                cpv[i] = z[i].X
            drc = al - cpv
            return maxerVal,subVal,tailVal,drc,cpv

fbca,Bca,SCLca = fbgen(.9,64.5)

Ta = 25
vP, fP, pP = 308.9, 626.8, 1044 # valley,flat,peak price of ele (yuan/MWh)
Pe = np.array([vP,pP,fP,pP],dtype=np.float64)
vle, fle, ple = .4, 1.0, 1.6
lE = np.array([vle,vle,fle,ple]) # load is in MW
lH = 3/5*lE
Png = 1.5 # price of NG (yuan/(kg/h * h))

U_chp = 4
stPchp, shPchp = 1044, 1044  # start and shut price of chp
beta_ng2e, eta_pchp, eta_hchp = .015, .45, .5  # units: MWh/kg,1,1
pA_chp, pB_chp, pC_chp, pD_chp, hA_chp, hB_chp, hC_chp, hD_chp = 4, 2.5, .8, 1.7, 0, 1.2, .5, 0  # chp region

Mamin,Mamax,Mhmin,Mhmax = 56053,94333,4145,37308 # dual soc
MaSCL,MhSCL = (Mamax-Mamin)/127,(Mhmax-Mhmin)/127
macB,mocB,madB,modB,mohB = 1.47/3.3,0.6385/3.3,2.39/5.5,1.0372/5.5,0.9075/2.4 # base val

learnRate = 5e6
nss = 4
MaxIte = 50000
cut = np.zeros((MaxIte,nss-1,2*Bca+U_chp+1)) # last stage do not need a cut
act = np.zeros((nss,2*Bca+U_chp))
val,vlb = np.zeros(nss),np.zeros(nss)
lbite = np.zeros(MaxIte)
ubite = 1000000*np.ones(MaxIte)
bestIte = 0
bestub=60000
for ite in range(MaxIte):
    dR = 0
    for t in range(nss):
        val[t],act[t],vlb[t] = cacase(dR,ite,t,act[t-1],0,0)
        if ite == 13:
            print('t(%1d)|%8g|%8g|'%(t,np.dot(act[t,:Bca],fbca),np.dot(act[t,Bca:(2*Bca)],fbca)),end='')
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
    else:
        lbite[ite] = lb

    for i in range(nss-1):
        print('%8g + ' % (val[i]),end='')
    print('%8g' % (val[nss-1]))
    # ----------- statistic -----------
    dR = 1
    for t in range(nss-1,0,-1):
        cutGen(dR, ite, t, act[t-1])

print(bestIte)
print('the ubite:')
for i in ubite[:ite+1]:
    print(i)
print('the lbite:')
for i in lbite[:ite+1]:
    print(i)
