import gurobipy as gp
from gurobipy import GRB
import numpy as np
# 4-stage, small-case, 20221011
# model setting:
# pure CAES schedule, with learnRate = 150000
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
def cutGen(dR,ite,t,al):
    global cut
    valt = cacase(dR, ite, t, al, 0, 0)
    pai = cacase(dR, ite, t, al, 1, 0)
    maxerVal, subVal, tailVal, drc, cpv = cacase(dR, ite, t, al, 2, pai)
    cut[ite, t - 1, 0] = np.append(pai, subVal)  # sB cut
    while 1:
        maxerVal, subVal, tailVal, drc, cpv = cacase(dR, ite, t, al, 2, pai)
        # print('\nf(pai)=%8g=%8g+%8g' % (maxerVal,subVal,tailVal))
        # print('drc is:',drc,'\ncpv is:',cpv)
        # print('t=%8d: maxer = %8g|stdVal = %8g' % (t,maxerVal,valt))
        if sum(abs(drc)) < 1 / 1024:  # or maxerVal + 1e-4 >= valt:
            # print('subGrad = 0, opt Attained!')
            cut[ite, t - 1, 1] = np.append(pai, subVal)
            # check tightness
            # print('t=%8d: maxer = %8g|stdVal = %8g' % (t,maxerVal,valt))
            if maxerVal + 1e-4 <= valt:
                print('gap = %8g' % (valt - maxerVal))
            break
        pai += 150000 * drc
def cacase(dR,ite,t,al,s,pai):
    m = gp.Model('water_t')
    if t != 0:
        zl = len(al)  # length of copy vector = length of lastAction
        z = m.addVars(zl, ub=1)
        if s != 2:
            m.addConstrs(z[i] == al[i] for i in range(zl))
    # -----model part-----
    Ma, Mh = m.addVars(Bca, vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1), m.addVars(Bca, vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1) # Mamin + MaSCL * (gdot(fbca, Ma))
    mac, moc, mad, mod = m.addVar(lb=.9, ub=1), m.addVar(lb=.6, ub=1.3), m.addVar(lb=.7, ub=1), m.addVar(lb=.7, ub=1.3)
    moh, psc, psd, hsh = m.addVar(lb=.8707, ub=1), m.addVar(lb=.8, ub=1), m.addVar(lb=.4, ub=1), m.addVar(lb=.6, ub=1)
    mac_bar, moc_bar, mad_bar, mod_bar, moh_bar, psc_bar, psd_bar, hsh_bar = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    musc, musd, mush = m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1),m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1),m.addVar(vtype=GRB.CONTINUOUS if s == 1 else GRB.BINARY,ub=1)

    m.addConstr(( (Mamin + Mamax) / 2 if t == 0 else (Mamin + MaSCL * gsum(fbca[i]*z[i] for i in range(Bca))) ) + 3600 * (macB * mac_bar - madB * mad_bar) - (Mamin + MaSCL * (gdot(fbca, Ma))) == [-MaSCL/2,MaSCL/2])
    m.addConstr(( (Mhmin + Mhmax) / 2 if t == 0 else (Mhmin + MhSCL * gsum(fbca[i]*z[Bca+i] for i in range(Bca))) ) + 3600 * (2 * mocB * moc_bar - 2 * modB * mod_bar - mohB * moh_bar) - (Mhmin + MhSCL * (gdot(fbca, Mh))) == [-MhSCL/2,MhSCL/2])

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

    m.addConstr(.4*hsh_bar + heatrd == .001*lH[t])
    m.addConstr(.5*psd_bar - .5*psc_bar + elerd == .001*lE[t])

    theObj = Pe[t]*(heatrd+elerd)
    # if t == nss-1: # penalty for not-satisfying daily balance
    #     endPen = m.addVar()
    #     m.addConstr(endPen >= 35-gdot(fbhtou,htou))
    #     theObj += 500*endPen

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            for typ in range(2):
                m.addConstr(tha >= cut[its,t,typ,-1]+gsum(cut[its,t,typ,i]*Ma[i] for i in range(Bca))+gsum(cut[its,t,typ,Bca+i]*Mh[i] for i in range(Bca)))

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
        Mathis,Mhthis = np.zeros(Bca),np.zeros(Bca)
        for i in range(Bca):
            Mathis[i], Mhthis[i] = Ma[i].X, Mh[i].X
        return m.ObjVal-tha.X,np.concatenate((Mathis,Mhthis)),m.ObjVal
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
lH = 0.3*Pe # load is in kW
lE = 0.5*Pe

Mamin,Mamax,Mhmin,Mhmax = 56053,94333,4145,37308 # dual soc
MaSCL,MhSCL = (Mamax-Mamin)/127,(Mhmax-Mhmin)/127
macB,mocB,madB,modB,mohB = 1.47,0.6385,2.39,1.0372,0.9075 # base val

nss = 4
MaxIte = 50000
cut = np.zeros((MaxIte,nss-1,2,2*Bca+1)) # last stage do not need a cut
act = np.zeros((nss,2*Bca))
val,vlb = np.zeros(nss),np.zeros(nss)

for ite in range(MaxIte):
    dR = 0
    for t in range(nss):
        val[t],act[t],vlb[t] = cacase(dR,ite,t,act[t-1],0,0)
        print('t(%1d)|%8g|%8g|'%(t,np.dot(act[t,:Bca],fbca),np.dot(act[t,Bca:],fbca)),end='')
    # ----------- statistic -----------
    lb = vlb[0]
    ub = sum(val)
    print('ite=%8d\lb=%12g<%12g='%(ite,lb,ub),end='')
    for i in range(nss-1):
        print('%8g + ' % (val[i]),end='')
    print('%8g' % (val[nss-1]))
    # ----------- statistic -----------
    dR = 1
    for t in range(nss-1,0,-1):
        cutGen(dR, ite, t, act[t-1])
