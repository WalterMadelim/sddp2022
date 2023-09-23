import copy

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import matplotlib.pyplot as plt
# lag cut _ but local maxima with gaps
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

fb,B,SCL = fbgen(.9,64.5)
Ta = 25
vP, fP, pP = 308.9, 626.8, 1044 # yuan/MWh
Pe = np.array([vP,pP,fP,pP],dtype=np.float64)
vle, fle, ple = .4, 1.0, 1.6
lE = np.array([vle,vle,fle,ple])
lH = 3/5*lE
Png = 1.5

U_chp = 4
stPchp, shPchp = 1044, 1044  # start and shut price of chp
beta_ng2e, eta_pchp, eta_hchp = .015, .45, .5  # units: MWh/kg,1,1
pA_chp, pB_chp, pC_chp, pD_chp, hA_chp, hB_chp, hC_chp, hD_chp = 4, 2.5, .8, 1.7, 0, 1.2, .5, 0  # chp region

Mamin,Mamax,Mhmin,Mhmax = 56053,94333,4145,37308 # dual soc
MaSCL,MhSCL = (Mamax-Mamin)/127,(Mhmax-Mhmin)/127
macB,mocB,madB,modB,mohB = 1.47/3.3,0.6385/3.3,2.39/5.5,1.0372/5.5,0.9075/2.4 # base val

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

vW, fW, pW = 200,1000,1600 # 200, 1000, 2190
lW = np.array([vW,pW,fW,pW])

tkCoef = 6/160
pipeNum = 6
pumpNum = 2
hrdNum = 3

nss = 4
MaxIte = 50000
# 2*B (CAES_a_h) + B (waterTank) + U_chp
cut = np.zeros((MaxIte,nss-1,3*B+U_chp+1)) # last stage do not need a cut
act = np.zeros((nss,3*B+U_chp))
val,vlb = np.zeros(nss),np.zeros(nss)

def cacase(dR,ite,t,al,s,pai):
    m = gp.Model('sMallcAse')
    if t != 0:
        zl = len(al)
        z = m.addVars(zl, ub=1)
        if s != 2:
            m.addConstrs(z[i] == al[i] for i in range(zl))
    # -----model part-----
    vtype = GRB.CONTINUOUS if s == 1 else GRB.BINARY
    mu_chp = m.addVars(U_chp,vtype=vtype,ub=1)  # 'on' state indicator of chp
    P_chp, Q_chp, h_chp, cst_chp, csh_chp, ng_chp = m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp),m.addVars(U_chp)
    for u in range(U_chp):
        m.addConstr(h_chp[u] <= hB_chp * mu_chp[u])
        m.addConstr(P_chp[u] <= pA_chp * mu_chp[u])
        m.addConstr(Q_chp[u] <= .5 * 3 ** .5 * pA_chp * mu_chp[u])  # m.addConstrs(np.cos((2*k-1)/6*np.pi) * P_chp[t,u] + np.sin((2*k-1)/6*np.pi) * Q_chp[t,u] - np.cos(np.pi/6) * S_chp[t,u] <= 0 for k in range(6))
        m.addConstr(Q_chp[u] + 3 ** 0.5 * (P_chp[u] - pA_chp) <= 0)  # with ub of Q
        m.addConstr(P_chp[u] - pA_chp - (pA_chp - pB_chp) / (hA_chp - hB_chp) * (h_chp[u] - hA_chp) <= 0)
        m.addConstr(P_chp[u] - pD_chp - (pD_chp - pC_chp) / (hD_chp - hC_chp) * (h_chp[u] - hD_chp) >= -1.8 * (1 - mu_chp[u]))
        m.addConstr(P_chp[u] - pB_chp - (pC_chp - pB_chp) / (hC_chp - hB_chp) * (h_chp[u] - hB_chp) >= 0)
    # ------------------------------
    Ma, Mh = m.addVars(B, vtype=vtype,ub=1),m.addVars(B, vtype=vtype,ub=1)
    mac, moc, mad, mod = m.addVar(lb=.9, ub=1), m.addVar(lb=.6, ub=1.3), m.addVar(lb=.7, ub=1), m.addVar(lb=.7, ub=1.3)
    moh, psc, psd, hsh = m.addVar(lb=.8707, ub=1), m.addVar(lb=.8, ub=1), m.addVar(lb=.4, ub=1), m.addVar(lb=.6, ub=1)
    mac_bar, moc_bar, mad_bar, mod_bar, moh_bar, psc_bar, psd_bar, hsh_bar = m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar(), m.addVar()
    musc, musd, mush = m.addVar(vtype=vtype,ub=1),m.addVar(vtype=vtype,ub=1),m.addVar(vtype=vtype,ub=1)
    if t == 0:
        m.addConstr(gdot(fb,Ma) - 64 - 6*3600/MaSCL*(macB * mac_bar - madB * mad_bar) == [-1/2,1/2])
        m.addConstr(gdot(fb,Mh) - 64 - 6 * 3600/MhSCL * (2 * mocB * moc_bar - 2 * modB * mod_bar - mohB * moh_bar) == [-1/2,1/2])
    else:
        m.addConstr(gdot(fb,Ma) - gsum(fb[i]*z[i] for i in range(B)) - 6*3600/MaSCL*(macB * mac_bar - madB * mad_bar) == [-1/2,1/2])
        m.addConstr(gdot(fb,Mh) - gsum(fb[i]*z[B+i] for i in range(B)) - 6 * 3600/MhSCL * (2 * mocB * moc_bar - 2 * modB * mod_bar - mohB * moh_bar) == [-1/2,1/2])
    m.addConstr(psc == fpsc(mac, Ta))
    m.addConstr(mac == fmac(moc, Ta))
    m.addConstr(hsh == fhsh(moh))
    xt_madlb = np.linspace(.7, 1.3, 5)
    yt_madlb = fmadlb(xt_madlb)
    l_madlb, f_madlb = m.addVar(ub=1), m.addVars(4, vtype=vtype,ub=1)
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
    fx_psd, fy_psd = m.addVar(vtype=vtype,ub=1), m.addVar(vtype=vtype,ub=1)
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
    # ------------------------------
    s_pump = m.addVars(pumpNum, 4, vtype=vtype, ub=1)  # represent 4 finite states
    q_pump, h_pump, p_pump = m.addVars(pumpNum, lb=100, ub=1100), m.addVars(pumpNum, lb=86, ub=229), m.addVars(pumpNum,lb=160,ub=556)
    l_pump = m.addVars(pumpNum, ub=1)

    s_pipe = m.addVars(pipeNum, 3, vtype=vtype, ub=1)
    q_pipe, h_pipe = m.addVars(pipeNum, ub=1500), m.addVars(pipeNum, ub=704)
    l0_pipe = m.addVars(pipeNum, ub=1)
    l1_pipe = m.addVars(pipeNum)
    hrd = m.addVars(hrdNum)
    htin = m.addVar()
    htou = m.addVars(B, vtype=vtype, ub=1)
    m.addConstr(gdot(fb, htou) == (gsum(fb[i]*z[2*B+i] for i in range(B)) if t > 0 else 35) + tkCoef * (q_pipe[2] - q_pipe[3]))
    m.addConstr(gdot(fb, htou) <= 70)
    for n in range(pumpNum):
        m.addConstr(gsum(s_pump[n, i] for i in range(4)) == 1)
        m.addConstr(
            l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000 * (1 - s_pump[n, 0]))
        m.addConstr(
            l_pump[n] * qt_pump[0, 0] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000 * (1 - s_pump[n, 0]))
        m.addConstr(
            l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] >= -1000 * (1 - s_pump[n, 1]))
        m.addConstr(
            l_pump[n] * qt_pump[0, 2] + (1 - l_pump[n]) * qt_pump[0, 1] - q_pump[n] <= 1000 * (1 - s_pump[n, 1]))
        m.addConstr(
            l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000 * (1 - s_pump[n, 2]))
        m.addConstr(
            l_pump[n] * qt_pump[1, 0] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000 * (1 - s_pump[n, 2]))
        m.addConstr(
            l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] >= -1000 * (1 - s_pump[n, 3]))
        m.addConstr(
            l_pump[n] * qt_pump[1, 2] + (1 - l_pump[n]) * qt_pump[1, 1] - q_pump[n] <= 1000 * (1 - s_pump[n, 3]))
        m.addConstr(
            l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143 * (1 - s_pump[n, 0]))
        m.addConstr(l_pump[n] * ht_pump[0, 0] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143 * (1 - s_pump[n, 0]))
        m.addConstr(
            l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] >= -143 * (1 - s_pump[n, 1]))
        m.addConstr(l_pump[n] * ht_pump[0, 2] + (1 - l_pump[n]) * ht_pump[0, 1] - h_pump[n] <= 143 * (1 - s_pump[n, 1]))
        m.addConstr(
            l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143 * (1 - s_pump[n, 2]))
        m.addConstr(l_pump[n] * ht_pump[1, 0] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143 * (1 - s_pump[n, 2]))
        m.addConstr(
            l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] >= -143 * (1 - s_pump[n, 3]))
        m.addConstr(l_pump[n] * ht_pump[1, 2] + (1 - l_pump[n]) * ht_pump[1, 1] - h_pump[n] <= 143 * (1 - s_pump[n, 3]))
        m.addConstr(
            l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396 * (1 - s_pump[n, 0]))
        m.addConstr(l_pump[n] * pt_pump[0, 0] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396 * (1 - s_pump[n, 0]))
        m.addConstr(
            l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] >= -396 * (1 - s_pump[n, 1]))
        m.addConstr(l_pump[n] * pt_pump[0, 2] + (1 - l_pump[n]) * pt_pump[0, 1] - p_pump[n] <= 396 * (1 - s_pump[n, 1]))
        m.addConstr(
            l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396 * (1 - s_pump[n, 2]))
        m.addConstr(l_pump[n] * pt_pump[1, 0] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396 * (1 - s_pump[n, 2]))
        m.addConstr(
            l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] >= -396 * (1 - s_pump[n, 3]))
        m.addConstr(l_pump[n] * pt_pump[1, 2] + (1 - l_pump[n]) * pt_pump[1, 1] - p_pump[n] <= 396 * (1 - s_pump[n, 3]))
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
    m.addConstr(q_pump[0] == q_pipe[2] + q_pipe[0])  # node 0
    m.addConstr(q_pipe[0] + q_pump[1] == q_pipe[1])  # node 1
    m.addConstr(q_pipe[1] == q_pipe[5])  # node 2
    m.addConstr(q_pipe[4] + q_pipe[5] == lW[t])
    m.addConstr(q_pipe[4] == q_pipe[3])  # node 4
    # pressure constrs
    m.addConstr(h_pump[0] - h_pipe[2] == htin)  # node tank
    m.addConstr(gdot(fb, htou) - h_pipe[3] - h_pipe[4] - hrd[2] == 0)  # node 6
    m.addConstr(h_pump[0] - h_pipe[0] - hrd[0] == h_pump[1])
    m.addConstr(h_pump[1] - h_pipe[1] - h_pipe[5] - hrd[1] == 0)

    activePM = 1e-3*gsum(p_pump[p] for p in range(pumpNum))

    heatrd = m.addVar()
    elerd = m.addVar()
    Qelerd = m.addVar()

    lcaes_e = .5
    lcaes_h = lcaes_e/2

    m.addConstr(lcaes_h*hsh_bar + .1*gsum(h_chp[i] for i in range(U_chp)) == lH[t] - heatrd)
    m.addConstr(lcaes_e*(psd_bar - psc_bar) + .1*gsum(P_chp[i] for i in range(U_chp)) == activePM + lE[t] - elerd)
    m.addConstr(.1*gsum(Q_chp[i] for i in range(U_chp)) == activePM*np.tan(np.arccos(.85)) + 1/4*lE[t] - Qelerd)

    for u in range(U_chp):
        m.addConstr(beta_ng2e * ng_chp[u] == 6*(.1*P_chp[u] / eta_pchp + .1*h_chp[u] / eta_hchp))  # how much ng consumed by chp
        m.addConstr(csh_chp[u] >= ((shPchp*(1 - mu_chp[u])) if t==0 else (shPchp*(z[3*B+u] - mu_chp[u]))))  # only at t = 0, shut down cost of chp
        m.addConstr(cst_chp[u] >= (0 if t==0 else (stPchp*(mu_chp[u] - z[3*B+u]))))

    theObj = Png*gsum(ng_chp[u] for u in range(U_chp)) + gsum(cst_chp[u] + csh_chp[u] for u in range(U_chp)) + 6*Pe[t]*(heatrd/2.4+elerd+Qelerd)
    if t == nss-1: # penalty for not-satisfying daily balance
        endPen1,endPen2,endPen3 = m.addVar(),m.addVar(),m.addVar()
        m.addConstr(endPen1 >= 100*(64-gdot(fb,Ma)))
        m.addConstr(endPen2 >= 100*(64-gdot(fb,Mh)))
        m.addConstr(endPen3 >= 100*(35-gdot(fb,htou)))
        theObj += endPen1+endPen2+MaSCL*endPen3

    tha = m.addVar()
    if t < nss-1:
        for its in range(ite+dR):
            tmp = cut[its,t]
            m.addConstr(tha >= tmp[-1]+gdot(tmp[:B],Ma)+gdot(tmp[B:2*B],Mh)+gdot(tmp[2*B:3*B],htou)+gdot(tmp[3*B:(3*B+U_chp)],mu_chp))
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
        Mathis,Mhthis,htouthis = np.zeros(B,dtype=np.int8),np.zeros(B,dtype=np.int8),np.zeros(B,dtype=np.int8)
        muchpthis = np.zeros(U_chp,dtype=np.int8)
        for i in range(B):
            Mathis[i], Mhthis[i], htouthis[i] = Ma[i].X, Mh[i].X, htou[i].X
        for i in range(U_chp):
            muchpthis[i] = mu_chp[i].X
        # print('%.6g|%.6g|%.6g' % (psc_bar.X, psd_bar.X, hsh_bar.X))
        # for u in range(U_chp):
        #     print('%.6g|%.6g|%.6g' % (P_chp[u].X, Q_chp[u].X, h_chp[u].X))
        # print('pump,t=%8d,p0=(%8g|%8g|%8g),p1=(%8g|%8g|%8g)'%(t,q_pump[0].X,h_pump[0].X,p_pump[0].X,q_pump[1].X,h_pump[1].X,p_pump[1].X))
        # print('pump,state=(%d|%d|%d|%d|__|%d|%d|%d|%d)'%(s_pump[0,0].X,s_pump[0,1].X,s_pump[0,2].X,s_pump[0,3].X,s_pump[1,0].X,s_pump[1,1].X,s_pump[1,2].X,s_pump[1,3].X))
        print('e_Q_H:',elerd.X, Qelerd.X, heatrd.X)
        return m.ObjVal-tha.X,np.concatenate((Mathis,Mhthis,htouthis,muchpthis)),m.ObjVal
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

bestub = 5e4
train = 0
if train:
    for ite in range(MaxIte):
        dR = 0
        for t in range(nss):
            val[t],act[t],vlb[t] = cacase(dR,ite,t,act[t-1],0,0)
        # ----------- statistic -----------
        lb = vlb[0]
        ub = sum(val)
        if bestub > ub:
            bestIte = ite
            bestub = ub
            bestact = copy.deepcopy(act)
        print('ite=%8d\lb=%12g<%12g<%12g='%(ite,lb,bestub,ub),end='')
        if lb + 1 >= bestub: # final convergent
            print('converge at ite = %8d' % ite)
            print('best at bestIte = %8d' % bestIte)
            np.save('ite.npy', np.array([bestIte]))
            np.save('cut.npy',cut)
            np.save('act.npy',bestact)
            break
        elif (ite+1) % 60 == 0: # manual input ctrl
            input()
            print('save param temporarily at ite = %8d' % ite)
            print('best at bestIte = %8d' % bestIte)
            np.save('ite.npy', np.array([bestIte]))
            np.save('cuttmp.npy', cut)
            np.save('acttmp.npy', bestact)
        for i in range(nss-1):
            print('%8g + ' % (val[i]),end='')
        print('%8g' % (val[nss-1]))
        # ----------- statistic -----------
        dR = 1
        for t in range(nss-1,0,-1):
            cutGen(dR, ite, t, act[t-1])
# test after convergent training
print('This is testing:')
# ite = np.load('ite.npy')[0]
cut = np.load('cut.npy') # training result
ite = 506
dR = 0
for t in range(nss):
    val[t], act[t], vlb[t] = cacase(dR, ite, t, act[t - 1], 0, 0)
print(val)
#
# print('total opt cost',sum(val))
# print('CAES, htou')
#
# for t in range(nss):
#     tmp = act[t]
#     print('t=',t)
#     print(np.dot(fb, tmp[:B])/127)
#     print(np.dot(fb, tmp[B:2*B])/127)
#     print(np.dot(fb, tmp[2*B:3*B])/70)
