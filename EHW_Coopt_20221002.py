import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from matplotlib import pyplot as plt
def getbran32PQ(P33,Q33): # node load 2 branch load
    Pb32, Qb32 = np.zeros_like(x32), np.zeros_like(x32)  # branch active power flow
    Pb32[5:17], Qb32[5:17] = P33[17], Q33[17]
    Pb32[24:32], Qb32[24:32] = P33[32], Q33[32]
    Pb32[2:5], Qb32[2:5] = P33[17] + P33[32], Q33[17] + Q33[32]
    Pb32[21:24], Qb32[21:24] = P33[24], Q33[24]
    Pb32[1], Qb32[1] = P33[17] + P33[32] + P33[24], Q33[17] + Q33[32] + Q33[24]
    Pb32[17:21], Qb32[17:21] = P33[21], Q33[21]
    Pb32[0], Qb32[0] = P33[17] + P33[32] + P33[24] + P33[21], Q33[17] + Q33[32] + Q33[24] + Q33[21]
    return Pb32,Qb32
def fmac(moc,Tam):
    return 0.322 + 0.5621 * moc -2.209e-5 * Tam ** 2 + 0.005004 * Tam
def fpsc(mac,Tam):
    return 0.1672 + 0.8257 * mac - 1.006e-5 * Tam ** 2 + 0.0005409 * Tam
def fmadlb(mod): # do 1-dim pwl
    return 1.9094 * mod ** 4 -8.5877 * mod ** 3 + 14.5357 * mod ** 2 -10.9975 * mod + 4
def fpsd(mad,mod): # do 2-dim pwl
    return -8.05+21.2095*mad+1.7761*mod -20.6123*mad**2 +3.3126*mad*mod -2.754*mod**2 +6.8186*mad**3 -1.0925*mod*mad**2 -0.4445*mad*mod**2 +0.8743*mod**3
def fhsh(moh):
    return 0.31265*moh+0.7 # fixed functions of CAES

T = 8 # scheduling periods
Tam = 25 # ambient temperature
vP, fP, pP = 308.9, 626.8, 1044 # valley,flat,peak price of ele (yuan/MWh)
Pe = np.array([vP,vP,vP,vP,vP,vP,vP,vP,pP,pP,pP,pP,fP,fP,fP,fP,fP,pP,pP,pP,pP,fP,fP,fP])
Png = 1.5 # price of NG (yuan/(kg/h * h))

m = gp.Model('co-gen')
# CHP (multi-units) block: physical-(p_chp,q_chp,h_chp) cost-(cst_chp,csh_chp,ng_chp), base = 1MW
U_chp = 8
stPchp,shPchp = 1044,1044 # start and shut price of chp
beta_ng2e,eta_pchp,eta_hchp = .015,.45,.5 # units: MWh/kg,1,1
pA_chp,pB_chp,pC_chp,pD_chp,hA_chp,hB_chp,hC_chp,hD_chp = 4,2.5,.8,1.7,0,1.2,.5,0 # chp region
mu_chp = m.addVars(T,U_chp,vtype=GRB.BINARY) # 'on' state indicator of chp
P_chp,Q_chp,h_chp,cst_chp,csh_chp,ng_chp = m.addVars(T,U_chp),m.addVars(T,U_chp),m.addVars(T,U_chp),m.addVars(T,U_chp),m.addVars(T,U_chp),m.addVars(T,U_chp)
for u in range(U_chp):# chp assumed to be 'on' initially, so no start up cost at t = 0
    m.addConstr(csh_chp[0,u] >= shPchp*(1-mu_chp[0,u])) # only at t = 0, shut down cost of chp
    for t in range(1,T):
        m.addConstr(cst_chp[t,u] >= stPchp*(mu_chp[t,u]-mu_chp[t-1,u]))
        m.addConstr(csh_chp[t,u] >= shPchp*(mu_chp[t-1,u]-mu_chp[t,u]))
    for t in range(T): # unit of ng_chp: (kg/h), unit of P_chp: (MW)
        m.addConstr(beta_ng2e * ng_chp[t,u] == P_chp[t,u] / eta_pchp + h_chp[t,u] / eta_hchp) # how much ng consumed by chp
        m.addConstr(h_chp[t,u] <= hB_chp * mu_chp[t,u])
        m.addConstr(P_chp[t,u] <= pA_chp * mu_chp[t,u])
        m.addConstr(Q_chp[t,u] <= .5*3**.5*pA_chp * mu_chp[t,u]) # m.addConstrs(np.cos((2*k-1)/6*np.pi) * P_chp[t,u] + np.sin((2*k-1)/6*np.pi) * Q_chp[t,u] - np.cos(np.pi/6) * S_chp[t,u] <= 0 for k in range(6))
        m.addConstr(Q_chp[t,u]+3**0.5*(P_chp[t,u]-pA_chp) <= 0) # with ub of Q
        m.addConstr(P_chp[t,u] - pA_chp - (pA_chp - pB_chp) / (hA_chp - hB_chp) * (h_chp[t,u] - hA_chp) <= 0)
        m.addConstr(P_chp[t,u] - pD_chp - (pD_chp - pC_chp) / (hD_chp - hC_chp) * (h_chp[t,u] - hD_chp) >= -5e3 * (1 - mu_chp[t,u]))
        m.addConstr(P_chp[t,u] - pB_chp - (pC_chp - pB_chp) / (hC_chp - hB_chp) * (h_chp[t,u] - hB_chp) >= 0) # -5e3 * (1 - mu_chp[t,u]))
# CAES
Mamin,Mamax,Mhmin,Mhmax = 56053,94333,4145,37308 # dual soc
macB,mocB,madB,modB,mohB = 1.47,0.6385,2.39,1.0372,0.9075 # base val
mac,moc,mad,mod = m.addVars(T,lb=.9,ub=1),m.addVars(T,lb=.6,ub=1.3),m.addVars(T,lb=.7,ub=1),m.addVars(T,lb=.7,ub=1.3)
moh,psc,psd,hsh = m.addVars(T,lb=.8707,ub=1),m.addVars(T,lb=.8,ub=1),m.addVars(T,lb=.4,ub=1),m.addVars(T,lb=.6,ub=1)
xlb_madlb,xub_madlb,segs_madlb = .7,1.3,4 # 1d-pwl of fmadlb(mod)
xt_madlb = np.linspace(xlb_madlb, xub_madlb, segs_madlb + 1)
yt_madlb = fmadlb(xt_madlb)
l_madlb,f_madlb = m.addVars(T,ub=1),m.addVars(T,segs_madlb,vtype=GRB.BINARY)
m.addConstrs(gp.quicksum( f_madlb[s,r] for r in range(segs_madlb) ) == 1 for s in range(T))
for s in range(T):
    for r in range(segs_madlb):
        m.addConstr((f_madlb[s,r] == 1) >> (mod[s] == l_madlb[s] * xt_madlb[r] + (1 - l_madlb[s]) * xt_madlb[r + 1]))
        m.addConstr((f_madlb[s,r] == 1) >> (mad[s] >= l_madlb[s] * yt_madlb[r] + (1 - l_madlb[s]) * yt_madlb[r + 1]))

xlb_psd,xub_psd,ylb_psd,yub_psd,segs_psd = .7,1,.7,1.3,4 # 2d-pwl of fpsd(mad,mod)
xt_psd,yt_psd = np.linspace(xlb_psd,xub_psd,segs_psd+1),np.linspace(ylb_psd,yub_psd,segs_psd+1)
tmp1,tmp2 = np.meshgrid(xt_psd,yt_psd) # this is an intermediate assign (temporary)
zt_psd = fpsd(tmp1,tmp2)
l1_psd,l2_psd,l3_psd,l4_psd = m.addVars(T,ub=1),m.addVars(T,ub=1),m.addVars(T,ub=1),m.addVars(T,ub=1)
fx_psd,fy_psd,f_psd = m.addVars(T,segs_psd,vtype=GRB.BINARY),m.addVars(T,segs_psd,vtype=GRB.BINARY),m.addVars(T,segs_psd**2,vtype=GRB.BINARY)
m.addConstrs(l1_psd[s]+l2_psd[s]+l3_psd[s]+l4_psd[s]==1 for s in range(T))
m.addConstrs(gp.quicksum( fx_psd[s,r] for r in range(segs_psd) ) == 1 for s in range(T))
m.addConstrs(gp.quicksum( fy_psd[s,r] for r in range(segs_psd) ) == 1 for s in range(T))
for s in range(T):
    for r in range(segs_psd):
        m.addConstr((fx_psd[s,r]==1)>>(mad[s] == (l1_psd[s]+l3_psd[s])*xt_psd[r]+(l2_psd[s]+l4_psd[s])*xt_psd[r+1]))
        m.addConstr((fy_psd[s,r]==1)>>(mod[s] == (l1_psd[s]+l2_psd[s])*yt_psd[r]+(l3_psd[s]+l4_psd[s])*yt_psd[r+1]))
for s in range(T):
    cnt = 0
    for rx in range(segs_psd):
        for ry in range(segs_psd):
            m.addGenConstrAnd(f_psd[s,cnt],[fx_psd[s,rx],fy_psd[s,ry]])
            m.addConstr((f_psd[s,cnt]==1)>>(psd[s]==l1_psd[s]*zt_psd[ry,rx] + l2_psd[s]*zt_psd[ry,rx+1] + l3_psd[s]*zt_psd[ry+1,rx] + l4_psd[s]*zt_psd[ry+1,rx+1]))
            cnt += 1 # pwl part

mac_bar,moc_bar,mad_bar,mod_bar,moh_bar,psc_bar,psd_bar,hsh_bar = m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T),m.addVars(T)
musc,musd,mush = m.addVars(T,vtype=GRB.BINARY),m.addVars(T,vtype=GRB.BINARY),m.addVars(T,vtype=GRB.BINARY)
for s in range(T): # CAES dispatch
    m.addConstr(musc[s] + musd[s] <= 1)
    m.addConstr(mac_bar[s] == musc[s] * mac[s])
    m.addConstr(moc_bar[s] == musc[s] * moc[s])
    m.addConstr(psc_bar[s] == musc[s] * psc[s])
    m.addConstr(mad_bar[s] == musd[s] * mad[s])
    m.addConstr(mod_bar[s] == musd[s] * mod[s])
    m.addConstr(psd_bar[s] == musd[s] * psd[s])
    m.addConstr(hsh_bar[s] == mush[s] * hsh[s])
    m.addConstr(moh_bar[s] == mush[s] * moh[s])
    m.addConstr(psc[s] == fpsc(mac[s],Tam))
    m.addConstr(mac[s] == fmac(moc[s], Tam))
    m.addConstr(hsh[s] == fhsh(moh[s]))
# storage management
Ma,Mh = m.addVars(T,lb=Mamin,ub=Mamax),m.addVars(T,lb=Mhmin,ub=Mhmax)
m.addConstr(Ma[0] == (Mamin+Mamax)/2 + 3600*(macB*mac_bar[0]-madB*mad_bar[0]))
m.addConstrs(Ma[s] == Ma[s-1] + 3600*(macB*mac_bar[s]-madB*mad_bar[s]) for s in range(1,T))
# m.addConstr(Ma[T-1] == (Mamin+Mamax)/2)
m.addConstr(Mh[0] == (Mhmin+Mhmax)/2 + 3600*(2*mocB*moc_bar[0]-2*modB*mod_bar[0]-mohB*moh_bar[0]))
m.addConstrs(Mh[s] == Mh[s-1] + 3600*(2*mocB*moc_bar[s]-2*modB*mod_bar[s]-mohB*moh_bar[s]) for s in range(1,T))
# m.addConstr(Mh[T-1] == (Mhmin+Mhmax)/2) # CAES: psc_bar[s], psd_bar[s], hsh_bar[s]

# H-Network CHP,heat pump with TES at node 0.l3,l4,l5: heat load represented as temperature difference (since mdot is fixed)
l3 = 1/3*np.array([67.062569,69.618591,72.980553,63.927284,66.37442,65.317696,56.202152,63.019165,65.380951,75.999237,91.044609,71.134369,56.963692,63.618328,56.504597,70.608322,70.137527,67.183121,49.154568,57.468884,55.974743,61.736359,60.441044,52.715538])
l4 = 1/3*np.array([56.612175,58.35405,61.938713,59.884624,46.24231,52.470699,52.947582,52.702362,51.537457,58.487629,57.184742,52.763432,53.91449,62.615383,61.054836,60.288395,56.089123,51.948872,57.203434,52.288155,64.545334,52.880962,71.769089,58.06044])
l5 = 1/3*np.array([60.14608,48.359856,60.912048,61.684097,53.240303,46.026188,51.088699,48.642979,57.088631,73.323219,76.385078,73.230881,60.699425,74.634453,77.290367,76.362175,78.169487,76.306686,62.259544,59.063469,63.923267,70.608986,68.152222,68.436165])
Cm = 4200
P_cir,Q_cir = 63.5e-6,39.4e-6 # circulation pump consumed power (MW)
hdf = np.exp(-0.12*1000*np.array([3.5,1.75,1.75,.75,1.75])/4200/np.array([49,25,15,10,24])) # heat dissipation factor
Ts,Tr = m.addVars(T,6),m.addVars(T,6)
P_hp, H_TESc, H_TESd, E_TES = m.addVars(T), m.addVars(T), m.addVars(T), m.addVars(T) # P_hp in (MW), H_TESc in (MW), E_TES[ini] = 5MWh
m.addConstr(E_TES[0] == .9 * 5 + .9 * H_TESc[0] - H_TESd[0] / .9)
m.addConstrs(E_TES[t] == .9 * E_TES[t-1] + .9 * H_TESc[t] - H_TESd[t] / .9 for t in range(1, T))
for t in range(T):
    m.addConstr(Ts[t,0] == 100) # the supply temperature of CHP unit assumed to be fixed at 100
    m.addConstr(Ts[t,1] == Tam + hdf[0]*(Ts[t,0]-Tam))
    m.addConstr(Ts[t,2] == Tam + hdf[1]*(Ts[t,1]-Tam))
    m.addConstr(Ts[t,3] == Tam + hdf[2]*(Ts[t,2]-Tam))
    m.addConstr(Ts[t,4] == Tam + hdf[4]*(Ts[t,1]-Tam))
    m.addConstr(Ts[t,5] == Tam + hdf[3]*(Ts[t,2]-Tam))
    m.addConstr(Tr[t,0] == Tam + hdf[0]*(Tr[t,1]-Tam))
    m.addConstr(Tr[t,1] == (25 * (Tam + hdf[1] * (Tr[t,2] - Tam)) + 24 * (Tam + hdf[4] * (Tr[t,4] - Tam))) / 49)
    m.addConstr(Tr[t,2] == (15 * (Tam + hdf[2] * (Tr[t,3] - Tam)) + 10 * (Tam + hdf[3] * (Tr[t,5] - Tam))) / 25)
    m.addConstr(Tr[t,3] == Ts[t,3] - l3[t])
    m.addConstr(Tr[t,4] == Ts[t,4] - l4[t])
    m.addConstr(Tr[t,5] == Ts[t,5] - l5[t])
# W-Network
lW = 2*np.array([368,395,395,401,380,358,357,319,314,343,357,335,334,324,316,319,331,341,338,331,398,422,385,352]) # 24h-water load
tkCoef = 3600/20
qt_pump = np.array([[1,7,11],[1,5.5,8.8]])*100/3600 # 3 points: 2 segments
ht_pump = np.array([[229,210,149],[138,122,86]])
pt_pump = np.array([[330,500,556],[160,233,262]]) # kW
qt_pipe = np.array([100,500,1000])/3600
ppK = np.array([1547.2653,1041.4285,833.1428,5207.1425,3644.99975,2395.28555])# h=Kq**2, q in (m3/s), h in (m)
ht_pipe = np.zeros((len(ppK),3),dtype=np.float64)
for i,e in enumerate(ppK):
    ht_pipe[i,:] = e * qt_pipe**2
pipesNum = 6
pumpNum = 2
s_pump = m.addVars(T,pumpNum,4,vtype=GRB.BINARY) # represent 4 finite states
q_pump,h_pump,p_pump = m.addVars(T,pumpNum),m.addVars(T,pumpNum),m.addVars(T,pumpNum)
l_pump = m.addVars(T, pumpNum, ub=1)
for t in range(T):
    for n in range(pumpNum):
        m.addConstr(gp.quicksum(s_pump[t, n, i] for i in range(4)) == 1)
        m.addConstr((s_pump[t, n, 0] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[0, 0] + (1 - l_pump[t, n]) * qt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 0] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[0, 0] + (1 - l_pump[t, n]) * ht_pump[0, 1]))
        m.addConstr((s_pump[t, n, 0] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[0, 0] + (1 - l_pump[t, n]) * pt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[0, 2] + (1 - l_pump[t, n]) * qt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[0, 2] + (1 - l_pump[t, n]) * ht_pump[0, 1]))
        m.addConstr((s_pump[t, n, 1] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[0, 2] + (1 - l_pump[t, n]) * pt_pump[0, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[1, 0] + (1 - l_pump[t, n]) * qt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[1, 0] + (1 - l_pump[t, n]) * ht_pump[1, 1]))
        m.addConstr((s_pump[t, n, 2] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[1, 0] + (1 - l_pump[t, n]) * pt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (q_pump[t, n] == l_pump[t, n] * qt_pump[1, 2] + (1 - l_pump[t, n]) * qt_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (h_pump[t, n] == l_pump[t, n] * ht_pump[1, 2] + (1 - l_pump[t, n]) * ht_pump[1, 1]))
        m.addConstr((s_pump[t, n, 3] == 1) >> (p_pump[t, n] == l_pump[t, n] * pt_pump[1, 2] + (1 - l_pump[t, n]) * pt_pump[1, 1]))
s_pipe = m.addVars(T,pipesNum,vtype=GRB.BINARY) # represent 3 finite states
q_pipe,h_pipe = m.addVars(T,pipesNum),m.addVars(T,pipesNum)
l_pipe = m.addVars(T,pipesNum)
for t in range(T):
    for n in range(pipesNum):
        m.addConstr((s_pipe[t,n]==0)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[0]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n]==0)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,0]+(1-l_pipe[t,n])*ht_pipe[n,1]))
        m.addConstr((s_pipe[t,n]==1)>>(q_pipe[t,n]==l_pipe[t,n]*qt_pipe[2]+(1-l_pipe[t,n])*qt_pipe[1]))
        m.addConstr((s_pipe[t,n]==1)>>(h_pipe[t,n]==l_pipe[t,n]*ht_pipe[n,2]+(1-l_pipe[t,n])*ht_pipe[n,1]))
hrd = m.addVars(T,3) # RDV
h0, h1, h2, h3, h4 = m.addVars(T), m.addVars(T), m.addVars(T), m.addVars(T), m.addVars(T)
htin,htou = m.addVars(T),m.addVars(T)
m.addConstr(htou[0] == 100 + tkCoef * (q_pipe[0, 2] - q_pipe[0, 3]))
m.addConstrs(htou[t] == htou[t-1] + tkCoef * (q_pipe[t, 2] - q_pipe[t, 3]) for t in range(1, T))
# m.addConstr(htou[T-1] >= 80)
for t in range(T):
    # flow constrs
    m.addConstr(q_pump[t,0] == q_pipe[t,2] + q_pipe[t,0]) # node 0
    m.addConstr(q_pipe[t,0]+q_pump[t,1] == q_pipe[t,1]) # node 1
    m.addConstr(q_pipe[t,1] == q_pipe[t,5]) # node 2
    m.addConstr(q_pipe[t,4]+q_pipe[t,5] == lW[t]/3600) # node 3: load level: 140~1000
    m.addConstr(q_pipe[t,4] == q_pipe[t,3]) # node 4
    # pressure constrs
    m.addConstr(h_pump[t,0] == h0[t]) # node 0_with pump
    m.addConstr(h_pump[t,1] == h1[t]) # node 1_with pump
    m.addConstr(h0[t]-hrd[t,0]-h_pipe[t,0] == h1[t])# node 1
    m.addConstr(h1[t]-h_pipe[t,1] == h2[t])# node 2
    m.addConstr(h2[t] - h_pipe[t,5] - hrd[t,1] == h3[t])# node 5
    m.addConstr(h0[t]-h_pipe[t,2] == htin[t])# node tank
    m.addConstr(htou[t] - h_pipe[t,3] == h4[t])# node 6
    m.addConstr(h4[t] - h_pipe[t, 4] - hrd[t, 2] == h3[t])# node 5
# E-Network: l17 are p.u. load
Vbase,Sbase = 10000,1e6
zbase = Vbase ** 2 / Sbase
l17 = 1*np.array([0.412,0.415,0.413,0.494,0.417,0.42,0.421,0.55,0.578,0.704,1.07,1.4,1.17,1.3,1.36,1.29,1.27,1.2,0.432,0.552,0.4,0.461,0.411,0.432]) # p.u.
l21 = 1*np.array([0.494,0.409,0.519,0.405,0.427,0.41,0.405,0.514,0.652,0.786,0.769,0.807,1.19,1.33,1.21,1.4,1.37,1.34,1.16,0.425,0.424,0.417,0.438,0.4])
l24 = 1*np.array([0.422,0.406,0.408,0.425,0.423,0.464,0.441,0.406,0.453,0.866,1.35,1.4,1.12,1.29,1.31,1.31,1.23,1.24,1.2,0.4,0.4,0.418,0.432,0.423])
l32 = 1*np.array([0.42,0.4,0.4,0.4,0.4,0.47,0.42,0.48,0.54,1.14,1.27,1.31,1.1,1.18,1.31,1.4,1.35,1.3,1.12,1.12,1.1,1.11,1.13,0.43])
lft = -1+np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12,13],[13,14],[14,15],[15,16],[16,17],[17,18],[2,19],[19,20],[20,21],[21,22],[3,23],[23,24],[24,25],[6,26],[26,27],[27,28],[28,29],[29,30],[30,31],[31,32],[32,33]])
r32 = 1/zbase*np.array([0.00922,0.0493,0.0366,0.0381,0.0819,0.0187,0.0712,0.103,0.104,0.0197,0.0374,0.147,0.0542,0.0591,0.0746,0.129,0.0732,0.0164,0.15,0.041,0.0709,0.0451,0.0898,0.0896,0.0203,0.0284,0.106,0.0804,0.0507,0.0975,0.0311,0.0341])
x32 = 1/zbase*np.array([0.0047,0.0251,0.0186,0.0194,0.0707,0.0619,0.0235,0.074,0.074,0.00651,0.013,0.115,0.0713,0.0526,0.0545,0.172,0.0574,0.0157,0.136,0.0478,0.0937,0.0308,0.0709,0.0707,0.0103,0.0145,0.0934,0.0701,0.0259,0.0963,0.0362,0.053])
Pb32,Qb32 = np.zeros((T,32),dtype=np.float64),np.zeros((T,32),dtype=np.float64)
for t in range(T): # generate the per-period branch PQ-load
    P33 = np.zeros(33,dtype=np.float64)
    P33[17],P33[21],P33[24],P33[32] = l17[t],l21[t],l24[t],l32[t]
    Q33 = copy.deepcopy(P33)
    Q33 /= 4
    Pb32[t],Qb32[t] = getbran32PQ(P33,Q33) # actually use only PQ at 4 points [17,21,24,32]
V33 = m.addVars(T,33,ub=1.05,lb=0.95) # node voltage per-period
for t in range(T): # calculate node voltage
    m.addConstr(V33[t,0]==1.05) # root node voltage keep const
    m.addConstrs(V33[t,lft[b,1]] == V33[t,lft[b,0]] - (Pb32[t,b]*r32[b] + Qb32[t,b]*x32[b]) for b in range(32)) # already checked

P_pur, Q_pur = m.addVars(T),m.addVars(T) # direct purchase from market
for t in range(T): # np.tan(np.arccos(0.85)) = 0.62
    activePM = 1e-3*gp.quicksum(p_pump[t,p] for p in range(2))
    m.addConstr(gp.quicksum(P_chp[t,u] for u in range(U_chp)) + P_pur[t] + psd_bar[t] == psc_bar[t] + Pb32[t,0] + activePM + P_hp[t] + P_cir)
    m.addConstr(gp.quicksum(Q_chp[t,u] for u in range(U_chp)) + Q_pur[t] == Qb32[t,0] + .62*activePM + Q_cir)
    m.addConstr(H_TESc[t] + 1e-6*Cm*49*(Ts[t,0]-Tr[t,0])==.5*hsh_bar[t] +H_TESd[t]+2.4*P_hp[t]+gp.quicksum(h_chp[t,u] for u in range(U_chp)))

C_CHP = gp.quicksum(cst_chp[t,u]+csh_chp[t,u]+Png*ng_chp[t,u] for u in range(U_chp) for t in range(T)) # already checked
rdvPen = 10*gp.quicksum(hrd[t,i] for i in range(3) for t in range(T))
ebuy = gp.quicksum(Pe[t]*(P_pur[t]+Q_pur[t]) for t in range(T))
m.setObjective(rdvPen+ebuy+C_CHP)

m.optimize()
if m.status != GRB.OPTIMAL:
    print('fail')
    exit(3)
print('\nObjVal:%12g' % m.ObjVal)

for u in range(U_chp):
    print('chp%d:'%u)
    for t in range(T):
        print('chp output:%8g|%8g|%8g'%(P_chp[t,u].X,Q_chp[t,u].X,h_chp[t,u].X))

for t in range(T):
    print('purchase:%8g|%8g'%(P_pur[t].X,Q_pur[t].X))

print('heat pump:')
for t in range(T):
    print(P_hp[t].X)

print('TES in DHN')
for t in range(T):
    print('%8g|%8g|%8g' % (H_TESc[t].X, H_TESd[t].X, E_TES[t].X))

print('TK in WSN')
for t in range(T):
    print('%8g' % (htou[t].X))

print('pump power')
for t in range(T):
    print('%8g | %8g' % (p_pump[t,0].X/1000,p_pump[t,1].X/1000))


print('\nCAES air  SoC:',end='')
for s in range(T):
    print('%8.4g,' % ((Ma[s].X - Mamin)/(Mamax-Mamin)),end='')

print('\nCAES heat Soc:',end='')
for s in range(T):
    print('%8.4g,' % ((Mh[s].X-Mhmin)/(Mhmax-Mhmin)),end='')

print('\nCAES charging:',end='')
for s in range(T):
    print('%8.4g,' % psc_bar[s].X,end='')

print('\nCAES discharg:',end='')
for s in range(T):
    print('%8.4g,' % psd_bar[s].X,end='')

print('\nCAES heat pwr:',end='')
for s in range(T):
    print('%8.4g,' % hsh_bar[s].X,end='')

