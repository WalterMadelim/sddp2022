import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import sys
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



Tam = 25
stages = 6
Mamin,Mamax,Mhmin,Mhmax = 56053,94333,4145,37308 # dual soc
macB,mocB,madB,modB,mohB = 1.47,0.6385,2.39,1.0372,0.9075 # base val

ePricev, epricef, epricep = 308.9, 626.8, 1044 # yuan/MWh

NgPrice = 1.5 # yuan/kg
stPricechp,shPricechp = 1044,1044 # start shut price
b_ng2e = .015 # MWh/kg
eta_pchp,eta_hchp = .45,.5
pA_chp,pB_chp,pC_chp,pD_chp,hA_chp,hB_chp,hC_chp,hD_chp = 4,2.5,.8,1.7,0,1.2,.5,0 # chp region
eload,hload,price = [1.8,2,3,1.8,3,4.8],[.8,.9,1.2,1.5,1.2,1.5],[ePricev, ePricev, epricef, epricep, epricef, epricep]
qload = [0.45,0.5,0.75,0.45,0.75,1.2]
bigM = 5000
m = gp.Model('myopt')

muchp = m.addVars(stages,vtype=GRB.BINARY)
Schp,pchp,qchp,hchp = m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages)
cstchp,cshchp = m.addVars(stages),m.addVars(stages) # cost of startup and shutdown, pre-schedule state = on
ngchp = m.addVars(stages)
m.addConstr(cstchp[0] >= stPricechp*(muchp[0]-1))
m.addConstr(cshchp[0] >= shPricechp*(1-muchp[0]))
m.addConstrs(cstchp[s] >= stPricechp*(muchp[s]-muchp[s-1]) for s in range(1,stages))
m.addConstrs(cshchp[s] >= shPricechp*(muchp[s-1]-muchp[s]) for s in range(1,stages))
for s in range(stages):
    m.addConstr(b_ng2e * ngchp[s] == pchp[s] / eta_pchp + hchp[s] / eta_hchp)
    m.addConstr(hchp[s] <= hB_chp * muchp[s])
    m.addConstr(pchp[s] <= pA_chp * muchp[s])
    m.addConstr(Schp[s] <= pA_chp * muchp[s])
    m.addConstrs(np.cos((2*k-1)/6*np.pi) * pchp[s] + np.sin((2*k-1)/6*np.pi) * qchp[s] <= np.cos(np.pi/6) * Schp[s] for k in range(6))
    m.addConstr(pchp[s] - pA_chp - (pA_chp-pB_chp)/(hA_chp-hB_chp)*(hchp[s]-hA_chp) <= 0)
    m.addConstr(pchp[s] - pD_chp - (pD_chp-pC_chp)/(hD_chp-hC_chp)*(hchp[s]-hD_chp) >= -bigM*(1-muchp[s]))
    m.addConstr(pchp[s] - pB_chp - (pC_chp-pB_chp)/(hC_chp-hB_chp)*(hchp[s]-hB_chp) >= -bigM*(1-muchp[s])) # CHP: pchp,qchp,hchp[s]

xlb_madlb,xub_madlb,segs_madlb = .7,1.3,4 # 1d-pwl of fmadlb(mod)
xt_madlb = np.linspace(xlb_madlb, xub_madlb, segs_madlb + 1)
yt_madlb = fmadlb(xt_madlb)
# fig,ax = plt.subplots(1,1,figsize=(2.5,2.4))
# ax.plot(xt_madlb,yt_madlb,color='Black')
# plt.show()
# exit(5)
l_madlb,f_madlb = m.addVars(stages,ub=1),m.addVars(stages,segs_madlb,vtype=GRB.BINARY)
x_madlb,y_madlb = m.addVars(stages,lb=-GRB.INFINITY),m.addVars(stages,lb=-GRB.INFINITY)
m.addConstrs(gp.quicksum( f_madlb[s,r] for r in range(segs_madlb) ) == 1 for s in range(stages))
for s in range(stages):
    for r in range(segs_madlb):
        m.addConstr((f_madlb[s,r] == 1) >> (x_madlb[s] == l_madlb[s] * xt_madlb[r] + (1 - l_madlb[s]) * xt_madlb[r + 1]))
        m.addConstr((f_madlb[s,r] == 1) >> (y_madlb[s] == l_madlb[s] * yt_madlb[r] + (1 - l_madlb[s]) * yt_madlb[r + 1]))

xlb_psd,xub_psd,ylb_psd,yub_psd,segs_psd = .7,1,.7,1.3,4 # 2d-pwl of fpsd(mad,mod)
xt_psd,yt_psd = np.linspace(xlb_psd,xub_psd,segs_psd+1),np.linspace(ylb_psd,yub_psd,segs_psd+1)
x_psd,y_psd = np.meshgrid(xt_psd,yt_psd) # this is an intermediate assign (temporary)
zt_psd = fpsd(x_psd,y_psd)
l1_psd,l2_psd,l3_psd,l4_psd = m.addVars(stages,ub=1),m.addVars(stages,ub=1),m.addVars(stages,ub=1),m.addVars(stages,ub=1)
fx_psd,fy_psd,f_psd = m.addVars(stages,segs_psd,vtype=GRB.BINARY),m.addVars(stages,segs_psd,vtype=GRB.BINARY),m.addVars(stages,segs_psd**2,vtype=GRB.BINARY)
x_psd,y_psd,z_psd = m.addVars(stages,lb=-GRB.INFINITY),m.addVars(stages,lb=-GRB.INFINITY),m.addVars(stages,lb=-GRB.INFINITY)
m.addConstrs(l1_psd[s]+l2_psd[s]+l3_psd[s]+l4_psd[s]==1 for s in range(stages))
m.addConstrs(gp.quicksum( fx_psd[s,r] for r in range(segs_psd) ) == 1 for s in range(stages))
m.addConstrs(gp.quicksum( fy_psd[s,r] for r in range(segs_psd) ) == 1 for s in range(stages))
for s in range(stages):
    for r in range(segs_psd):
        m.addConstr((fx_psd[s,r]==1)>>(x_psd[s] == (l1_psd[s]+l3_psd[s])*xt_psd[r]+(l2_psd[s]+l4_psd[s])*xt_psd[r+1]))
        m.addConstr((fy_psd[s,r]==1)>>(y_psd[s] == (l1_psd[s]+l2_psd[s])*yt_psd[r]+(l3_psd[s]+l4_psd[s])*yt_psd[r+1]))
for s in range(stages):
    cnt = 0
    for rx in range(segs_psd):
        for ry in range(segs_psd):
            m.addGenConstrAnd(f_psd[s,cnt],[fx_psd[s,rx],fy_psd[s,ry]])
            m.addConstr((f_psd[s,cnt]==1)>>(z_psd[s]==l1_psd[s]*zt_psd[ry,rx] + l2_psd[s]*zt_psd[ry,rx+1] + l3_psd[s]*zt_psd[ry+1,rx] + l4_psd[s]*zt_psd[ry+1,rx+1]))
            cnt += 1 # pwl part

mac,moc,mad,mod = m.addVars(stages,lb=.9,ub=1),m.addVars(stages,lb=.6,ub=1.3),m.addVars(stages,lb=.7,ub=1),m.addVars(stages,lb=.7,ub=1.3)
moh,psc,psd,hsh = m.addVars(stages,lb=.8707,ub=1),m.addVars(stages,lb=.8,ub=1),m.addVars(stages,lb=.4,ub=1),m.addVars(stages,lb=.6,ub=1)

mac_bar,moc_bar,mad_bar,mod_bar,moh_bar,psc_bar,psd_bar,hsh_bar = m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages),m.addVars(stages)
musc,musd,mush = m.addVars(stages,vtype=GRB.BINARY),m.addVars(stages,vtype=GRB.BINARY),m.addVars(stages,vtype=GRB.BINARY)
for s in range(stages): # CAES dispatch
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
    m.addConstr(psd[s] == z_psd[s])
    m.addConstr(mad[s] == x_psd[s])
    m.addConstr(mod[s] == y_psd[s])
    m.addConstr(mad[s] >= y_madlb[s])
    m.addConstr(x_madlb[s] == mod[s])
    m.addConstr(hsh[s] == fhsh(moh[s]))
# storage management
Ma,Mh = m.addVars(stages,lb=Mamin,ub=Mamax),m.addVars(stages,lb=Mhmin,ub=Mhmax)
m.addConstr(Ma[0] == (Mamin+Mamax)/2 + 3600*(macB*mac_bar[0]-madB*mad_bar[0]))
m.addConstrs(Ma[s] == Ma[s-1] + 3600*(macB*mac_bar[s]-madB*mad_bar[s]) for s in range(1,stages))
m.addConstr(Ma[stages-1] == (Mamin+Mamax)/2)
m.addConstr(Mh[0] == (Mhmin+Mhmax)/2 + 3600*(2*mocB*moc_bar[0]-2*modB*mod_bar[0]-mohB*moh_bar[0]))
m.addConstrs(Mh[s] == Mh[s-1] + 3600*(2*mocB*moc_bar[s]-2*modB*mod_bar[s]-mohB*moh_bar[s]) for s in range(1,stages))
m.addConstr(Mh[stages-1] == (Mhmin+Mhmax)/2) # CAES: psc_bar[s], psd_bar[s], hsh_bar[s]
ph,hh = m.addVars(stages,ub=2),m.addVars(stages,ub=5)
m.addConstrs(hh[s] == 2.4 * ph[s] for s in range(stages)) # heat pump: hh[s], ph[s]

# grid purchase and bus-balance
pG = m.addVars(stages,ub=5)
pQ = m.addVars(stages,ub=5)
m.addConstrs(hh[s] + hchp[s] + 0.5*hsh_bar[s] == hload[s] for s in range(stages))
m.addConstrs(pG[s] + pchp[s] + psd_bar[s] == psc_bar[s] + ph[s] + eload[s] for s in range(stages)) # hub equ in named value
m.addConstrs(pQ[s] + qchp[s] == qload[s] for s in range(stages))

m.setObjective(gp.quicksum(price[s] * (pG[s] + pQ[s]) + NgPrice * ngchp[s] + cstchp[s] + cshchp[s] for s in range(stages)))
m.optimize()
if m.status != GRB.OPTIMAL:
    print('>> Q(x,xi) eval problem: No solution!')
    sys.exit(3)
print('objval is %g' % m.ObjVal)

print('\nCAES air  SoC:',end='')
for s in range(stages):
    print('%8.4g,' % ((Ma[s].X - Mamin)/(Mamax-Mamin)),end='')

print('\nCAES heat Soc:',end='')
for s in range(stages):
    print('%8.4g,' % ((Mh[s].X-Mhmin)/(Mhmax-Mhmin)),end='')

print('\nCAES charging:',end='')
for s in range(stages):
    print('%8.4g,' % psc_bar[s].X,end='')

print('\nCAES discharg:',end='')
for s in range(stages):
    print('%8.4g,' % psd_bar[s].X,end='')

print('\nCAES heat pwr:',end='')
for s in range(stages):
    print('%8.4g,' % hsh_bar[s].X,end='')

print('\nchp power P:',end='')
for s in range(stages):
    print('%8.4g,' % pchp[s].X, end='')

print('\nchp power Q:',end='')
for s in range(stages):
    print('%8.4g,' % qchp[s].X, end='')

print('\nchp power h:',end='')
for s in range(stages):
    print('%8.4g,' % hchp[s].X, end='')

print('\nhp e power:',end='')
for s in range(stages):
    print('%8.4g,' % ph[s].X, end='')

print('\nhp h power:',end='')
for s in range(stages):
    print('%8.4g,' % hh[s].X, end='')


