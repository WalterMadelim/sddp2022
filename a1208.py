import copy
import numpy as np
import matplotlib.pyplot as plt
from gurobipy import GRB
import gurobipy as gp

def ev2ms(m,v):
    return np.log(m**2/(v+m**2)**.5),np.log(v/m**2+1)**.5
def gsum(a):
    return gp.quicksum(a)

np.random.seed(1)
mu = np.array([67,135,105,80,79,72])
nss = len(mu)
samples = 3000 # a sample is a full rand process
xi = np.zeros((samples,2,nss)) # scene,e/d,t
for sce in range(samples):
    for t in range(nss):
        m,s = ev2ms(1,.25)
        xi[sce,0,t] = np.random.lognormal(m,s)
        m,s = ev2ms(mu[t],((t+1)/20*mu[t])**2)
        xi[sce,1,t] = np.random.lognormal(m,s)


perfInfoObj = np.zeros(samples) # = 17499.933728562126 (rdseed = 1, samples = 3000);
# PerfInfoObj_LPrelax = 5984.709179172252;

# used to generate initial pai
for sm in range(samples):
    m = gp.Model('extensive')
    Y = m.addVars(nss)
    p = m.addVars(nss,ub=897)
    n = m.addVars(nss)
    x = m.addVars(nss)
    D = m.addVars(nss)
    o = m.addVars(nss,ub=201.75)
    y = m.addVars(nss,ub=1,vtype=GRB.BINARY) # vtype=GRB.BINARY

    m.setObjective(gsum(15*p[t]+30*n[t]+6456*y[t]+100*o[t] for t in range(nss)) + 120*n[nss-1])

    m.addConstr(n[0] - p[0] - D[0] + 60 == 0)
    m.addConstr(-Y[0] +.6 * 1 + .4 * xi[sm, 0, 0] == 0)  # random param
    for t in range(1,nss): # 1,2,3,4,5
        m.addConstr(n[t]-p[t]-D[t]+p[t-1]-n[t-1]+x[t-1] == 0)
        m.addConstr(-Y[t] + .6*Y[t-1] + .4*xi[sm,0,t] == 0) # random param

    for t in range(nss):
        m.addConstr(D[t] == .6*mu[t]*Y[t] + .4*xi[sm,1,t]) # random param
    for t in range(nss):
        m.addConstr(22.4 * y[t] + x[t] - o[t] <= 807)
        m.addConstr(538 * y[t] - x[t] >= 0)
        m.addConstr(p[t] + x[t] <= 897)

    m.setParam('OutputFlag',0)
    m.optimize()
    if m.status != 2:
        print('m.status =',m.status)
        exit(3)
    perfInfoObj[sm] = m.ObjVal
    # for t in range(nss):
    #     print('t=%2d:y %8g|x %8g|n %8g|p %8g|Y %8g|D %8g|o %8g' % (t,y[t].X,x[t].X,n[t].X,p[t].X,Y[t].X,D[t].X,o[t].X))
    # print('\n')

beta = np.zeros(96)

for bite in range(300000000):
    # begin conv opt program
    # scenario_time decomposition
    valMat = np.zeros((samples, nss + 1))
    grad = np.zeros_like(beta)
    for sm in range(samples):
        addTerm = 0 # additional term (without x)
        t = 0
        m = gp.Model('extensive')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)
        # next stage base --> use conditional Expectation
        base = np.array([1,xi[sm,0,0],xi[sm,1,0]])
        pai_0,pai_1 = np.dot(base,beta[:3]),np.dot(base,beta[3:6])
        pen = pai_0*(n-p-D) + pai_1*-Y
        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],1,mu[t+1]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[6:11]),np.dot(baseNxt,beta[11:16])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y

        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+pen+penNxt)/samples)

        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag', 0)
        # m.setParam('DualReductions', 0)
        m.optimize()
        if m.status != 2:
            print('m.status =', m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[:3] += base * (n.X-p.X-D.X) / samples
        grad[3:6] += base * (-Y.X) / samples
        grad[6:11] += baseNxt * (p.X - n.X + x.X) / samples
        grad[11:16] += baseNxt * (.6 * Y.X) / samples
        addTerm += (pai_0*60 + pai_1*(.6 + .4*xi[sm,0,t]) + pai_Nxt_1*.4*xi[sm,0,t+1])/samples
        grad[:3] += base*60/ samples
        grad[3:6] += base*(.6 + .4*xi[sm,0,t])/ samples
        grad[11:16] += baseNxt * .4*xi[sm,0,t+1]/ samples

        t = 1
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1]])
        pai_0,pai_1 = np.dot(base,beta[6:11]),np.dot(base,beta[11:16])
        pen = pai_0*(n-p-D) + pai_1*-Y
        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],1,mu[t+1]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[16:23]),np.dot(baseNxt,beta[23:30])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y
        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)

        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t])  # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[6:11] += base * (n.X-p.X-D.X) / samples
        grad[11:16] += base * (-Y.X) / samples
        grad[16:23] += baseNxt * (p.X - n.X + x.X) / samples
        grad[23:30] += baseNxt * (.6 * Y.X) / samples
        addTerm += (pai_1*.4*xi[sm,0,t] + pai_Nxt_1*.4*xi[sm,0,t+1])/samples
        grad[11:16] += base * .4*xi[sm,0,t] / samples
        grad[23:30] += baseNxt * .4*xi[sm,0,t+1] / samples

        t = 2
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2]])
        pai_0,pai_1 = np.dot(base,beta[16:23]),np.dot(base,beta[23:30])
        pen = pai_0*(n-p-D) + pai_1*-Y
        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],1,mu[t+1]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[30:39]),np.dot(baseNxt,beta[39:48])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y

        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t])  # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[16:23] += base * (n.X-p.X-D.X) / samples
        grad[23:30] += base * (-Y.X) / samples
        grad[30:39] += baseNxt * (p.X - n.X + x.X) / samples
        grad[39:48] += baseNxt * (.6 * Y.X) / samples
        addTerm += (pai_1*.4*xi[sm,0,t] + pai_Nxt_1*.4*xi[sm,0,t+1])/samples
        grad[23:30] += base * .4*xi[sm,0,t] / samples
        grad[39:48] += baseNxt * .4*xi[sm,0,t+1] / samples

        t = 3
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1, xi[sm, 0, 0], xi[sm, 1, 0], xi[sm, 0, 1], xi[sm, 1, 1], xi[sm, 0, 2], xi[sm, 1, 2], xi[sm, 0, 3], xi[sm, 1, 3]])
        pai_0, pai_1 = np.dot(base, beta[30:39]), np.dot(base, beta[39:48])
        pen = pai_0 * (n - p - D) + pai_1 * -Y
        baseNxt = np.array([1, xi[sm, 0, 0], xi[sm, 1, 0], xi[sm, 0, 1], xi[sm, 1, 1], xi[sm, 0, 2], xi[sm, 1, 2], xi[sm, 0, 3], xi[sm, 1, 3], 1, mu[t+1]])
        pai_Nxt_0, pai_Nxt_1 = np.dot(baseNxt, beta[48:59]), np.dot(baseNxt, beta[59:70])
        penNxt = pai_Nxt_0 * (p - n + x) + pai_Nxt_1 * .6 * Y

        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t])  # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)
        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[30:39] += base * (n.X - p.X - D.X) / samples
        grad[39:48] += base * (-Y.X) / samples
        grad[48:59] += baseNxt * (p.X - n.X + x.X) / samples
        grad[59:70] += baseNxt * (.6 * Y.X) / samples
        addTerm += (pai_1 * .4 * xi[sm, 0, t] + pai_Nxt_1 * .4 * xi[sm, 0, t + 1]) / samples
        grad[39:48] += base * .4 * xi[sm, 0, t] / samples
        grad[59:70] += baseNxt * .4 * xi[sm, 0, t + 1] / samples


        t = 4
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1, xi[sm, 0, 0], xi[sm, 1, 0], xi[sm, 0, 1], xi[sm, 1, 1], xi[sm, 0, 2], xi[sm, 1, 2], xi[sm, 0, 3],xi[sm, 1, 3],xi[sm, 0, 4],xi[sm, 1, 4]])
        pai_0, pai_1 = np.dot(base, beta[48:59]), np.dot(base, beta[59:70])
        pen = pai_0 * (n - p - D) + pai_1 * -Y
        baseNxt = np.array([1, xi[sm, 0, 0], xi[sm, 1, 0], xi[sm, 0, 1], xi[sm, 1, 1], xi[sm, 0, 2], xi[sm, 1, 2], xi[sm, 0, 3],xi[sm, 1, 3],xi[sm, 0, 4],xi[sm, 1, 4], 1, mu[t+1]])
        pai_Nxt_0, pai_Nxt_1 = np.dot(baseNxt, beta[70:83]), np.dot(baseNxt, beta[83:96])
        penNxt = pai_Nxt_0 * (p - n + x) + pai_Nxt_1 * .6 * Y

        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t])  # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)
        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[48:59] += base * (n.X - p.X - D.X) / samples
        grad[59:70] += base * (-Y.X) / samples
        grad[70:83] += baseNxt * (p.X - n.X + x.X) / samples
        grad[83:96] += baseNxt * (.6 * Y.X) / samples
        addTerm += (pai_1 * .4 * xi[sm, 0, t] + pai_Nxt_1 * .4 * xi[sm, 0, t+1])/samples
        grad[59:70] += base * .4 * xi[sm, 0, t] / samples
        grad[83:96] += baseNxt * .4 * xi[sm, 0, t+1] / samples


        t = 5
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1, xi[sm, 0, 0], xi[sm, 1, 0], xi[sm, 0, 1], xi[sm, 1, 1], xi[sm, 0, 2], xi[sm, 1, 2], xi[sm, 0, 3],xi[sm, 1, 3],xi[sm, 0, 4],xi[sm, 1, 4],xi[sm, 0, 5],xi[sm, 1, 5]])
        pai_0, pai_1 = np.dot(base, beta[70:83]), np.dot(base, beta[83:96])
        pen = pai_0 * (n - p - D) + pai_1 * -Y

        cost = 15*p+150*n+6456*y+100*o
        m.setObjective((cost+pen)/samples)
        m.addConstr(D == .6 * mu[t] * Y + .4 * xi[sm, 1, t])  # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)
        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm, t] = m.ObjVal
        grad[70:83] += base * (n.X - p.X - D.X) / samples
        grad[83:96] += base * (-Y.X) / samples
        addTerm += (pai_1 * .4 * xi[sm, 0, t]) / samples
        grad[83:96] += base * .4 * xi[sm, 0, t] / samples


    val = sum(sum(valMat))
    print(val)
    beta += 1e-12*grad


