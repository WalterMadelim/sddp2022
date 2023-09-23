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

# dataseq = xi[:,1,5]
# h,v = np.histogram(dataseq,bins=280)
# wid = v[1]-v[0]
# p = v + wid/2
# p = p[:-1]
# fig,ax = plt.subplots()
# ax.bar(p,h,width=wid)
# plt.show()
# print(np.mean(dataseq))

# eval pi_6_2, with dim = 13
A10,b10 = np.zeros((samples, 3)),np.zeros(samples)
A11,b11 = np.zeros((samples, 3)),np.zeros(samples)
A20,b20 = np.zeros((samples, 5)),np.zeros(samples)
A21,b21 = np.zeros((samples, 5)),np.zeros(samples)
A30,b30 = np.zeros((samples, 7)),np.zeros(samples)
A31,b31 = np.zeros((samples, 7)),np.zeros(samples)
A40,b40 = np.zeros((samples, 9)),np.zeros(samples)
A41,b41 = np.zeros((samples, 9)),np.zeros(samples)
A50,b50 = np.zeros((samples, 11)),np.zeros(samples)
A51,b51 = np.zeros((samples, 11)),np.zeros(samples)
A60,b60 = np.zeros((samples, 13)),np.zeros(samples)
A61,b61 = np.zeros((samples, 13)),np.zeros(samples)
for sm in range(samples):
    A60[sm, 0] = 1
    A60[sm, 1] = xi[sm, 0, 0]
    A60[sm, 2] = xi[sm, 1, 0]
    A60[sm, 3] = xi[sm, 0, 1]
    A60[sm, 4] = xi[sm, 1, 1]
    A60[sm, 5] = xi[sm, 0, 2]
    A60[sm, 6] = xi[sm, 1, 2]
    A60[sm, 7] = xi[sm, 0, 3]
    A60[sm, 8] = xi[sm, 1, 3]
    A60[sm, 9] = xi[sm, 0, 4]
    A60[sm, 10] = xi[sm, 1, 4]
    A60[sm, 11] = xi[sm, 0, 5]
    A60[sm, 12] = xi[sm, 1, 5]
A61 = copy.deepcopy(A60)
for sm in range(samples):
    A50[sm, 0] = 1
    A50[sm, 1] = xi[sm, 0, 0]
    A50[sm, 2] = xi[sm, 1, 0]
    A50[sm, 3] = xi[sm, 0, 1]
    A50[sm, 4] = xi[sm, 1, 1]
    A50[sm, 5] = xi[sm, 0, 2]
    A50[sm, 6] = xi[sm, 1, 2]
    A50[sm, 7] = xi[sm, 0, 3]
    A50[sm, 8] = xi[sm, 1, 3]
    A50[sm, 9] = xi[sm, 0, 4]
    A50[sm, 10] = xi[sm, 1, 4]
A51 = copy.deepcopy(A50)
for sm in range(samples):
    A40[sm, 0] = 1
    A40[sm, 1] = xi[sm, 0, 0]
    A40[sm, 2] = xi[sm, 1, 0]
    A40[sm, 3] = xi[sm, 0, 1]
    A40[sm, 4] = xi[sm, 1, 1]
    A40[sm, 5] = xi[sm, 0, 2]
    A40[sm, 6] = xi[sm, 1, 2]
    A40[sm, 7] = xi[sm, 0, 3]
    A40[sm, 8] = xi[sm, 1, 3]
A41 = copy.deepcopy(A40)
for sm in range(samples):
    A30[sm, 0] = 1
    A30[sm, 1] = xi[sm, 0, 0]
    A30[sm, 2] = xi[sm, 1, 0]
    A30[sm, 3] = xi[sm, 0, 1]
    A30[sm, 4] = xi[sm, 1, 1]
    A30[sm, 5] = xi[sm, 0, 2]
    A30[sm, 6] = xi[sm, 1, 2]
A31 = copy.deepcopy(A30)
for sm in range(samples):
    A20[sm, 0] = 1
    A20[sm, 1] = xi[sm, 0, 0]
    A20[sm, 2] = xi[sm, 1, 0]
    A20[sm, 3] = xi[sm, 0, 1]
    A20[sm, 4] = xi[sm, 1, 1]
A21 = copy.deepcopy(A20)
for sm in range(samples):
    A10[sm, 0] = 1
    A10[sm, 1] = xi[sm, 0, 0]
    A10[sm, 2] = xi[sm, 1, 0]
A11 = copy.deepcopy(A10)

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
    y = m.addVars(nss,ub=1) # vtype=GRB.BINARY

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
    g = m.getConstrs()
    b10[sm] = -g[0].Pi
    b11[sm] = -g[1].Pi
    b20[sm] = -g[2].Pi
    b21[sm] = -g[3].Pi
    b30[sm] = -g[4].Pi
    b31[sm] = -g[5].Pi
    b40[sm] = -g[6].Pi
    b41[sm] = -g[7].Pi
    b50[sm] = -g[8].Pi
    b51[sm] = -g[9].Pi
    b60[sm] = -g[10].Pi
    b61[sm] = -g[11].Pi

r = np.linalg.lstsq(A10, b10, rcond=None)[0] # initial vector for pi_6_2
r = np.concatenate((r,np.linalg.lstsq(A11, b11, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A20, b20, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A21, b21, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A30, b30, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A31, b31, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A40, b40, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A41, b41, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A50, b50, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A51, b51, rcond=None)[0]))
r = np.concatenate((r,np.linalg.lstsq(A60, b60, rcond=None)[0]))
beta = np.concatenate((r,np.linalg.lstsq(A61, b61, rcond=None)[0])) # dim=96
beta = np.zeros(96)
print(beta[:3])
print(beta[3:6])
print(beta[6:11])
print(beta[11:16])
for bite in range(1):
    # begin conv opt program
    # scenario_time decomposition
    valMat = np.zeros((samples, nss + 1))
    grad = np.zeros_like(beta)
    for sm in range(samples):
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
        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],1,mu[1]])
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
        print(m.ObjVal)

        grad[:3] += base * (n.X-p.X-D.X) / samples
        grad[3:6] += base * (-Y.X) / samples
        grad[6:11] += baseNxt * (p.X - n.X + x.X) / samples
        grad[11:16] += baseNxt * (.6 * Y.X) / samples

        print(grad[:16])
        exit(12)

        addTerm = 0 # additional term (without x)
        t = 1
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],2,mu[2]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[10:17]),np.dot(baseNxt,beta[17:24])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y
        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1]])
        pai_0,pai_1 = np.dot(base,beta[:5]),np.dot(base,beta[5:10])
        pen = pai_0*(n-p-D) + pai_1*-Y
        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        addTerm += pai_1 * .4 * xi[sm,0,t]
        m.addConstr(D == mu[t]/5 * Y + .8 * xi[sm,1,t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[10:17] += baseNxt * (p.X - n.X + x.X) / samples
        grad[17:24] += baseNxt * .6 * Y.X / samples
        grad[:5] += base * (n.X-p.X-D.X) / samples
        grad[5:10] += base * -Y.X / samples

        grad[5:10] += base * .4 * xi[sm, 0, t] / samples


        t = 2
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],3,mu[3]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[24:33]),np.dot(baseNxt,beta[33:42])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y
        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2]])
        pai_0,pai_1 = np.dot(base,beta[10:17]),np.dot(base,beta[17:24])
        pen = pai_0*(n-p-D) + pai_1*-Y
        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        addTerm += pai_1 * .4 * xi[sm, 0, t]
        m.addConstr(D == mu[t]/5 * Y + .8 * xi[sm,1,t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[24:33] += baseNxt * (p.X - n.X + x.X) / samples
        grad[33:42] += baseNxt * .6 * Y.X / samples
        grad[10:17] += base * (n.X-p.X-D.X) / samples
        grad[17:24] += base * -Y.X / samples

        grad[17:24] += base * .4 * xi[sm, 0, t] / samples

        t = 3
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],xi[sm,0,3],xi[sm,1,3],4,mu[4]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[42:53]),np.dot(baseNxt,beta[53:64])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y
        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],xi[sm,0,3],xi[sm,1,3]])
        pai_0,pai_1 = np.dot(base,beta[24:33]),np.dot(base,beta[33:42])
        pen = pai_0*(n-p-D) + pai_1*-Y
        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        addTerm += pai_1 * .4 * xi[sm, 0, t]
        m.addConstr(D == mu[t]/5 * Y + .8 * xi[sm,1,t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[42:53] += baseNxt * (p.X - n.X + x.X) / samples
        grad[53:64] += baseNxt * .6 * Y.X / samples
        grad[24:33] += base * (n.X-p.X-D.X) / samples
        grad[33:42] += base * -Y.X / samples

        grad[33:42] += base * .4 * xi[sm, 0, t] / samples

        t = 4
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        baseNxt = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],xi[sm,0,3],xi[sm,1,3],xi[sm,0,4],xi[sm,1,4],5,mu[5]])
        pai_Nxt_0,pai_Nxt_1 = np.dot(baseNxt,beta[64:77]),np.dot(baseNxt,beta[77:90])
        penNxt = pai_Nxt_0*(p-n+x) + pai_Nxt_1*.6*Y
        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],xi[sm,0,3],xi[sm,1,3],xi[sm,0,4],xi[sm,1,4]])
        pai_0,pai_1 = np.dot(base,beta[42:53]),np.dot(base,beta[53:64])
        pen = pai_0*(n-p-D) + pai_1*-Y
        cost = 15*p+30*n+6456*y+100*o
        m.setObjective((cost+penNxt+pen)/samples)
        addTerm += pai_1 * .4 * xi[sm, 0, t]
        m.addConstr(D == mu[t]/5 * Y + .8 * xi[sm,1,t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[64:77] += baseNxt * (p.X - n.X + x.X) / samples
        grad[77:90] += baseNxt * .6 * Y.X / samples
        grad[42:53] += base * (n.X-p.X-D.X) / samples
        grad[53:64] += base * -Y.X / samples

        grad[53:64] += base * .4 * xi[sm, 0, t] / samples

        t = 5
        m = gp.Model('scena')
        Y = m.addVar()
        p = m.addVar(ub=897)
        n = m.addVar()
        x = m.addVar()
        D = m.addVar()
        o = m.addVar(ub=201.75)
        y = m.addVar(vtype=GRB.BINARY)

        base = np.array([1,xi[sm,0,0],xi[sm,1,0],xi[sm,0,1],xi[sm,1,1],xi[sm,0,2],xi[sm,1,2],xi[sm,0,3],xi[sm,1,3],xi[sm,0,4],xi[sm,1,4],xi[sm,0,5],xi[sm,1,5]])
        pai_0,pai_1 = np.dot(base,beta[64:77]),np.dot(base,beta[77:90])
        pen = pai_0*(n-p-D) + pai_1*-Y
        cost = 15*p+150*n+6456*y+100*o
        m.setObjective((cost+pen)/samples)
        addTerm += pai_1 * .4 * xi[sm, 0, t]
        m.addConstr(D == mu[t]/5 * Y + .8 * xi[sm,1,t]) # random param
        m.addConstr(22.4 * y + x - o <= 807)
        m.addConstr(538 * y - x >= 0)
        m.addConstr(p + x <= 897)

        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != 2:
            print('m.status =',m.status)
            exit(3)
        valMat[sm,t] = m.ObjVal
        grad[64:77] += base * (n.X-p.X-D.X) / samples
        grad[77:90] += base * -Y.X / samples

        grad[77:90] += base * .4 * xi[sm, 0, t] / samples

        valMat[sm,t+1] = addTerm/samples

    val = sum(sum(valMat))
    print(val)
    beta += 1e-14*grad


