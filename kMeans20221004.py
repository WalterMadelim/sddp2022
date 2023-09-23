import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from matplotlib import pyplot as plt
import xlrd

# k-means clustering: generate xi for my sddp
# 20221004
# temperature /= 10
# elec load (total): 1.6-5.4MW
# heat load (total): 3-5MW
# water load (total): 3.5-5.5 hundred_t/h

T = 24
Ss,S = 80,240
clsN = 3 # classfied into 3 classes
feaN = 4 # dimension of features
xi = np.zeros((3,T,feaN)) # used directly in sddp
bk = xlrd.open_workbook('xi.xls')
sh = bk.sheet_by_name("tems")
tems = np.zeros((T, S), dtype=np.float64)
for t in range(T):
    for s in range(S):
        tems[t,s] = sh.cell_value(t+1,s+1)
sh = bk.sheet_by_name("elds")
elds = np.zeros((T, S), dtype=np.float64)
for t in range(T):
    for s in range(S):
        elds[t,s] = sh.cell_value(t+1,s+1)
sh = bk.sheet_by_name("hlds")
hlds = np.zeros((T, S), dtype=np.float64)
for t in range(T):
    for s in range(S):
        hlds[t,s] = sh.cell_value(t+1,s+1)
sh = bk.sheet_by_name("wlds")
wlds = np.zeros((T, S), dtype=np.float64)
for t in range(T):
    for s in range(S):
        wlds[t,s] = sh.cell_value(t+1,s+1)
d = np.zeros((feaN, T, S),dtype=np.float64)
d[0],d[1],d[2],d[3] = tems,elds,hlds,wlds
xik = np.zeros((clsN,feaN),dtype=np.float64) # 3 centers, every vector 4 components

for t in range(T):
    for dim in range(feaN):
        tmp = np.sort(d[dim, t, :])
        for c in range(clsN):
            xik[c,dim] = np.average(tmp[c*Ss:(c+1)*Ss])
    for g in range(10):
        totDist = 0
        dist = np.zeros(clsN,dtype=np.float64) # a sample point's dist to 3 ref points
        cls = np.zeros(S) # record each sample's class Num
        for s in range(S):
            sp = d[:, t, s] # sample point
            for c in range(clsN):
                dist[c] = np.linalg.norm(sp - xik[c,:])
            cls[s] = np.argmin(dist)
            totDist += np.min(dist)
        N0,N1,N2 = 0,0,0
        sm0 = np.zeros(4,dtype=np.float64)
        sm1 = np.zeros(4,dtype=np.float64)
        sm2 = np.zeros(4,dtype=np.float64)
        for s in range(S):
            if cls[s] == 0:
                N0 += 1
                sm0 += d[:,t,s]
            elif cls[s] == 1:
                N1 += 1
                sm1 += d[:,t,s]
            elif cls[s] == 2:
                N2 += 1
                sm2 += d[:,t,s]
        # print('%4d|%8g'%(g,totDist))
        # print(xik)
        # print(N0, N1, N2)
        xik[0,:] = sm0/N0
        xik[1,:] = sm1/N1
        xik[2,:] = sm2/N2

    for level in range(clsN):
        xi[level,t,:] = xik[level,:]
print(xi)
for t in range(T):
    print('t=',t)
    for level in range(clsN):
        print(xi[level,t,:])

np.save('myxi.npy',xi)