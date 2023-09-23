import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np

# generate Markov Chain's transition matrices: 23 Ps
# 20221004

T = 24
S = 240
ls = 3 # levels
xi = np.load('myxi.npy') # [level,t,feaN]
d = np.load('xisamples.npy') # (feaN, T, S)


P = np.zeros((T-1,ls,ls),dtype=np.float64)
for stage in range(T-1):
    cls = np.zeros((2,S)) # t=0 -> t=1
    for s in range(S):
        for t in range(2):
            sp = d[:,stage+t,s]
            dist = np.zeros(ls,dtype=np.float64)
            for l in range(ls):
                dist[l] = np.linalg.norm(sp-xi[l, stage+t, :])
            cls[t,s] = np.argmin(dist)
    tmp = np.zeros(ls,dtype=np.float64)
    for c in range(ls):
        a = cls[1][np.where(cls[0] == c)]
        tmp[0] = len(np.where(a==0)[0]) / len(a)
        tmp[1] = len(np.where(a==1)[0]) / len(a)
        tmp[2] = len(np.where(a==2)[0]) / len(a)
        P[stage,:,c] = tmp

for stage in range(T-1):
    print(P[stage])

np.save('MCPs.npy',P)





# T = 24
# clsN = 3
# for t in range(T):
#     print('t=',t)
#     for level in range(clsN):
#         print(xi[level,t,:])
