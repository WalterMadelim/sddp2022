import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
# 20221109
# variable metric bundle method
# 1-dim min convex f(x)

D1 = 4
def q(pai):
    m = gp.Model('extensive')
    y0 = m.addVar(vtype=GRB.BINARY)
    x0 = m.addVar()
    p1,n1 = m.addVar(),m.addVar()
    m.setObjective(y0+x0+2*p1+4*n1+pai*(p1 - n1 - x0 + D1))
    m.addConstr(x0<=7*y0)
    m.setParam('OutputFlag',0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        if m.status == GRB.UNBOUNDED:
            return -9e9,999
        else:
            print('opt Fail >>>>>>>>>>>>>')
        exit(3)
    return m.ObjVal, p1.X-n1.X-x0.X+D1
def f(pai):
    tmp = q(pai)
    return -tmp[0],-tmp[1]
m1,m2 = .5,.8
x = -1.26 # initial prox center
k,n = 0,0
mu = 1 # mu(0) = 1
MaxIte = 10
cut = np.zeros((MaxIte+1,3)) # (f(x),g,x)
v,g = f(x)
cut[k] = v,g,x
xL = 999*np.ones(MaxIte+1)
vL = copy.deepcopy(xL)
xL[n],vL[n] = x,v
for ite in range(MaxIte):
    t = 1 # initially, t = 1
    tL, tR = 0, 9e9
    while 1:  # CS, input a t
        m = gp.Model('Fk_with_prox_center_x')
        F,y = m.addVar(lb = -GRB.INFINITY),m.addVar(lb = -GRB.INFINITY)
        m.setObjective(F + mu/(2*t)*(y-x)**2)
        for c in range(k+1):
            tmp = cut[c]
            m.addConstr(F >= tmp[0] + tmp[1]*(y-tmp[2]))
        m.setParam('OutputFlag',0)
        m.optimize()
        if m.status != GRB.OPTIMAL:
            if m.status == GRB.UNBOUNDED:
                print('opt unbounded >>>>>>>>>>>>')
            else:
                print('opt Fail >>>>>>>>>>>>>')
            exit(3)
        v_,g_ = f(y.X)
        if v_ > 8e9:
            print('y gened is not valid, please increase prox term!')
            exit(7)
        y = y.X
        Gcheck = mu/t*(x-y)
        deltat = v-m.ObjVal
        echeck = v-(F.X+Gcheck*(x-y))
        if abs(Gcheck) <= 1e-6 and echeck <= 1e-6:
            print('f(x) Minimum attained after %d iterates.'%ite)
            print('the proximal centers are:')
            print(xL[:n + 1])
            print(vL[:n + 1])
            print('cuts are')
            print(cut[:k + 1])
            for i in cut[:k + 1]:
                print(i[0] - i[1] * i[2])
            exit(98)
        if v-v_ >= m1 * deltat: # step 2, if (31) satisfied
            tL = t
            if g_*(y-x) >= -m2*deltat: # step 4, if (34) satisfied
                csFlag = 1 # descent
                break
            hold36 = abs(Gcheck) <= 1e-6 or Gcheck*(y-x) >= -echeck # m4 = 1
            if tR == 9e9 and hold36: # step 5
                csFlag = 2 # cutting plane
                break
        else: # step 3, armijo is NOT satisfied
            tR = t
            et = v - (v_ + g_*(x-y))
            if et <= deltat and tL == 0: # (33) satisfied, m3 = 1
                csFlag = 0 # null
                break
        # step 6
        if tR == 9e9:
            t *= 1.5
        else:
            t = (tL+tR)/2
    if csFlag: # if not Null
        xnext,Gn = y,Gcheck
        if csFlag == 1:
            xi = xnext - x # 4
            v = g_ - g # 7
            mu = 1/(1/mu + v*xi/v**2) # .6
        x,v = y,v_ # update current best point and best value
        n += 1
        xL[n], vL[n] = x, v # record
    k += 1
    cut[k] = v_,g_,y


