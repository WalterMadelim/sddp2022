import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
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
# one_dim convex optimization algorithm
# 20221109
MaxIte = 20
x = -1.26 # initial prox center
v,g = f(x)
n = 0
print('n = %8d, f(%8g)=%8g'%(n,x,v))
cut = np.zeros((MaxIte+1,3)) # (f(x),g,x)
cut[0] = v,g,x
cutIdx = 1
tmin = .05
t = 1
tol_e,tol_g = 1e-6,1e-7
notFirstDesent = 0
for ite in range(MaxIte):
    m = gp.Model('Fk_with_prox_center_x')
    F, y = m.addVar(lb=-GRB.INFINITY), m.addVar(lb=-GRB.INFINITY)
    m.setObjective(F + 1 / (2 * t) * (y - x) ** 2)
    for c in range(cutIdx):
        tmp = cut[c]
        m.addConstr(F >= tmp[0] + tmp[1] * (y - tmp[2]))
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        if m.status == GRB.UNBOUNDED:
            print('opt unbounded >>>>>>>>>>>>')
        else:
            print('opt Fail >>>>>>>>>>>>>')
        exit(3)
    y = y.X
    v_,g_ = f(y)
    gcheck = (x-y)/t
    echeck = v - (F.X+gcheck*(x-y))
    if echeck <= tol_e and abs(gcheck) <= tol_g:
        print('after %d ites, min attained with e_err %g.' % (ite,echeck))
        print('f(%g)=%g is optimal!'%(y,v_))
        break
    taux = t if abs(g_-g) < 1e-4 else t * (1 + (g_ - g) * (y - x) / (g_ - g) ** 2)
    if v-v_ >= .5*(v-m.ObjVal): # armijo test
        if notFirstDesent:
            cut[0] = F.X,gcheck,y
        else: # only the 1st iteration we do this
            notFirstDesent = 1
        cut[1] = v_,g_,y
        cutIdx = 2
        x,v,n = y,v_,n+1 # update current prox_center and best value
        print('n = %8d, f(%8g)=%8g' % (n, x, v))
        t = max(t,min(taux,10*t))
    else: # null step
        cut[cutIdx] = v_,g_,y
        cutIdx += 1
        t = min(t, max(taux,t/2,tmin))
