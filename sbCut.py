import gurobipy as gp
from gurobipy import GRB
import sys

def getPai_in_Sb(x1,x2):
    m = gp.Model('getPaiTrial(LP)')
    y, z1, z2 = m.addVar(), m.addVar(), m.addVar()
    m.addConstr(z1 == x1)
    m.addConstr(z2 == x2)
    m.addConstr(y + z1 / 4 + z2 / 2 >= 2.6)
    m.addConstr(y <= 4)
    m.addConstr(z1 <= 1)
    m.addConstr(z2 <= 1)
    m.setObjective(4 * y)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        sys.exit(3)
    l = m.getConstrs()  # Y,p,n,x,const
    pai1 = l[0].Pi
    pai2 = l[1].Pi
    return pai1,pai2
def cutSb(x1,x2):
    pai1,pai2 = getPai_in_Sb(x1,x2)
    m = gp.Model('getSbCutCoeff(IP)')
    y,z1,z2 = m.addVar(vtype=GRB.INTEGER),m.addVar(vtype=GRB.INTEGER),m.addVar(vtype=GRB.INTEGER)
    m.addConstr(y + z1 / 4 + z2 / 2 >= 2.6)
    m.addConstr(y <= 4)
    m.addConstr(z1 <= 1)
    m.addConstr(z2 <= 1)
    m.setObjective(4 * y - pai1 * z1 - pai2 * z2)
    m.setParam('OutputFlag', 0)
    m.optimize()
    if m.status != GRB.OPTIMAL:
        print('opt Fail >>>>>>>>>>>>>')
        sys.exit(3)
    return m.ObjVal,pai1,pai2
x1,x2 = 0,0 # last stage action
cst,cx1,cx2 = cutSb(x1,x2) # get sb cut
print(cst,cx1,cx2)
