import gurobipy as gp
import matplotlib.pyplot as plt
from gurobipy import GRB
import numpy as np
import sys
import xlrd
import xlwt
T = 24
S = 240
a = np.zeros((T,S),dtype=np.float64)
wb = xlwt.Workbook()
ws = wb.add_sheet('elds')

sh = xlrd.open_workbook('ele.xls').sheet_by_name("Sheet1")
for t in range(T):
    for s in range(S):
        ws.write(1+t,1+s,sh.cell_value(s,t))

wb.save('elestd.xls')
