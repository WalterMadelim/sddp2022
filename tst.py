import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import matplotlib.pyplot as plt

pai = np.array([1,3,4])
print(pai[:3])

# exit(2)
# def cut1(x):
#     return 5.04-4*(x+1.26)
# def cut2(x):
#     return 0.22+3*(x-2.74)
# def quad(x):
#     return -2.96-1/2*(x+1.26)**2
# x = np.array([-2,-1,0,1,2,3,4,5])
# xquad = np.linspace(-2,5,50)
# fig,ax = plt.subplots()
# ax.plot(x,cut1(x))
# ax.plot(x,cut2(x))
# ax.plot(xquad,quad(xquad))
# ax.grid()
# plt.show()