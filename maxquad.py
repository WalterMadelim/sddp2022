import copy
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
def f(x):
    f = np.zeros(5)
    for o in range(5):
        a = A[o]
        bb = b[o]
        f[o] = sum(x[i]*a[i,j]*x[j] for i in range(10) for j in range(10)) + sum(bb[i]*x[i] for i in range(10))
    return f[0]
A = np.zeros((5,10,10))
for j in range(1,5+1):
    for i in range(1,10+1):
        for k in range(1,10+1):
            if k > i:
                A[j-1,i-1,k-1] = np.exp(i/k)*np.cos(i*k)*np.sin(j)
for j in range(1,5+1):
    for i in range(1,10+1):
        for k in range(1,10+1):
            if k < i:
                A[j-1,i-1,k-1] = A[j-1,k-1,i-1]
for j in range(1,5+1):
    for i in range(1,10+1):
        tmp = A[j-1,i-1]
        tmp = sum(abs(tmp))
        A[j-1,i-1,i-1] = abs(np.sin(j))*i/10 + tmp
b = np.zeros((5,10))
for j in range(1,5+1):
    for i in range(1,10+1):
        b[j-1,i-1] = np.exp(i/j)*np.sin(i*j)
print(minimize(f, np.ones(10)))