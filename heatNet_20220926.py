import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
from matplotlib import pyplot as plt

# CHP at node 0 and heat pump with TES at node 2
hppl = 1000*np.array([3.5,1.75,1.75,.75,1.75])
hppmdot = np.array([49,25,15,10,24])
hdf = np.exp(-0.12*hppl/4200/hppmdot) # heat dissipation factor
Tam = 20
Ts0 = 100
Cm = 4200
hLoad = np.array([0,0,0,15*30*Cm,24*30*Cm,10*30*Cm])

# Ts1 = Tam + hdf[0] * (Ts0 - Tam)
# Ts4 = Tam + hdf[4] * (Ts1 - Tam)
# Ts2 = Tam + hdf[1] * (Ts1 - Tam)
# Ts3 = Tam + hdf[2] * (Ts2 - Tam)
# Ts5 = Tam + hdf[3] * (Ts2 - Tam)
Tr3 = Ts3 - hLoad[3] / Cm / 15
Tr4 = Ts4 - hLoad[4] / Cm / 24
Tr5 = Ts5 - hLoad[5] / Cm / 10

Tr2 = (15 * (Tam + hdf[2] * (Tr3 - Tam)) + 10 * (Tam + hdf[3] * (Tr5 - Tam)))/25
Tr1 = (25 * (Tam + hdf[1] * (Tr2 - Tam)) + 24 * (Tam + hdf[4] * (Tr4 - Tam)))/49
Tr0 = Tam + hdf[0] * (Tr1 - Tam)

hGen = Cm * 49 * (Ts0 - Tr0)
print('heat gened by chp')
print(hGen)
print('heat load total')
print( sum(hLoad) )
print(Tr0,Tr1,Tr2)


