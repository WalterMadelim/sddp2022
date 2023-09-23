import gurobipy as gp
from gurobipy import GRB
import numpy as np
import sys
from matplotlib import pyplot as plt
ppl = np.array([400,1300,3700,1000,4300,2600,800,5000,3500,2300])
ppd = np.array([250,250,300,250,300,250,250,250,250,250])/1000
ppA = np.pi*(ppd/2)**2
# flow region for generating ppK: 1-1000 t/h
ppK = np.array([416.5714,1353.85705,1547.2653,1041.4285,1798.1732,2707.7141,833.1428,5207.1425,3644.99975,2395.28555])
# need to approx hydro character through PWL
