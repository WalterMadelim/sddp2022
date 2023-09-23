import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-2,2,400)
print(x)
y = abs(x)**(3/2)
print(y)
dy = 3/2*(abs(x))**(1/2)
fig,ax = plt.subplots()
ax.plot(x,dy)
plt.show()