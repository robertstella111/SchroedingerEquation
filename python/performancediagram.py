import time
from Schroedinger import Schroedinger
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import math



xkoord = [1,2,4,6,8,12,14,16]
y = [10.3,5.216, 3.294, 2.820,2.894,3.586,3.8,5.639]

plt.rc('grid', linestyle=':', color='g', linewidth=0.5)
fig, ax = plt.subplots()
ax.plot(xkoord,y, "r-o")

plt.grid()
plt.xlabel("Anzahl OpenMP Threads")
plt.ylabel("Zeit in s")
plt.title("Performance für 1008 Stützpunkte")


plt.show()


