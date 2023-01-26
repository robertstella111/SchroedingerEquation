from Schroedinger import Schroedinger
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import math


length = 1002
b = 4
a = 0
obj = Schroedinger.Schroedinger1DParallel()
obj.setNumPoints(length)
obj.setDomain(a,b)
obj.setNumThreads(8)
obj.solve(5, 1e-4)
n = 3

eigenvalues = obj.getEigenValues()

eigenvector = obj.getEigenvektor(2)
error = []
xkoord = []

error1 = []
xkoord2 = obj.getXKoord()  

#error for all eigenvalues 
for i in range(length-2):
    fehler = (eigenvalues[i]-(((i+1)*math.pi)/(b-a))**2)/(((i+1)*math.pi/(b-a)))**2
    error.append(100*fehler)
    xkoord.append(i+1)
#error for eigenvector, sign can vary because sign of eigenvector can have -sin and +sin
for i in range(len(eigenvector)):
    fehler = eigenvector[i] + np.sqrt(2/(b-a))*np.sin((n*math.pi*(xkoord2[i]-a))/(b-a))
    error1.append(fehler)


plt.rc('grid', linestyle=':', color='g', linewidth=0.5)
fig, ax = plt.subplots()
ax.plot(xkoord2,error1)
plt.grid()
plt.xlabel("x in m")
plt.ylabel("Amplitude")
plt.title("Lösung der Schrödinger Gleichung")

plt.show()